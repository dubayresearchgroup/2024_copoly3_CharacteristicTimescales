import chainanalysis
import pandas as pd
import numpy as np
import os
import subprocess as sb
import sys
import io
#from sklearn.cluster import DBSCAN
#from sklearn.cluster import OPTICS
import hdbscan
from sklearn.neighbors import NearestNeighbors
import copy

#Cantor function for bijective mapping of cluster ids when merging clusters
def cantor(i,j):
    c = int(0.5*(i+j)*(i+j+1) + j)
    return c

#Functions for extracting specific timesteps of traj and bond files - ZZ wrote these
def fetch(name, step, dest_name):
    with open(dest_name, 'w') as fp:
        fp.write(bi_fetch(name, step))

#%%
def bi_fetch(name, step, counts=1):
    with open(name, 'rb') as fp:
        #mm = mmap.mmap(fp.fileno(), 0)
        def get_ps(c):
            i = c
            step = None
            fp.seek(i)
            while True:
                line = fp.readline()
                if line.startswith(b'ITEM: TIMESTEP'):
                    break
                if line == b'':
                    break
                i = fp.tell()
            if line != b'':
                secondln = fp.readline()
                step = int(secondln)
            return i, step
        size = fp.seek(0, io.SEEK_END)
        left = 0
        left_step = get_ps(0)[1]
        left2, left2_step = get_ps(1)
        right, right_step = get_ps(size)
        pos = None
        while True:
            m = (left + right) // 2
            m, m_step = get_ps(m)
            #print(left_step, right_step)
            if step <= left_step:
                pos = left 
                break
            if step <= left2_step:
                pos = left2
                break
            if step == right_step:
                pos = right
                break
            if step <= m_step:
                #right, right_step = get_ps(m)
                right = m
                right_step = m_step
                left = left2
                left_step = left2_step
            else:
                #left, left_step = get_ps(m)
                left = m
                left_step = m_step
            left2, left2_step = get_ps(left+1)
                
            if right <= left2 and left2 < step:
                break
        if pos is None:
            return ''
        else:
            fp.seek(pos)
            txt = fp.readline()
            n = counts
            while n > 0:
                line = fp.readline()
                if line == b'':
                    break
                if line.startswith(b'ITEM: TIMESTEP'):
                    n -= 1
                txt += line
            return txt.decode()
            
    
#Clustering analysis    
def cluster_alg(filename,step,cutoff,criteria,S_weight,box_size=50.0):
    reader = chainanalysis.BondReader(filename)
    reader.calc()
    allseqs = reader.frames[-1]['numseq']
    numheaderlines = 9
    inputfile = open(filename+'.lammpstrj','r')
    tempfile = open('temp.csv','w')
    inputlines = inputfile.readlines()
    try:
        tempfile.write('atomid,type,x,y,z\n')
        for line in inputlines[numheaderlines:]:
            split = str.split(line)
            tempfile.write(split[0]+','+split[1]+','+str(box_size*float(split[2]))+','+str(box_size*float(split[3]))+','+str(box_size*float(split[4]))+'\n')
    finally:
        inputfile.close()
        tempfile.close()


    print('Generating crystallite clusters.')
    df = pd.read_csv('temp.csv')
    df = df.sort_values(by=['atomid'])
    df.to_csv('temp.csv', index=False)
    X = df.values
    X_centers = X[1::3,:]

    #Filtering by neighbor density before clustering
    #criteria = 12 #number of monos to be considered aggregated
    neighbors_model = NearestNeighbors(radius = cutoff)
    neighbors_model.fit(X_centers[:,2:])
    neighborhoods = neighbors_model.radius_neighbors(X_centers[:,2:],return_distance=False)
    X_filter = []
    filter_centers = []
    filter_ids = []
    num_neighbors = []
    orient_vecs = []
    S_order = []

    #Calculating Slocal as a feature for DBSCAN
    for i,neighbor_list in enumerate(neighborhoods):
        num_neighbors.append(len(neighbor_list))
        if len(neighbor_list) > criteria:
            X_filter.append(X_centers[i,:].tolist())
            filter_centers.append(int(X_centers[i,0]))
            filter_ids.append(int(X_centers[i,0]))
            filter_ids.append(int(X_centers[i,0])-1)
            filter_ids.append(int(X_centers[i,0])+1)
            vec = np.array([X[(3*i+1),2] - X[(3*i+1)+1,2], X[(3*i+1),3] - X[(3*i+1)+1,3], X[(3*i+1),4] - X[(3*i+1)+1,4]])
            for neighbor in neighbor_list:
                orient_vecs.append([X[(3*neighbor+1),2] - X[(3*neighbor+1)+1,2], X[(3*neighbor+1),3] - X[(3*neighbor+1)+1,3], X[(3*neighbor+1),4] - X[(3*neighbor+1)+1,4]])
            orient_vecs = np.array(orient_vecs)
            n_vec = np.array([np.mean(orient_vecs[:,0]),np.mean(orient_vecs[:,1]),np.mean(orient_vecs[:,2])])
            cos_angle = (vec[0]*n_vec[0] + vec[1]*n_vec[1] + vec[2]*n_vec[2])/(((vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5)*((n_vec[0]**2 + n_vec[1]**2 + n_vec[2]**2)**0.5))
            S = (3*(cos_angle**2)-1)/2
            S_order.append(S)
            orient_vecs = []


    #weighting the order parameter and adding it to the feature space
    #S_weight = 1.
    if S_weight != 0.: 
        X_filter = np.hstack((X_filter,S_weight*(np.array(S_order).reshape(-1,1))))
    else:
        X_filter = np.array(X_filter)


    #HDBSCAN algorithm
    model = hdbscan.HDBSCAN(min_cluster_size = criteria,min_samples=2*criteria,allow_single_cluster=True)
    model.fit(X_filter[:,2:])
    cluster = model.labels_
    
    # S_order_all = np.zeros((len(X[:,0]),), dtype=int)
    cluster_all = -1*np.ones((len(X[:,0]),), dtype=int)
    cluster_merge = np.copy(cluster_all)
    cluster_all_unmerge = np.copy(cluster_all)
    for i,value in enumerate(cluster):
        cluster_all[int(X_filter[i,0] - 1)] = value
        cluster_all[int(X_filter[i,0])] = value
        cluster_all[int(X_filter[i,0] - 2)] = value
        # S_order_all[int(X_filter[i,0] - 1)] = S_order[i]
        # S_order_all[int(X_filter[i,0])] = S_order[i]
        # S_order_all[int(X_filter[i,0] - 2)] = S_order[i]



    #Forcing all monomers in a given chain to belong to a single cluster
    cluster_all += 1 #IMPORTANT - ASSIGNS NEGATIVE VALUES TO ZERO
    cluster_seq = copy.deepcopy(allseqs)
    cluster_all_unmerge = np.copy(cluster_all)
    first_choices = []
    second_choices = []
    dovetail = []

    for i,seq in enumerate(allseqs):
        multi_cluster = True
        for j,monoid in enumerate(seq):
            cluster_seq[i][j] = cluster_all[monoid-1]
        cluster_val, votes = np.unique(cluster_seq[i],return_counts=True)
        top_choice = cluster_val[np.argmax(votes)]
        if top_choice == 0:
            try:
                top_choice = cluster_val[[votes.tolist().index(k) for k in sorted(votes.tolist(), reverse=True)][1]]
            except IndexError:
                pass
        for j,monoid in enumerate(seq):
            if cluster_all[monoid-1] == 0:
                cluster_all[monoid-1] = top_choice
                cluster_seq[i][j] = top_choice
        cluster_val, votes = np.unique(cluster_seq[i],return_counts=True)
        top_choice = cluster_val[np.argmax(votes)]
        first_choices.append(top_choice)
        try:
            second_choice = cluster_val[[votes.tolist().index(k) for k in sorted(votes.tolist(), reverse=True)][1]]
            second_choices.append(second_choice)
        except IndexError:
            multi_cluster = False
            second_choices.append(top_choice)    
        if multi_cluster:
            dovetail.append(sorted([top_choice,second_choice]))
        else:
            dovetail.append([top_choice,top_choice])

        for monoid in seq:
            cluster_all[monoid-1] = top_choice
            
            
    for pair in dovetail:
        pair.sort()
    unique_pairs = sorted(set(sorted([tuple(sorted(x)) for x in dovetail])))

    for pair in unique_pairs:
        if pair[0] != pair[1]:        
            newid = cantor(pair[0],pair[1])
            for i,value in enumerate(cluster_all):
                if ((value==pair[0]) or (value==pair[1])) and value>0:
                    cluster_all[i] = newid
                    
                    
                    
    #Communicating number and size of clusters generated
    values, counts = np.unique(cluster_all_unmerge, return_counts=True)
    ctr=0
    print('Unmerged!')
    print('{} clusters generated.'.format(len(values)))
    for i,value in enumerate(values):
        print('Cluster {}: {} monomers'.format(ctr,int(counts[i]/3)))
        ctr+=1        
            
    values, counts = np.unique(cluster_all, return_counts=True)
    ctr=0
    print('Merged!')
    print('{} clusters generated.'.format(len(values)))
    for i,value in enumerate(values):
        print('Cluster {}: {} monomers'.format(ctr,int(counts[i]/3)))
        ctr+=1
        
    #Updating dataframe and saving results
    df['cluster'] = cluster_all_unmerge
    df.to_csv('./'+filename+'_clusteringv6-unmerge_clusters.csv', index=False)
    df['cluster'] = cluster_all
    df.to_csv('./'+filename+'_clusteringv6-all_clusters.csv', index=False)
    os.remove('temp.csv')
    
    
    
    
# For running from command line
def calc(simname, step, cutoff, neigh, S_weight):
    fetch(simname+'.lammpstrj',step,simname+'_{}snap.lammpstrj'.format(step))
    fetch(simname+'.lammpsvel',step,simname+'_{}snap.lammpsvel'.format(step))
    fetch(simname+'-bond',step,simname+'_{}snap-bond'.format(step))
    #### These will break if the snap is the last step in a trajectory just comment
    sb.call('head -n -1 '+simname+'_{}snap.lammpstrj'.format(step)+' > ./tmp && mv ./tmp '+simname+'_{}snap.lammpstrj'.format(step),shell=True)
    sb.call('head -n -1 '+simname+'_{}snap.lammpsvel'.format(step)+' > ./tmp && mv ./tmp '+simname+'_{}snap.lammpsvel'.format(step),shell=True)
    sb.call('head -n -1 '+simname+'_{}snap-bond'.format(step)+' > ./tmp && mv ./tmp '+simname+'_{}snap-bond'.format(step),shell=True)
    #### this block out and run for that last snap only
    sb.call('cp '+simname+'-type ./'+simname+'_{}snap-type'.format(step),shell=True)
    sb.call('cp '+simname+'_param.json ./'+simname+'_{}snap_param.json'.format(step),shell=True)
    cluster_alg(simname+'_{}snap'.format(step),step,cutoff, neigh, S_weight)
    
    
    
    
if __name__ == '__main__':
    if len(sys.argv) < 10:
        print('Usage: crystallite_clustering_v6 [FLAGS] filename. No extensions on filename, and filename of trajectory and bond file must match. \n')
        print('Necessary flags: -step step to take snapshot, -cut cutoff distance for DBSCAN, -neigh minimum neighbor density for initial clusters, -Sweight weight for Slocal feature in DBSCAN. \n')
        exit()
    i = 1
    while i < len(sys.argv):
        if sys.argv[i]=='-step':
            step = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-cut':
            cutoff = float(sys.argv[i+1])
            i += 2
            continue    
        if sys.argv[i]=='-neigh':
            neigh = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-Sweight':
            S_weight = float(sys.argv[i+1])
            i += 2
            continue
        #print(sys.argv)
        calc(sys.argv[i],step,cutoff,neigh,S_weight)
        i += 1
