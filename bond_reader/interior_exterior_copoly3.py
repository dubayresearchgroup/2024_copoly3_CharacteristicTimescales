import chainanalysis
import pandas as pd
import numpy as np
import os
import subprocess as sb
import sys
import io
import hdbscan
from sklearn.neighbors import NearestNeighbors
import copy
import pickle


    
#Clustering analysis    
def calculate_int_ext(filename,criteria_int,criteria_ext,box_size=50.0):
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


    df = pd.read_csv('temp.csv')
    df = df.sort_values(by=['atomid'])
    df.to_csv('temp.csv', index=False)
    X = df.values
    X_centers = X[1::3,:]
    
    # edge_to_center_map = {}
    # for center_id in X_centers[:,0]:
        # edge_to_center_map[int(center_id-1)] = int(center_id)
        # edge_to_center_map[int(center_id+1)] = int(center_id)
        

    #Filtering by neighbor density before clustering
    #criteria = 12 #number of monos to be considered aggregated
    neighbors_model = NearestNeighbors(radius = 2.5)
    neighbors_model.fit(X_centers[:,2:])
    neighborhoods = neighbors_model.radius_neighbors(X_centers[:,2:],return_distance=False)
    X_int = []
    int_centers = []
    int_ids = []
    X_ext = []
    ext_centers = []
    ext_ids = []
    X_both = []
    both_centers = []
    both_ids = []
    num_neighbors = []
    for i,neighbor_list in enumerate(neighborhoods):
        num_neighbors.append(len(neighbor_list))
        if len(neighbor_list) > criteria_int:
            X_int.append(X_centers[i,:].tolist())
            int_centers.append(int(X_centers[i,0]))
            X_both.append(X_centers[i,:].tolist())
            both_centers.append(int(X_centers[i,0]))
        elif len(neighbor_list) > criteria_ext and len(neighbor_list) <= criteria_int: 
            X_ext.append(X_centers[i,:].tolist())
            ext_centers.append(int(X_centers[i,0]))
            X_both.append(X_centers[i,:].tolist())
            both_centers.append(int(X_centers[i,0]))
    




    ends_list = [entry for pair in reader.frames[-1]['chainends'] for entry in pair]
    ends = len(ends_list)
    end_neigh = []
    ends_list_no_mono = []
    end_neigh_no_mono = []
    for end in ends_list:
        # center_id = edge_to_center_map[end]
        center_id = end
        center_index = int((center_id-2)/3)
        end_neigh.append(num_neighbors[center_index])
        for seq in allseqs:
            if (end in seq) and (len(seq)>3):
                ends_list_no_mono.append(end)
    for end in ends_list_no_mono:
        center_id = end
        center_index = int((center_id-2)/3)
        end_neigh_no_mono.append(num_neighbors[center_index])
    
        
    
    


    if X_int != []:
        X_int = np.array(X_int)
        mono_int = len(X_int[:,1])
        int_As = 0
        int_Bs = 0
        int_ends = 0
        for i,center in enumerate(X_int[:,1]):
            if int(center)==2:
                int_As += 1
            else:
                int_Bs += 1
        
            if int(X_int[i,0]) in ends_list:
                int_ends += 1
    
        FA_int = int_As/mono_int
        FB_int = int_Bs/mono_int
        Fends_int = int_ends/ends
    else:
        FA_int = float("nan")
        FB_int = float("nan")
        Fends_int = float("nan")
        mono_int = float("nan")
    
    
    
    if X_ext != []:
        X_ext = np.array(X_ext)
        mono_ext = len(X_ext[:,1])
        ext_As = 0
        ext_Bs = 0
        ext_ends = 0
        for i,center in enumerate(X_ext[:,1]):
            if int(center)==2:
                ext_As += 1
            else:
                ext_Bs += 1
            if int(X_ext[i,0]) in ends_list:
                ext_ends += 1
        FA_ext = ext_As/mono_ext
        FB_ext = ext_Bs/mono_ext
        Fends_ext = ext_ends/ends
    else:
        FA_ext = float("nan")
        FB_ext = float("nan")
        Fends_ext = float("nan")
        mono_ext = float("nan")
    
    
    
    
    if X_both != []:
        X_both = np.array(X_both)
        mono_both = len(X_both[:,1])
        both_As = 0
        both_Bs = 0
        both_ends = 0
        for i,center in enumerate(X_both[:,1]):
            if int(center)==2:
                both_As += 1
            else:
                both_Bs += 1
            if int(X_both[i,0]) in ends_list:
                both_ends += 1
        FA_both = both_As/mono_both
        FB_both = both_Bs/mono_both
        Fends_both = both_ends/ends

    else:
        FA_both = float("nan")
        FB_both = float("nan")
        Fends_both = float("nan")
        mono_both = float("nan")
    




    print('Internal FA {}'.format(FA_int))
    print('Internal FB {}'.format(FB_int))
    print('External FA {}'.format(FA_ext))
    print('External FB {}'.format(FB_ext))
    print('Full structure FA {}'.format(FA_both))
    print('Full structure FB {}'.format(FB_both))


    print('Internal ends {}'.format(Fends_int))
    print('External ends {}'.format(Fends_ext))
    print('Full structure ends {}'.format(Fends_both))
    print('Average end neighbors {}'.format(np.mean(np.array(end_neigh))))

    # monos_ratio = [int_As,int_Bs,mono_int,ext_As,ext_Bs,mono_ext,both_As,both_Bs,mono_both]

    # with open('./jar/'+filename+'_int-ext.pickle','wb') as f:
        # pickle.dump(monos_ratio,f)
        
    # monos_ratio = [int_As,int_Bs,mono_int,ext_As,ext_Bs,mono_ext,both_As,both_Bs,mono_both]

    dataset = [ends_list,end_neigh,Fends_int,Fends_ext,Fends_both,mono_int,mono_ext,num_neighbors,ends_list_no_mono,end_neigh_no_mono]

    with open('./'+filename+'_end_int_ext.pickle','wb') as f:
        pickle.dump(dataset,f)

    os.remove('temp.csv')





if __name__ == "__main__":
    usage = '''Usage: interior_exterior filename cutoff_int cutoff_ext
No extension for filename. Cutoffs are optional, must be integer of monos within 2.5sigma and int > ext.
'''
    if len(sys.argv) < 2:
       print(usage)
    elif len(sys.argv) < 3:
       filename = sys.argv[1]
       criteria_int = int(60)
       criteria_ext = int(12)
       calculate_int_ext(filename,criteria_int,criteria_ext)
    elif len(sys.argv) < 4:
       filename = sys.argv[1]
       criteria_int = int(sys.argv[2])
       criteria_ext = int(12)
       calculate_int_ext(filename,criteria_int,criteria_ext)
    else: 
       filename = sys.argv[1]
       criteria_int = int(sys.argv[2])
       criteria_ext = int(sys.argv[3])
       calculate_int_ext(filename,criteria_int,criteria_ext)
