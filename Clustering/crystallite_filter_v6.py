import chainanalysis
import pandas as pd
import numpy as np
import os
import subprocess as sb
import sys

def calc(filename,clusterstep,Nfilter,csvfile,nskip):
    #Filenames and properties for I/O
    #filename = '210412_copoly3_standard-4'
    #Nfilter = 12 #minimum number of monomers in cluster for splitting
    inputfile = open(filename+'.data')
    inputlines = inputfile.readlines()
    numheaderlines = 33 #DONT TOUCH THIS: 33 for standard 21 for early output
    typefile = open(filename+'-type')
    typelines = typefile.readlines()
    typefile.close()
    numatoms = int(str.split(inputlines[2])[0])
    numbonds = int(str.split(inputlines[4])[0])
    numangles = int(str.split(inputlines[6])[0])
    numskips = 3



    #df = pd.read_csv('./'+filename+'-all_clusters.csv')
    df = pd.read_csv(csvfile)
    cluster = df['cluster'].values


    #Splitting out crystallite data files based on clustering analysis
    values, counts = np.unique(cluster, return_counts=True)
    ctr=0
    print('{} clusters generated.'.format(len(values)))
    if nskip > 0:
        for value in values[:nskip]:
            print('Skippping cluster {}.'.format(value))
    for i,value in enumerate(values[nskip:]):
        string = str(i)
        df_temp = df[df.cluster==value]
        #df_temp.to_csv('cluster'+string+'_coords.csv')
        ids = df_temp['atomid'].values
        id_strings = [str(id) for id in ids]



        dfile = open(filename+'.lammpstrj')
        dout = open('{}snap_cluster'.format(clusterstep)+string+'.lammpstrj','w')
        vfile = open(filename+'.lammpsvel')
        vout = open('{}snap_cluster'.format(clusterstep)+string+'.lammpsvel','w')
        bfile = open(filename+'-bond')
        bout = open('{}snap_cluster'.format(clusterstep)+string+'-bond','w')
        outputfile = open('{}snap_cluster'.format(clusterstep)+string+'.data','w')
        typeout = open('{}snap_cluster'.format(clusterstep)+string+'-type','w')
        clusteratoms = 0
        clusterbonds = 0
        clusterangles = 0
        entrynum = 0
        entries = []
        print('Splitting cluster {}.'.format(i))
        
        
        #THIS BLOCK SHOULD SKIP SMALL CLUSTERS - ADJUST SIZE BY NFILTER
        if len(id_strings) < (3*Nfilter):
            print('Cluster less than {} monomers, ignoring.'.format(Nfilter))
            typeout.close()
            bfile.close()
            bout.close()
            dfile.close()
            dout.close()
            vfile.close()
            vout.close()
            inputfile.close()
            outputfile.close()
            
            #Cleanup
            os.remove('{}snap_cluster'.format(clusterstep)+string+'.data')
            os.remove('{}snap_cluster'.format(clusterstep)+string+'.lammpstrj')
            os.remove('{}snap_cluster'.format(clusterstep)+string+'.lammpsvel')
            os.remove('{}snap_cluster'.format(clusterstep)+string+'-bond')
            os.remove('{}snap_cluster'.format(clusterstep)+string+'-type')
            ctr+=1

        
        else:
                #JSON file
            print('Copying JSON file.')
            sb.call(["cp",filename+'_param.json','{}snap_cluster'.format(clusterstep)+string+'_param.json'])
        
        
                #Data file
            print('Splitting data file.')
        
        
            #Header
            for line in inputlines[:numheaderlines]:
                outputfile.write(line)
                
            # #Making sure all atoms in every chain are included in the cluster
            # for line in inputlines[2*numatoms+numheaderlines+2*numskips:2*numatoms+numheaderlines+2*numskips+numbonds]:
                # split = str.split(line)
                # if split[1] != '1':
                    # if split[2] in id_strings and split[3] not in id_strings:
                        # id_strings.append(split[3])
                        # missing_id = np.array([int(split[3])])
                        # ids = np.append(ids,missing_id)
                    # elif split[2] not in id_strings and split[3] in id_strings:
                        # id_strings.append(split[2])
                        # missing_id = np.array([int(split[2])])
                        # ids = np.append(ids,missing_id)


            # #Gotta run through this a second time in case we had an extra side atom rather than were short one (first loop would catch center, second loop then catches opposite side, all monomers should then be complete)
            # for line in inputlines[2*numatoms+numheaderlines+2*numskips:2*numatoms+numheaderlines+2*numskips+numbonds]:
                # split = str.split(line)
                # if split[1] != '1':
                    # if split[2] in id_strings and split[3] not in id_strings:
                        # id_strings.append(split[3])
                        # missing_id = np.array([int(split[3])])
                        # ids = np.append(ids,missing_id)
                    # elif split[2] not in id_strings and split[3] in id_strings:
                        # id_strings.append(split[2])
                        # missing_id = np.array([int(split[2])])
                        # ids = np.append(ids,missing_id)                    


                    
                
            #Atoms
            for line in inputlines[numheaderlines:numheaderlines+numatoms]:
                split = str.split(line)
                if split[0] in id_strings:
                    outputfile.write(split[0]+' '+split[2]+' '+split[2]+' '+split[3]+' '+split[4]+' '+split[5]+'\n')
                    clusteratoms += 1
            #Velocities
            for line in inputlines[numheaderlines+numatoms:numheaderlines+numatoms+numskips]:
                outputfile.write(line)
            for line in inputlines[numheaderlines+numatoms+numskips:2*numatoms+numheaderlines+numskips]:
                split = str.split(line)
                if split[0] in id_strings:
                    outputfile.write(line)
                    
                    
            #Bonds
            for line in inputlines[2*numatoms+numheaderlines+numskips:2*numatoms+numheaderlines+2*numskips]:
                outputfile.write(line)
            for line in inputlines[2*numatoms+numheaderlines+2*numskips:2*numatoms+numheaderlines+2*numskips+numbonds]:
                split = str.split(line)
                if split[2] in id_strings and split[3] in id_strings:
                    outputfile.write(line)
                    clusterbonds += 1

            #Angles
            for line in inputlines[2*numatoms+numheaderlines+2*numskips+numbonds:2*numatoms+numheaderlines+3*numskips+numbonds]:
                outputfile.write(line)
            for line in inputlines[2*numatoms+numheaderlines+3*numskips+numbonds:]:
                split = str.split(line)
                if len(split) < 3:
                    pass
                elif split[2] in id_strings and split[3] in id_strings and split[4] in id_strings:
                    outputfile.write(line)
                    clusterangles += 1
            
        
            inputfile.close()
            outputfile.close()
        
        
        
            #Correcting header
            if clusteratoms!= len(ids):
                print('This is bad! Length of IDs vector: {} Clusteratoms:{} on Cluster {}'.format(len(ids),clusteratoms,value))
            sb.call(["sed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ atoms:'+"{}".format(clusteratoms)+'\ atoms:','{}snap_cluster'.format(clusterstep)+string+'.data'])
            sb.call(["sed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ bonds:'+"{}".format(clusterbonds)+'\ bonds:','{}snap_cluster'.format(clusterstep)+string+'.data'])
            sb.call(["sed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ angles:'+"{}".format(clusterangles)+'\ angles:','{}snap_cluster'.format(clusterstep)+string+'.data'])

        
        
        
        
            #Typefile
            print('Splitting type file.')
        
            for line in typelines[:3]:
                typeout.write(line)
            for line in typelines[3:4]:
                typeout.write(str(clusteratoms)+'\n')
            for line in typelines[4:9]:
                typeout.write(line)
            for line in typelines[9:]:
                split = str.split(line)
                if split[0] in id_strings:
                    typeout.write(line)
            typeout.close()
        
        
            #Bond dump
            print('Splitting bond dump.')
            for i,line in enumerate(bfile):
                split = str.split(line)
                if line.startswith('ITEM: TIMESTEP'):
                    bout.write(line)
                    if i > 0:
                        entries.append(entrynum) 
                        entrynum = 0
     
                
                elif line.startswith('ITEM: NUMBER OF ENTRIES'):
                    bout.write(line)
                    bout.write('x\n')
                    next(bfile,None)
                    continue
                elif len(split) == 3:
                    if split[1] in id_strings and split[2] in id_strings:
                        bout.write(line)
                        entrynum+=1
                else:
                    bout.write(line)
        
            entries.append(entrynum)
            bfile.close()
            bout.close()

        
            print('Correcting entry numbers.')
            for entry in entries:
                sb.call(["sed","-i",'1,/x/{s/x/'+"{}".format(entry)+'/}','{}snap_cluster'.format(clusterstep)+string+'-bond'])
        
        
        
        
        
            #Trajectory file
            print('Splitting trajectory file.')
        
            for line in dfile:
                split = str.split(line)
                if len(split) == 5:
                    if split[0] in id_strings:
                        dout.write(line)    
                elif line.startswith('ITEM: NUMBER OF ATOMS'):
                    dout.write(line)
                    dout.write('x\n')
                    next(dfile,None)
                    continue
                else:
                    dout.write(line)
            
            dfile.close()
            dout.close()
            
            #Correcting trajectory entry numbers
            sb.call(["sed","-i",'/^x/ s/x*$/'+"{}".format(clusteratoms)+'/g','{}snap_cluster'.format(clusterstep)+string+'.lammpstrj'])
            
            
            
            #Velocity file
            print('Splitting velocities file.')
        
            for line in vfile:
                split = str.split(line)
                if len(split) == 5:
                    if split[0] in id_strings:
                        vout.write(line)    
                elif line.startswith('ITEM: NUMBER OF ATOMS'):
                    vout.write(line)
                    vout.write('x\n')
                    next(vfile,None)
                    continue
                else:
                    vout.write(line)
            
            vfile.close()
            vout.close()
            
            #Correcting trajectory entry numbers
            sb.call(["sed","-i",'/^x/ s/x*$/'+"{}".format(clusteratoms)+'/g','{}snap_cluster'.format(clusterstep)+string+'.lammpsvel'])
        
            ctr+=1
    




if __name__ == '__main__':
    if len(sys.argv) < 8:
        print('Usage: crystallite_clustering_v6 [FLAGS] filename. No extensions on filename, and filename of trajectory and bond file must match.\n')
        print('Use _#snap_ filename for single timestep extraction or simulation filename for extraction from full simulation.\n')
        print('Necessary flags: -step timestep of clustering snapshot, -filter smallest cluster size extracted (in # of monos), -csv filename for cluster csv.\n')
        print('Optional flag: -skip number of clusters to skip before extracting (can be used for restarting an interrupted extraction).\n')
        exit()
    i = 1
    nskip = 0
    while i < len(sys.argv):
        if sys.argv[i]=='-step':
            step = int(sys.argv[i+1])
            i += 2
        if sys.argv[i]=='-filter':
            filter = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-csv':
            csvfile = sys.argv[i+1]
            i += 2
            continue
        if sys.argv[i]=='-skip':
            nskip = int(sys.argv[i+1])
            i+=2
            continue
        calc(sys.argv[i],step,filter,csvfile,nskip)
        i += 1
