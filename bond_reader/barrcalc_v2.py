import json
import numpy as np
import bz2
import sys
import os

# barrier = 0, bond in A = 0.4, bond in B = 0.4 by default
BARRIER = 0
shiftA = 0
shiftB = 0
KA = 5
KB = 5

def getbonds(fname):
    if os.path.isfile(fname):
        fp = open(fname, 'rb')
    elif os.path.isfile(fname+'.bz2'):
        fp = bz2.open(fname+'.bz2', mode='rb')
    else:
        raise Exception("bond file not exist.")
    while True:
        line = fp.readline()
        if line == b'':
            break
        if line.startswith(b'ITEM: TIMESTEP'):
            step = int(fp.readline().decode())
        if line.startswith(b'ITEM: NUMBER OF ENTRIES'):
            number = int(fp.readline().decode())
        if line.startswith(b'ITEM: ENTRIES'):
            bonds = np.genfromtxt(fp, max_rows=number, dtype=int)
            yield {'step': step, 'bonds': bonds}
    fp.close()

def getatoms(fname):
    if os.path.isfile(fname):
        fp = open(fname, 'rb')
    elif os.path.isfile(fname+'.bz2'):
        fp = bz2.open(fname+'.bz2', mode='rb')
    else:
        raise Exception("atom file not exist.")
    while True:
        line = fp.readline()
        if line == b'':
            break
        if line.startswith(b'ITEM: TIMESTEP'):
            step = int(fp.readline().decode())
        if line.startswith(b'ITEM: NUMBER OF ATOMS'):
            number = int(fp.readline().decode())
        if line.startswith(b'ITEM: BOX BOUNDS'):
            bound = np.genfromtxt(fp, max_rows=3)
        if line.startswith(b'ITEM: ATOMS'): # ITEM: ATOMS id type xs ys zs
            atoms = np.genfromtxt(fp, max_rows=number)
            yield {'step': step, 'bound': bound, 'atoms':atoms}
    fp.close()

def lj(r, eps, sigma, shift):
    r -= shift
    temp = (sigma / r)**6
    if r/sigma > 2**(1/6):
        return 0
    else:
        return 4*eps*(temp*(temp - 1)) + eps

def harmonic(r, K):
    return K * r * r

def wrap(p):
    return (p + np.array([0.5, 0.5, 0.5])) % np.array([1,1,1]) - np.array([0.5, 0.5, 0.5])

def dist(p1, p2, bound):
    # calculate distance, p1 and p2 are scaled from 0 to 1.
    r = wrap(p1 - p2)
    scale = bound[:, 1] - bound[:, 0]
    r *= scale
    return np.sqrt(np.sum(r * r))
    
def get_theta(ps):
    p1 = ps[0]
    p2 = ps[1]
    p3 = ps[2]
    v21 = p1 - p2
    v23 = p3 - p2
    angle = np.arccos((np.dot(v21,v23)/(np.linalg.norm(v21)*np.linalg.norm(v23))))
    return angle

def costheta(p1, p2):
    vs = []
    for p in (p1, p2):
        v = (p + np.array([0.5, 0.5, 0.5])) % np.array([1,1,1]) - np.array([0.5, 0.5, 0.5])
        v /= np.sqrt(np.sum(p * p))
        vs.append(v)
    #print(vs)
    return np.sum(vs[0] * vs[1])

def make_fixed_map(fixed):
    m = {}
    for bond in fixed:
        pair = tuple(bond)
        m[pair[0]] = pair[1]
        m[pair[1]] = pair[0]
    return m

def calc(name):
    bonds_gen = getbonds(name+'_reac-bond')
    atoms_gen = getatoms(name+'_reac.lammpstrj')

    fixed_bonds = None
    atype = None
    typeshift = {(2,2):shiftA, (2,3):(shiftA+shiftB)/2, (3,2):(shiftA+shiftB)/2, (3,3):shiftB}
    record = []
    try:
        old_bonds = set()
        while True:
            bonds = next(bonds_gen)
            atoms = next(atoms_gen)
            if fixed_bonds is None:
                fixed_bonds = set([frozenset([i[1], i[2]]) for i in bonds['bonds'] if i[0] != 1])
                fixed_map = make_fixed_map(fixed_bonds)
            if atype is None:
                atype = {int(i[0]): int(i[1]) for i in atoms['atoms']}

            all_bonds = set([frozenset([i[1], i[2]]) for i in bonds['bonds'] if i[0] == 1])
            new_bonds = all_bonds - old_bonds
            #print(new_bonds)
            all_pos = atoms['atoms']
            bound = atoms['bound']
            step = bonds['step']
            id2pos = {int(i[0]): i[2:5] for i in all_pos}
            for bond_x in new_bonds:
                bond = tuple(bond_x)
                #print(bond)
                ids = [fixed_map[bond[0]], bond[0], bond[1], fixed_map[bond[1]]]
                
                if fixed_map[bond[0]] > bond[0]:                                    
                    mono1_ids = [fixed_map[bond[0]]+1,fixed_map[bond[0]],bond[0]]  
                else:                                                               
                    mono1_ids = [fixed_map[bond[0]]-1,fixed_map[bond[0]],bond[0]]   
                mono1_pos = [id2pos[i] for i in mono1_ids]                          
                if fixed_map[bond[1]] > bond[1]:                                    
                    mono2_ids = [fixed_map[bond[1]]+1,fixed_map[bond[1]],bond[1]]   
                else:                                                               
                    mono2_ids = [fixed_map[bond[1]]-1,fixed_map[bond[1]],bond[1]]   
                mono2_pos = [id2pos[i] for i in mono1_ids]                          
                if fixed_map[bond[1]] > bond[1]:                                    
                    mono2_ids = [fixed_map[bond[1]]+1,fixed_map[bond[1]],bond[1]]   
                else:                                                               
                    mono2_ids = [fixed_map[bond[1]]-1,fixed_map[bond[1]],bond[1]]   
                mono2_pos = [id2pos[i] for i in mono1_ids]                          
                if atype[ids[0]] == 2:                                              
                    e_ang1 = harmonic(get_theta(mono1_pos)-np.pi,KA)                
                else:                                                               
                    e_ang1 = harmonic(get_theta(mono1_pos)-np.pi,KB)                
                if atype[ids[3]] == 2:                                              
                    e_ang2 = harmonic(get_theta(mono2_pos)-np.pi,KA)                
                else:                                                               
                    e_ang2 = harmonic(get_theta(mono2_pos)-np.pi,KB)                
                
                
                
                
                atom_pos = [id2pos[i] for i in ids]
                ps = [wrap(atom_pos[i+1] - atom_pos[i]) for i in range(3)] # scaled, from 0 to 1
                d = [
                    dist(atom_pos[0], atom_pos[1], bound),
                    dist(atom_pos[0], atom_pos[3], bound),
                    dist(atom_pos[2], atom_pos[3], bound)
                ]
                shift = [
                    typeshift[(atype[ids[0]], atype[ids[0]])], 
                    typeshift[(atype[ids[3]], atype[ids[3]])], 
                    typeshift[(atype[ids[0]], atype[ids[3]])], 
                ]
                e_lj = lj(d[1], 4, 1, shift[2])
                e_bond = harmonic(d[0]-0.4-shift[0]/2, 2000) + harmonic(d[2]-0.4-shift[1]/2, 2000)
                e_barr = BARRIER
                e_ang = e_ang1 + e_ang2                         
                e_total = e_lj + e_bond + e_barr + e_ang        
                cos = costheta(ps[0], ps[2])
                #print(d, [atype[ids[0]], atype[ids[3]]], e_bond, e_lj)
                record.append({
                    'step': step,
                    'atoms': [i.tolist() for i in atom_pos],
                    'bound': bound.tolist(),
                    'cos': cos,
                    'type': [atype[ids[0]], atype[ids[3]]],
                    'energy': {
                        'lj': e_lj,
                        'bond': e_bond,
                        'barr': e_barr,
                        'ang': e_ang,       
                        'total': e_total
                    }
                })
                old_bonds = all_bonds
    except StopIteration:
        pass
    eAAt = np.array([i['energy']['total'] for i in record if i['type']==[2,2]])
    eABt = np.array([i['energy']['total'] for i in record if i['type']==[2,3] or i['type']==[3,2]])
    eBBt = np.array([i['energy']['total'] for i in record if i['type']==[3,3]])
    print('''total barrier energy between : \nAA %g (%d x)\nAB %g (%d x)\nBB %g (%d x)\n''' % (eAAt.mean(), len(eAAt), eABt.mean(), len(eABt), eBBt.mean(), len(eBBt)))
    with open(name+'_reacbarr_v2.json', 'w') as fp:
        json.dump(record, fp, indent=2)
    sys.exit()
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('input filename, or -barr as barrier height, -shiftA and -shiftB for potential shift. (enter arguments before the file name.)')
        exit()
    i = 1
    while i < len(sys.argv):
        if sys.argv[i]=='-barr':
            BARRIER = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-shiftA':
            shiftA = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-shiftB':
            shiftB = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-KA':
            KA = int(sys.argv[i+1])
            i += 2
            continue
        if sys.argv[i]=='-KB':
            KB = int(sys.argv[i+1])
            i += 2
            continue
        calc(sys.argv[i])
        i += 1
