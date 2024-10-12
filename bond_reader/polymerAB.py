import re
from collections import Counter
from collections import deque
import numpy as np


# default configuration:
# DFT_poly_conf = {
#     "atom_names": {1:'', 2:'A', 3:'B'},
#     "newbond_type": 1
# }

class Polymer:
    """A class that encapsulates a set of polymers from a simulation run at any given simulation snapshot.  This method takes in raw bond data as well as atom type mappings and naming rules for atom types and bond types 
    to produce a polymer class that can calculate several values of interest including pAA/pBB/pAB, degree of polymerization, PDI, etc.  This class is used by chainanalysis.py for analysing and storing simulation results. 
    
    run .setN()
    (change settings)
    run .calc()
    output using .getrecord() and .gettraj()

    """
    def __init__(self, raw, atomtypes, conf):
        self.raw = raw
        self.tmap = atomtypes
        self.conf = conf
        self.chains = []
        self.seqs = []

    @classmethod
    # atomtypes: {id:type}; conf need to have 'atom_names': {atomtype: atomnames}
    def countAB(cls, atomtypes, conf):
        """Takes mapping of atom ID to atom type along with polymer configuration dictionary (conf) and creats a collections.Counter object from the unique naming of each monomer.

        Parameters
        ----------
        atomtypes : dict of int
            A dictionary mapping the atom ID (keys) to it's corresponding atom type (values)
        conf : dict
            A configuration dictionary for the polymer class containing an element with a dictionary mapping of atom names (i.e. 'A','B', etc.) to the corresponding atom types (i.e. 1,2, etc.) called 'atom_names'

        Returns
        -------
        Counter
            A collections.Counter object which keeps track of the number of monomers of each type in an efficient manner
        """
        atom2name = {i:conf["atom_names"][atomtypes[i]] for i in atomtypes}
        return Counter(atom2name.values())

    @classmethod
    # atomtypes: {id:type}; conf need to have 'atom_names': {atomtype: atomnames}
    def setN(cls, atomtypes, conf):
        """Sets values for the number of monomers of type A and B in a simulation using the values in an atomtypes dictionaries as well as a configuration dictionary.  
           This function sets these values as class variables instead of as instance variables and therefore must be accessed with class methods.

        Parameters
        ----------
        atomtypes : {atomID : atomType}
            A dictionary containing a mapping from every atom ID in the simulation and its corresponding atom type.
        conf : {'atom_names' : {atomType : atomalias},'newbond_type' : bondtype_of_new_bond}
            A dictionary containing an entry 'atom_names' that contains a dictionary with mapping between the integer atom types of the simulation and their corresponding readable alias (i.e. 'A', 'B', etc.).
            The dictionary can also contain an entry 'newbond_type' which gives the bond type of new bonds formed in the simulation.  
        """
        atomnames = [conf["atom_names"][j] for j in atomtypes.values()]
        counter = Counter(atomnames)
        cls.numA = counter['A']
        cls.numB = counter['B']
        cls.N = cls.numA + cls.numB

    @classmethod
    def getN(cls):
        """Returns number of monomers of A and B in a simulation previously set as class variables using PolymerAB.setN

        Returns
        -------
        {'numA' : numA, 'numB' : numB, 'Nmono' : N}
            A dictionary containing the number of monomers of type 'A', type 'B', and the sum of these two numbers.
        """
        return {
            "numA": cls.numA,
            "numB": cls.numB,
            "Nmono": cls.N
        }
    
    def buildchain(self):
        """Creates a list of polymer chains in a simulation.  Each chain in the list is another list of atom IDs of each atom in the chain 
        arranged in bonding order from chain terminal to chain tereminal.  
        The list of chains is then assigned to the instance of the PolymerAB object as self.chains
        """
        self.chains = []

        #bonding neighbors
        self.neighs = {i:[] for i in self.tmap}
        for i in range(self.raw.shape[0]):
            atom1, atom2 = self.raw[i][1:3]
            self.neighs[atom1].append(atom2)
            self.neighs[atom2].append(atom1)

        #self.neighs = {atomID : [atomID_1,atomID_2,atomID_3,...]} 

        #self.mark = set({bonded_atomID1,bonded_atomID2,...})
        #chain = [bonded_atomID1,bonded_atomID2,...]
        #self.chains = [[chain1],[chain2],[chain3]]
        self.mark = set()
        for i in self.tmap:
            if i in self.mark:
                continue
            neigh_i = self.neighs[i]
            self.mark.add(i)
            if len(neigh_i) == 0:
                self.chains.append([i])
            else:
                right = self.neighs[i][0]
                self.mark.add(right)
                chain = self.growchain(i, right)
                if len(neigh_i) == 2:
                    left = self.neighs[i][1]
                    if left not in self.mark:
                        self.mark.add(left)
                        chainleft = self.growchain(i, left)
                        chainleft.reverse()
                        chain = chainleft[:-1] + chain
                self.chains.append(chain)

    def growchain(self, a, b):
        """Creates an ordered list of all bonded atoms in a polymer chain from atom 'a' to the terminal polymer monomer given by the direction of monomer 'a'->monomer 'b'

        Parameters
        ----------
        a : int
            The atom ID of the monomer to begin growing the chain with.
        b : int
            The atom ID of a second monomer connected to monomer 'a', the chain grown by growchain will be grown along the direction of the bond from a->b

        Returns
        -------
        ls : list of int
            A list of bonded atom IDs listed in order from given monomer 'a' to a terminal chain monomer, along the path that includes monomer 'b'
        """
        # grow the chain in the direction from a to b.
        ls = [a, b]
        while True:
            if len(self.neighs[ls[-1]]) < 2:
                break
            last1, last2 = self.neighs[ls[-1]]
            if last1 in self.mark and last2 in self.mark:
                break
            extend = last1 if last2 in self.mark else last2
            self.mark.add(extend)
            ls.append(extend)
        return ls

    def literalchain(self):
        """Converts each chain in self.chains to a string of letters (i.e. [1,2,3,4,5,...]->'ABBBA...') the result is then saved as an instance variable as self.seqs
        """
        #build self.seqs (chains with A B letters) from self.chains (chains by atom id)
        self.seqs = []
        for chain in self.chains:
            trans = [self.conf["atom_names"][self.tmap[i]] for i in chain]
            seq = ''.join(trans)
            #if the chain is a ring, move the cutting point for the ring so that ..AAA|BBB... (else homopolymer chains)
            if (len(self.neighs[chain[0]]) == 2):  
                m = re.search(seq[0] + '+$', seq)   #'A+$' or 'B+$'
                if (m):
                    seq = m.group() + seq[:m.start()]
                seq += 'c'   # 'c' for cycle
            self.seqs.append(seq)

    def chainlength(self):
        """Calculates the length of all chains as well as the length of all continuous 'A' and 'B' monomer sequences and returns the result as three Counter objects with the count of all length values.
        These values are started as object instance variables as self.chaindist, self.bloackAdist, and self.blockBdist respectively.
        """
        chains = []
        blockAs = []
        blockBs = []
        for chain in self.seqs:
            l = len(chain)
            if chain[-1] == 'c':
                chains.append(l-1)
            else:
                chains.append(l)
            blockAs += list(map(len, re.findall('A+', chain)))
            blockBs += list(map(len, re.findall('B+', chain)))
        self.chaindist = list(zip(*sorted(list(Counter(chains).items()))))
        self.blockAdist = list(zip(*sorted(list(Counter(blockAs).items()))))
        self.blockBdist = list(zip(*sorted(list(Counter(blockBs).items()))))

    def makeseq(self):
        """Creates a list of chain sequences from self.seqs that have a chain length>1.
        The resulting list is saved under the instance variable self.seqlist which is sorted from longest sequence to shortest.  
        The number of monomers that are excludued from the sequence list are saved as instance variables self.monoAcount and self.monoBcount.
        """
        counter = Counter(self.seqs)
        self.monoAcount = counter['A']
        self.monoBcount = counter['B']
        self.seqlist = [seq for seq in self.seqs if len(seq) > 1]
        self.seqlist.sort(key=len, reverse=True)


    def reactions(self):
        """Calculates relevant information on the current state of the polymerization reaction and saves them as instance variables, this includes:
           
           -Number average of polymer chains
           -Mass average of polymer chains
           -Polydispersity Index (Number Average/Mass Average)
           -Extent of Reaction (#Bonds/#Monomers)
        """
        # reaction kinetics. 
        chaindop = list(map(len, self.seqs))
        nchain = len(chaindop)
        sum1 = sum(chaindop)
        sum2 = sum(map(lambda x: x*x, chaindop))
        self.numave = sum1 / nchain
        self.massave = sum2 / sum1
        self.pdi = self.massave / self.numave
        bonds = sum1 - nchain
        self.extent = bonds / self.N

    def countpairs(self):
        """Uses the string chain sequences generated by literalchain to generate statistics on the number of AA, AB, and BB pairs. 
        """
        countAA = 0
        countBB = 0
        countAB = 0
        for i in self.seqs:
            if (len(i) >= 2):
                countAA += len(re.findall('(?=AA)', i))
                countBB += len(re.findall('(?=BB)', i))
                countAB += len(re.findall('(?=AB|BA)', i))
                if (i[-1]=='c'):
                    if (i[0] == 'A' and i[-2] == 'A'):
                        countAA += 1
                    elif (i[0] == 'B' and i[-2] == 'B'):
                        countBB += 1
                    else:
                        countAB += 1
        countAll = countAA + countBB + countAB
        if countAll == 0:
            self.pAA = self.pBB = self.pAB = self.pAABB = None
        else:
            self.pAA = countAA / countAll
            self.pBB = countBB / countAll
            self.pAB = countAB / countAll
            self.pAABB = self.pAA + self.pBB

    def calc(self):
        """Runs class functions that calculates polymer statistics for polymer chains and reaction extent and saves the results as instance variables. 
        """
        self.buildchain()
        self.literalchain()
        self.chainlength()
        self.makeseq()
        self.reactions()
        self.countpairs()
    
    def getrecord(self):
        """Returns all calculated instance variables generated by the calc method and returns them as a dictionary.

        Returns
        -------
        dict
            A dictionary containing all calculated polymer statistics generated by the calc method.
        """
        return {
            "pAA": self.pAA,
            "pBB": self.pBB,
            "pAABB": self.pAABB,
            "pAB": self.pAB,
            "PDI": self.pdi,
            "seq": self.seqlist,
            "monoAcount": self.monoAcount,
            "monoBcount": self.monoBcount,
            "extent": self.extent,
            "NumAve": self.numave,
            "MassAve": self.massave,
            "chaindist": self.chaindist,
            "blockAdist": self.blockAdist,
            "blockBdist": self.blockBdist
        }

    def gettraj(self):
        """Returns the a subset of calculated instance variables generated by the calc method function related 
        to the trajectory of the polymerization reaction (i.e. reaction extent, PDI, number average, ...)
 
        Returns
        -------
        dict
            A dictionary containing a subset of calculated instance variables generated by the calc method.
        """
        return {
            "pAA": self.pAA,
            "pBB": self.pBB,
            "pAABB": self.pAABB,
            "pAB": self.pAB,
            "PDI": self.pdi,
            "extent": self.extent,
            "NumAve": self.numave,
            "MassAve": self.massave
        }
