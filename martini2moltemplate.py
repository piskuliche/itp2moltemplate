import numpy as np

class atom:
    def __init__(self,line):
        line = line.strip().split()
        self.name   =   str(line[0])
        self.M      =   float(line[1])
        self.Q      =   float(line[2])
        self.ptype  =   str(line[3])
        self.flj    =   float(line[4])
        self.fqq    =   float(line[5])

class nbparam:
    def __init__(self,line):
        line = line.strip().split()
        self.name = str(line[1])
        self.sig  = float(line[2])
        self.eps  = float(line[3])
        self.pairs = []
        self.funcs = []
    def add_interaction(self, atom1, atom2,func):
        self.pairs.append([atom1,atom2])
        self.funcs.append(func)
    def convert_to_real(self):
        self.sig = self.sig*10.0

class molecule:
    def __init__(self,line):
        line = line.strip().split()
        self.name = str(line[0])
        self.nrexcl = int(line[1])
        self.atoms={"i":[],"type":[],"resnr":[],"residue":[],"atom":[],"cgnr":[],"Q":[]}
    def add_atom(self,line):
        line = line.strip().split()
        c=0
        intkey=["i","resnr","cgnr"]
        fltkey=["Q"]
        strkey=["type","residue","atom"]
        for key in self.atoms:
            if key in intkey: self.atoms[key].append(int(line[c]))
            if key in fltkey: self.atoms[key].append(float(line[c]))
            if key in strkey: self.atoms[key].append(str(line[c]))
            c = c+1
        

def is_not_blank(s):
    return bool(s and not s.isspace())

def read_main_itp(filename):
    with open(filename,'r') as f:
        lines=f.readlines()
        flag_atoms, flag_nb = 0,0
        atoms,nbparams = {},{}
        newtype=0
        # Looks for Definitions
        for line in lines:
            if "#define" in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    nb = nbparam(l)
                    if nb.name in nbparams:
                        exit("%s is doubly defined!" % nb.name)
                    nbparams[nb.name] = nb
        # Loops and Saves Info
        for line in lines:
            # Finds atom type definitions
            if flag_atoms == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    atm = atom(l)
                    atoms[atm.name] = atm
            # Finds pairwise interactions, and if necessary, sorts them
            if flag_nb == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    l = l.strip().split()
                    if l[3] in nbparams:
                        nbparams[l[3]].add_interaction(str(l[0]),str(l[1]),int(l[2]))
                    else:
                        sig, eps = float(l[3]), float(l[4])
                        found=0
                        for param in nbparams:
                            if sig == nbparams[param].sig and eps == nbparams[param].eps:
                                nbparams[param].add_interaction(str(l[0]),str(l[1]),str(l[2]))
                                found = 1
                        if found == 0:
                            nbparams[str(newtype)]=nbparam("#define %s %10.5f %10.5f" % (str(newtype),sig,eps))
                            newtype += 1
            if "[ atomtypes ]" in line:
                flag_atoms  =   1
                flag_nb     =   0
            if "[ nonbond_params ]" in line:
                flag_nb     =   1
                flag_atoms  =   0

    print("Read primary ITP file!")
    print("There are %d types of atoms" % len(atoms))
    print("There are %d types of pairwise interactions" % len(nbparams))
    return

def read_mol_itp(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
        flag_atoms, flag_bonds, flag_angles, flag_mol = 0,0,0,0
        definitions = {}
        molecules = {}
        currentmol=None
        for line in lines:
            # This stores the definitions in a dictionary for later
            if "#define" in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    definitions[str(l.strip().split()[1])] = l
        for line in lines:
            if flag_mol == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    newmol = molecule(l)
                    molecules[newmol.name] = newmol
                    currentmol = newmol.name
            if flag_atoms == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    molecules[currentmol].add_atom(l)
            if "[ moleculetype ]" in line:
                flag_atoms, flag_bonds, flag_angles,flag_mol = 0,0,0,1
            if "[ atoms ]" in line:
                 flag_atoms, flag_bonds, flag_angles,flag_mol = 1,0,0,0
            if "[ bonds ]" in line:
                 flag_atoms, flag_bonds, flag_angles,flag_mol = 0,1,0,0
            if "[ angles ]" in line:
                 flag_atoms, flag_bonds, flag_angles,flag_mol = 0,0,1,0
        print("There are %d molecules" % len(molecules))

                

read_main_itp("Dry-Martini/dry_martini_v2.1.itp")
read_mol_itp("Dry-Martini/dry_martini_v2.1_cholesterol.itp")
