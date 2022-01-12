import numpy as np
import os
import sys,copy

TROUBLETYPE = -1
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
        self.eps = self.eps/4.184

class bondparam:
    def __init__(self,line):
        line = line.strip().split()
        self.a1 = str(line[0])
        self.a2 = str(line[1])
        self.R  = float(line[2])
        self.K  = float(line[3])
        self.name="%s-%s"%(self.a1,self.a2)
    def convert_to_real(self):
        self.K = self.K/418.4
        self.K = self.K/2 # for lammps it doesn't include 1/2 factor
        self.R = self.R*10

class anglparam:
    def __init__(self,line):
        line=line.strip().split()
        self.a1 = str(line[0])
        self.a2 = str(line[1])
        self.a3 = str(line[2])
        self.name = "%s-%s-%s"%(self.a1,self.a2,self.a3)
        self.K  = float(line[4])
        self.THETA = float(line[3])
    def convert_to_real(self):
        self.K = self.K/4.184
        self.K = self.K/2

class molecule:
    def __init__(self,line):
        line = line.strip().split()
        self.name = str(line[0])
        self.nrexcl = int(line[1])
        self.atoms={"i":[],"type":[],"resnr":[],"residue":[],"atom":[],"cgnr":[],"Q":[],"x":[],"y":[],"z":[]}
        self.bonds=[]
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
            if key in ["x","y","z"]: self.atoms[key].append(0)
            c = c+1
    def add_bond(self,line):
        line = line.strip().split()
        a1=int(line[0])-1
        a2=int(line[1])-1
        self.bonds.append([self.atoms["atom"][a1],self.atoms["atom"][a2]])
        return
    def write(self,ffname):
        with open("ltfiles/"+self.name+".lt",'w') as g:
            g.write('import "%s.lt"\n\n'%ffname)
            g.write('%s inherits %s {\n' % (self.name,ffname))
            g.write('write("Data Atoms") {\n')
            for i in range(len(self.atoms["i"])):
                g.write("$atom:%s $mol:. @atom:%s % 10.8f % 10.8f % 10.8f % 10.8f\n" % (self.atoms["atom"][i],self.atoms["type"][i],self.atoms["Q"][i],self.atoms["x"][i],self.atoms["y"][i],self.atoms["z"][i]))
            g.write("}\n\n")
            g.write('write("Data Bond List") {\n')
            c=0
            for bond in self.bonds:
                g.write("$bond:b%d $atom:%s $atom:%s\n" % (c,bond[0],bond[1]))
                c=c+1
            g.write("}\n\n")
            #g.write("%s.scale(10)\n}\n"%self.name)
            g.write("\n}\n"%self.name)
        return
    def add_coords(self,grofile):
        with open(grofile,'r') as g:
            g.readline()
            g.readline()
            for atom in range(len(self.atoms["i"])):
                line = g.readline().strip().split()
                self.atoms["x"][atom] = float(line[3])*10
                self.atoms["y"][atom] = float(line[4])*10
                self.atoms["z"][atom] = float(line[5])*10
            print("Wrote coords, assumed nanometer and converted to A")




def is_not_blank(s):
    return bool(s and not s.isspace())

def read_main_itp(filename):
    with open(filename,'r') as f:
        lines=f.readlines()
        flag_atoms, flag_nb = 0,0
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
                        if str(l[0]) == "Q0" and str(l[1]) == "Q0": print("test",l[3])
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
            if "[atomtypes]" in line.replace(" ", ""):
                flag_atoms  =   1
                flag_nb     =   0
            if "[nonbond_params]" in line.replace(" ", ""):
                flag_nb     =   1
                flag_atoms  =   0

    print("Read primary ITP file!")
    print("There are %d types of atoms" % len(atoms))
    print("There are %d types of pairwise interactions" % len(nbparams))
    return atoms, nbparams



def read_mol_itp(filename,TROUBLETYPE):
    print("Beginning read of %s" % filename)
    with open(filename,'r') as f:
        lines = f.readlines()
        flag_atoms, flag_bonds, flag_angles, flag_mol = 0,0,0,0
        definitions = {}
        currentmol=None
        newmol=None
        for line in lines:
            # This stores the definitions in a dictionary for later
            if "#define" in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    definitions[str(l.strip().split()[1])] = l
        for line in lines:
            # Starts reading molecule
            if flag_mol == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    newmol = molecule(l)
                    #molecules[newmol.name] = newmol
                    currentmol = newmol.name
            # Starts reading atom section
            if flag_atoms == 1 and "[" not in line:
                l = line.split(";")[0]
                if is_not_blank(l):
                    newmol.add_atom(l)
            # Starts reading bonds section            
            if flag_bonds == 1 and "[" not in line:
                l = line.split(";")[0]
                l = l.split("#")[0]
                if is_not_blank(l):
                    newmol.add_bond(l)
                    l=line.strip().split()
                    r,k=0,0
                    atom1, atom2 = newmol.atoms["type"][int(l[0])-1],newmol.atoms["type"][int(l[1])-1]
                    print(l)
                    if l[3] in definitions:
                        r,k=float(definitions[l[3]].strip().split()[2]),float(definitions[l[3]].strip().split()[3])
                    else:
                        r,k=float(l[3]),float(l[4])
                    name = "%s-%s"%(atom1,atom2)
                    if name in bonds:
                        if bonds[name].R != r:
                            print("Trouble with %s redefining bond length of %s from %10.8f to %10.8f"%(currentmol,name,bonds[name].R,r))
                            if currentmol not in troublemolecs:
                                TROUBLETYPE += 1
                                troublemolecs.append(currentmol)
                        if bonds[name].K != k:
                            print("Trouble with %s redefining bond const  of %s from %10.8f to %10.8f"%(currentmol,name,bonds[name].K,k))
                            if currentmol not in troublemolecs:
                                TROUBLETYPE +=1
                                troublemolecs.append(currentmol)
                        if currentmol in troublemolecs:
                            name = "%s-%s"%(atom1+"f"+str(TROUBLETYPE),atom2+"f"+str(TROUBLETYPE))
                            bonds[name]=bondparam("%s %s %10.8f %10.8f"%(atom1+"f"+str(TROUBLETYPE),atom2+"f"+str(TROUBLETYPE),r,k))
                    else:
                        bonds[name]=bondparam("%s %s %10.8f %10.8f"%(atom1,atom2,r,k))
            # Starts reading angles section
            if flag_angles == 1 and "[" not in line:
                l = line.split(";")[0]
                l = l.split("#")[0]
                if is_not_blank(l):
                    l=line.strip().split()
                    theta,k = 0,0
                    atom1 = newmol.atoms["type"][int(l[0])-1]
                    atom2 = newmol.atoms["type"][int(l[1])-1]
                    atom3 = newmol.atoms["type"][int(l[2])-1]
                    if l[4] in definitions:
                        theta, k = float(definitions[l[4]].strip().split()[2]),float(definitions[l[4]].strip().split()[3])
                    else:
                        theta, k = float(l[4]),float(l[5])
                    name = "%s-%s-%s"%(atom1,atom2,atom3)
                    if name in angles:
                        if angles[name].THETA != theta:
                            print("Trouble with %s redefining angle of %s from %10.8f to %10.8f"%(currentmol,name,angles[name].THETA,theta))
                            if currentmol not in troublemolecs:
                                TROUBLETYPE += 1
                                troublemolecs.append(currentmol)
                        if angles[name].K != k:
                            print("Trouble with %s redefining angle const  of %s from %10.8f to %10.8f"%(currentmol,name,angles[name].K,k))
                            if currentmol not in troublemolecs:
                                TROUBLETYPE += 1
                                troublemolecs.append(currentmol)
                        if currentmol in troublemolecs:
                            name = "%s-%s-%s"%(atom1+"f"+str(TROUBLETYPE),atom2+"f"+str(TROUBLETYPE),atom3+"f"+str(TROUBLETYPE))
                            angles[name]=anglparam("%s %s %s %10.8f %10.8f"%(atom1+"f"+str(TROUBLETYPE),atom2+"f"+str(TROUBLETYPE),atom3+"f"+str(TROUBLETYPE),theta,k))
                    else:
                        angles[name]=anglparam("%s %s %s %10.8f %10.8f"%(atom1,atom2,atom3,theta,k))

            if "[moleculetype]" in line.replace(" ", ""):
                if newmol is not None: molecules[newmol.name]=newmol
                newmol = None
                flag_atoms, flag_bonds, flag_angles,flag_mol = 0,0,0,1
            if "[atoms]" in line.replace(" ", ""):
                flag_atoms, flag_bonds, flag_angles,flag_mol = 1,0,0,0
            if "[bonds]" in line.replace(" ", ""):
                flag_atoms, flag_bonds, flag_angles,flag_mol = 0,1,0,0
            if "[angles]" in line.replace(" ", ""):
                flag_atoms, flag_bonds, flag_angles,flag_mol = 0,0,1,0
            if "[constraints]" in line.replace(" ","") or "[dihedrals]" in line.replace(" ",""):
                flag_atoms, flag_bonds, flag_angles, flag_mol = 0,0,0,0
        if newmol is not None: molecules[newmol.name]=newmol
        print("There are %d molecules" % len(molecules))
        print("There are %d bond types" % len(bonds))
        print("There are %d angle types" % len(angles))
    return TROUBLETYPE

def write_mass(f,atoms):
    f.write('write_once("Data Masses") {\n')
    for atm in atoms:
        f.write("@atom:%s %10.6f\n" % (atoms[atm].name, atoms[atm].M))
    f.write("}\n\n")
    return

def write_pair(f, paircoeffs):
    f.write("\n")
    f.write('write_once("In Settings") {\n')
    for p in paircoeffs:
        pairtype = paircoeffs[p]
        pairtype.convert_to_real()
        for pair in pairtype.pairs:
            f.write("pair_coeff @atom:%s @atom:%s lj/gromacs/coul/gromacs %10.8f %10.8f\n" % (pair[0],pair[1],pairtype.eps,pairtype.sig))
            if pair[0] in convertname:
                for elem in convertname[pair[0]]:
                    f.write("pair_coeff @atom:%s @atom:%s lj/gromacs/coul/gromacs %10.8f %10.8f\n" % (elem,pair[1],pairtype.eps,pairtype.sig))
            if pair[1] in convertname:
                for elem in convertname[pair[1]]:
                    f.write("pair_coeff @atom:%s @atom:%s lj/gromacs/coul/gromacs %10.8f %10.8f\n" % (pair[0],elem,pairtype.eps,pairtype.sig))
            if pair[0] in convertname and pair[1] in convertname:
                for elem1 in convertname[pair[0]]:
                    for elem2 in convertname[pair[1]]:
                        f.write("pair_coeff @atom:%s @atom:%s lj/gromacs/coul/gromacs %10.8f %10.8f\n" % (elem1,elem2,pairtype.eps,pairtype.sig))

    f.write("}\n\n")
    return


def write_bond(f,bonds):
    f.write("\n")
    f.write('write_once("In Settings") {\n')
    for b in bonds:
        bondtype = bonds[b]
        bondtype.convert_to_real()
        f.write("bond_coeff @bond:%s harmonic %10.8f %10.8f\n" % (bondtype.name, bondtype.K, bondtype.R))
    f.write("}\n\n")
    
    f.write('write_once("Data Bonds By Type") {\n')
    for b in bonds:
        bondtype = bonds[b]
        f.write("@bond:%s @atom:%s @atom:%s\n"%(bondtype.name,bondtype.a1,bondtype.a2))
    f.write("}\n\n")
    return

def write_angl(f,angles):
    f.write("\n")
    f.write('write_once("In Settings") {\n')
    for a in angles:
        angtype = angles[a]
        angtype.convert_to_real()
        f.write("angle_coeff @angle:%s cosine/squared %10.8f %10.8f\n" % (angtype.name, angtype.K, angtype.THETA))
    f.write("}\n\n")

    f.write('write_once("Data Angles By Type") {\n')
    for a in angles:
        angtype = angles[a]
        f.write("@angle:%s @atom:%s @atom:%s @atom:%s\n"%(angtype.name,angtype.a1,angtype.a2,angtype.a3))
    f.write("}\n\n")
    return

def write_sett(f):
    f.write("\n")
    f.write('write_once("In Init") {\n')
    f.write("units real\natom_style full\nbond_style hybrid harmonic\n")
    f.write("angle_style hybrid cosine/squared\n")
    f.write("pair_style hybrid lj/gromacs/coul/gromacs 9 12 0.000001 12\n")
    f.write("special_bonds lj/coul 0.0 1.0 1.0\n")
    f.write("dielectric 15.0\n")
    f.write("}\n\n")



def write_ff(fname,atoms, paircoeffs,bonds,angles):
    with open (fname+".lt",'w') as f:
        f.write("%s {\n" % fname)
        write_mass(f,atoms)
        write_pair(f, paircoeffs)
        write_bond(f, bonds)
        write_angl(f, angles)
        write_sett(f)
        f.write("\n}\n")
    return

def fix_trouble(troublemolecs):
    # Fixes molecule class details
    # Adds new atom types
    if len(troublemolecs) == 0:
        return
    c = 0
    for molec in troublemolecs:
        print("Fixing molecule %s" % molec)
        mol = molecules[molec]
        newmol = copy.deepcopy(mol)
        # Fix molecule itself
        for i in range(len(mol.atoms["type"])):
            newmol.atoms["type"][i] = mol.atoms["type"][i]+"f"+str(c)
            atoms[mol.atoms["type"][i]+"f"+str(c)]=copy.deepcopy(atoms[mol.atoms["type"][i]]) # copies to new entry
            atoms[mol.atoms["type"][i]+"f"+str(c)].name=mol.atoms["type"][i]+"f"+str(c)
            if mol.atoms["type"][i] not in convertname:
                convertname[mol.atoms["type"][i]]=[]
                convertname[mol.atoms["type"][i]].append(mol.atoms["type"][i]+"f"+str(c))
            else:
                if mol.atoms["type"][i]+"f"+str(c) not in  convertname[mol.atoms["type"][i]]:
                    convertname[mol.atoms["type"][i]].append(mol.atoms["type"][i]+"f"+str(c))
        for b in range(len(mol.bonds)):
            newmol.bonds[b][0] = mol.bonds[b][0]+"f"+str(c)
            newmol.bonds[b][1] = mol.bonds[b][1]+"f"+str(c)
        c=c+1
        molecules[molec] = newmol
        
if len(sys.argv) < 2:
    print("Usage martini2moltemplate.py ffname [molecule to write] [molecule gro file]")
    exit()

fname = str(sys.argv[1])
mol2write, grofile = None,None
if len(sys.argv) == 4:
    mol2write=str(sys.argv[2])
    grofile = str(sys.argv[3])
if not os.path.exists("ltfiles"):
    os.makedirs("ltfiles")
atoms,nbparams,molecules, bonds, angles = {},{},{},{},{}
troublemolecs,convertname=[],{}
#read_main_itp("Dry-Martini/dry_martini_v2.1.itp")
read_main_itp("Dry-Martini/dry_martini_v2.1.itp")
TROUBLETYPE=read_mol_itp("Dry-Martini/dry_martini_v2.1_lipids.itp",TROUBLETYPE)
TROUBLETYPE=read_mol_itp("Dry-Martini/dry_martini_v2.1_solvents.itp",TROUBLETYPE)
TROUBLETYPE=read_mol_itp("Dry-Martini/dry_martini_v2.1_cholesterol.itp",TROUBLETYPE)
TROUBLETYPE=read_mol_itp("Dry-Martini/dry_martini_v2.1_ions.itp",TROUBLETYPE)
TROUBLETYPE=read_mol_itp("Dry-Martini/addmol.itp",TROUBLETYPE)
fix_trouble(troublemolecs)
write_ff(fname,atoms,nbparams,bonds,angles)
if len(sys.argv) == 4:
    molecules[mol2write].add_coords(grofile)
    molecules[mol2write].write(fname)

