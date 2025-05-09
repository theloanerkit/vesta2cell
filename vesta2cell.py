import argparse

# set up command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s","--seed",help="seedname for input .vesta file",required=True)
parser.add_argument("-sp","--spin",help="treatment of spin",default="nospin",choices=["nospin","collinear","noncollinear"])
parser.add_argument("-o","--out",help="seedname for output .cell file (optional)")

class symop:
    """holds a single symmetry operation from a .vesta SYMOP block
    Args:
        vector (list[float]) : vector symmetry operation
        matrix (list[list[float]]) : matrix symmetry operation
    """
    def __init__(self,vector:list[float],matrix:list[list[float]]):
        self.vector = vector
        self.matrix = matrix

    def apply(self,vector:list[float])->list[float]:
        """applys fully symmetry operation (matrix and vector) to the input vector
        Args:
            vector (list[float]): vector to be transformed
        Returns:
            list[float]: transformed vector
        """
        temp = [0,0,0]
        for i in range(3):
            for j in range(3):
                temp[i] += self.matrix[i][j]*vector[j]
            temp[i] += self.vector[i]
        return temp

    def apply_mat(self,vector:list[float])->list[float]:
        """applys the matrix symmetry operation to the input vector

        Args:
            vector (list[float]): vector to be transformed

        Returns:
            list[float]: transformed vector
        """
        temp = [0,0,0]
        for i in range(3):
            for j in range(3):
                temp[i] += self.matrix[i][j]*vector[j]
        return temp

class vesta_symop:
    """holds all symmetry operations from a .vesta SYMOP block
    Args:
        data (list[str]) : SYMOP block from .vesta file
    Attributes:
        symops (list[symop]) : 
    """
    def __init__(self,data:list[str]):
        self.data = data
        self.symops = []
        self.parse()

    def parse(self):
        """parses a .vesta SYMOP block
        """
        for line in self.data:
            sep = line.split()
            vec = list(map(float,sep[0:3]))
            mat = [list(map(float,sep[3:6])),
                   list(map(float,sep[6:9])),
                   list(map(float,sep[9:12]))]
            s = symop(vec,mat)
            self.symops.append(s)

class vesta_cellp:
    """holds all cell parameters from a .vesta CELLP block
    Args:
        data (list[str]) : CELLP block from .vesta file
    Attributes:
        params (list[float]) : magnitude of the lattice vectors a,b,c
        angles (list[float]) : angles between the lattice vectors α,β,γ
    """
    def __init__(self,data:list[str]):
        self.data = data
        self.params = []
        self.angles = []
        self.parse()

    def parse(self):
        """parse a .vesta CELLP block
        """
        self.params = list(map(float,self.data[0].split()[0:3]))
        self.angles = list(map(float,self.data[0].split()[3:]))

class vesta_struc:
    """holds all atomic positions from a .vesta STRUC block
    Args:
        data (list[str]) : STRUC block from a .vesta file
    Attributes:
        atoms (list[str]) : element symbols
        positions (list[list[float]]) : atomic positions
        symops (list[int]) : how many symops apply to each atomic position
    """
    def __init__(self,data):
        self.data = data
        self.atoms = []
        self.positions = []
        self.symops = []
        self.parse()

    def parse(self):
        """parses a .vesta STRUC block
        """
        for line in self.data:
            if "." not in line.split()[0] and line[0] != "0":
                self.atoms.append(line.split()[1])
                self.positions.append(list(map(float,line.split()[4:7])))
                self.symops.append(int(line.split()[7]))

class vesta_vectr:
    """holds all spin vectors from a .vesta VECTR block
    Args:
        data (list[str]) : VECTR block from a .vesta file
    Attributes:
        vectors (list[list[float]]) : spin vectors
    """
    def __init__(self,data):
        self.data = data
        self.vectors = []
        self.parse()

    def parse(self):
        """parses a .vesta VECTR block
        """
        for i in range(0,len(self.data)-2,3):
            self.vectors.append(list(map(float,self.data[i].split()[1:4])))

# map vesta keywords to class
vesta_block = {"SYMOP": vesta_symop,
               "CELLP": vesta_cellp,
               "STRUC": vesta_struc,
               "VECTR": vesta_vectr}

class castep_lattice_abc:
    def __init__(self,params,angles):
        self.params = params
        self.angles = angles

    def __repr__(self):
        l = get_max_length(self.params,self.angles)
        string = ""
        string += "%block lattice_abc\n"
        for m in self.params:
            string += f"{extend_str(str(m),l)} "
        string += "\n"
        for a in self.angles:
            string += f"{extend_str(str(a),l)} "
        string += "\n"
        string += "%endblock lattice_abc\n"
        return string

class castep_positions_frac:
    def __init__(self,atoms,positions):
        self.atoms = atoms
        self.positions = positions
        self.vectors = None

    def __repr__(self):
        l = get_max_length(*self.positions)
        if self.vectors:
            lv = get_max_length(*self.vectors)
        string = ""
        string += "%block positions_frac\n"
        for i in range(len(self.atoms)):
            string += f"{extend_str(self.atoms[i],3)}"
            for j in range(len(self.positions[i])):
                string += f"{extend_str(str(self.positions[i][j]),l)} "
            if self.vectors:
                string += f"spin= "
                for j in range(len(self.vectors[i])):
                    string += f"{extend_str(str(self.vectors[i][j]),lv)} "
            string += "\n"
        string += "%endblock positions_frac"
        return string

def get_max_length(*arrays):
    arr = []
    l = 0
    for a in arrays:
        arr.extend(a)
    for elem in arr:
        temp = str(elem)
        if len(temp) > l:
            l = len(temp)
    return l

def extend_str(string,length):
    l = len(string)
    return string + " "*(length-l)

def apply_symops_pos(struct,symop):
    atoms = []
    positions = []
    for i in range(len(struct.atoms)):
        for j in range(struct.symops[i]):
            pos = symop.symops[j].apply(struct.positions[i])
            for k in range(len(pos)):
                if pos[k] < 0:
                    pos[k] = 1 + pos[k]
                elif pos[k] > 1:
                    pos[k] -= 1
            atoms.append(struct.atoms[i])
            positions.append(pos)
    return atoms,positions

def apply_symops_vec(vectors,struct,symop,atom_pos):
    atom_pos.vectors = []
    for i in range(len(vectors.vectors)):
        for j in range(struct.symops[i]):
            pos = symop.symops[j].apply_mat(vectors.vectors[i])
            atom_pos.vectors.append(pos)

def read_vesta(fname):
    with open(f"{fname}.vesta","r") as file:
        data = [line.strip() for line in file.readlines()]
    keyword = ""
    data_block = []
    vesta_objs = {}
    for line in data:
        if line.isupper() and not line[0].isnumeric():
            if keyword in vesta_block.keys():
                obj = vesta_block[keyword](data_block)
                vesta_objs[keyword] = obj
            keyword = line
            data_block = []
        else:
            data_block.append(line)
    return vesta_objs

def vesta_to_castep(vesta_objs,spin):
    lattice = castep_lattice_abc(vesta_objs["CELLP"].params,vesta_objs["CELLP"].angles)
    atoms,positions = apply_symops_pos(vesta_objs["STRUC"],vesta_objs["SYMOP"])
    atom_pos = castep_positions_frac(atoms,positions)
    if spin=="noncollinear":
        apply_symops_vec(vesta_objs["VECTR"],vesta_objs["STRUC"],vesta_objs["SYMOP"],atom_pos)
    castep_objs = [lattice,atom_pos]
    return castep_objs

def write_castep(fname,castep_objs):
    with open(f"{fname}.cell","w") as file:
        for obj in castep_objs:
            file.write(repr(obj))
            file.write("\n")

def convert(f_in,f_out,spin="nospin"):
    vesta_objs = read_vesta(f_in)
    castep_objs = vesta_to_castep(vesta_objs,spin)
    write_castep(f_out,castep_objs)

if __name__=="__main__":
    args = parser.parse_args()
    f_in = args.seed
    spin = args.spin
    if args.out:
        f_out = args.out
    else:
        f_out = f_in
    convert(f_in,f_out,spin)