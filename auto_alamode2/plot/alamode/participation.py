import numpy as np

class Participation():
    """
    Variables
    -----------
    .nk : integer
        # of k-points
    .nbands : integer
        # of bands
    kpoints : ndarray, float, shape=(nk,3)
        k-points
    ratio : ndarray, float, shape=(nk,nband)
        Participation ratio of each modes
    """
    def __init__(self, file=None):
        self.nk = None
        self.nbands = None
        self.kpoints = None
        self.ratio = None
        if file is not None:
            self.set_pr_ratio(file)
    
    def set_numbers(self, FILE):
        self.nk, self.nbands = get_numbers(FILE)
        

    def set_pr_ratio(self, FILE):
        if self.nk is None or self.nbands is None:
            self.set_numbers(FILE)
          
        ifs = open(FILE, "r")
        nline = sum(1 for line in open(FILE))
        self.kpoints = np.zeros((self.nk, 3))
        self.ratio = np.zeros((self.nk, self.nbands))
        for il in range(nline):
            line = ifs.readline(); data = line.split()
            if len(data) == 0:
                continue
            if line[0] == "#":
                if data[2] == "xk":
                    ik = int(data[1]) - 1
                    for j in range(3):
                        self.kpoints[ik,j] = data[4+j]
                continue
            ik = int(data[0]) - 1
            ib = int(data[1]) - 1
            self.ratio[ik,ib] = float(data[2])
        ifs.close()

def get_numbers(FILE):
    ifs = open(FILE, "r")
    lines = ifs.readlines()
    nline = len(lines)
    for il in range(nline-1,0,-1):
        data = lines[il].split()
        if len(data) == 0:
            continue
        nk = int(data[0])
        nbands = int(data[1])
        break 
    return nk, nbands

