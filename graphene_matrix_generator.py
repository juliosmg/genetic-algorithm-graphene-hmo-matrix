import os
import numpy as np
import pandas as pd


# size: number of holes in the pattern, e.g. 2x2 -> 2, 3x3 -> 3, etc

def TotalNumber(size):
    return 2*size*(size+2)


def NonSequentials(size):
    NStart = size*2
    NStep = size*2 + 2
    return list(range(NStart, TotalNumber(size), NStep))


def ExceptionalBonds(size):
    oof = CornerBonds(size)[1]
    return [oof-1,TotalNumber(size)-oof]


def CornerBonds(size):
    return [0] + [i+1 for i in NonSequentials(size) if i != NonSequentials(size)[-1]]


def NitrogenPositions(size):    
    positions = []
    cb = CornerBonds(size)
    ns = NonSequentials(size)
    total = TotalNumber(size)
    for k in range(len(cb)):
        for i in range(cb[k],ns[k]+1,2):
            positions.append(i)

    for k in range(total - ns[0],total,2):
        positions.append(k)
    
    return positions

def UnitCells(size):
    positions = []
    except_ = ExceptionalBonds(size)
    for k in NitrogenPositions(size):
        if k in except_:
            positions.append((k, except_[1]))
        else:
            positions.append((k, k+1))
    
    return positions



def SecondaryBonds(size):    
    TBonds = []
    cb = CornerBonds(size)
    ns = NonSequentials(size)
    for i in range(len(cb)):
        oof = []
        if i == 0:
            for j in range(cb[i],ns[i]+2,2):
                oof.append(j)
            # print(oof)
        else:
            for j in range(cb[i],ns[i],2):
                oof.append(j)
            # print(oof)
        
        TBonds.append(oof)
    return TBonds

def ConnectedBonds(size):
    TBonds = SecondaryBonds(size)
    Conds =  []
    for j in range(len(TBonds)):
        oof = []
        if j == 0 or j == len(TBonds)-1:
            bind = 2*(size + 1)
            for k in TBonds[j]:
                oof.append(k+bind)
        else:
            bind = 2*(size + 1) + 1
            for k in TBonds[j]:
                oof.append(k+bind)
        Conds.append(oof)
    
    return Conds

def HexagonalLattice(size):
    NumberOfAtoms = TotalNumber(size)
    ListOfNonSequentialBonds = NonSequentials(size)
    
    TriangularMatrix = np.zeros((NumberOfAtoms,NumberOfAtoms))

    for i in range(NumberOfAtoms):
        if i not in ListOfNonSequentialBonds:
            try:
                TriangularMatrix[i,i+1] = 1
            except:
                pass
            
    TBonds = SecondaryBonds(size)
    for j in range(len(TBonds)):
        if j == 0 or j == len(TBonds)-1:
            bind = 2*(size + 1)
            for k in TBonds[j]:
                TriangularMatrix[k,k+bind] = 1
        else:
            bind = 2*(size + 1) + 1
            for k in TBonds[j]:
                TriangularMatrix[k,k+bind] = 1
        
            
    return TriangularMatrix + TriangularMatrix.T


if __name__ == '__main__':
    size = int(input('ZigZag Hexagonal Lattice Dimension, e.g.: 2x2 -> 2, 3x3 -> 3, etc..., must be greater than 1 (one): '))
    assert size > 1
    NumberOfAtoms = TotalNumber(size)
    path = os.getcwd()
    print('matrix saved in',path)
    pd.DataFrame(HexagonalLattice(size)).to_excel(path + f'//graphene {size}x{size}, {NumberOfAtoms} atoms.xlsx')