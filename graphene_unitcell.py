
import os
from collections import Counter
from scipy.linalg import eigh, norm,block_diag
from numpy.random import choice, randint
import numpy as np
import pandas as pd

# HMO PARAMETERS
hb = -1.0
hn = 1.5
hc = 0.0
kcc = 1.0
kbn = 0.9
kbc = 0.7
kcn = 0.8

import graphene_matrix_generator as gmg
from random import choices
from sympy.solvers import solve
from sympy import Symbol
import sys


def rand_int(size):
    # Determine the number of bytes needed
    num_bytes = size.bit_length() // 8 + 1
    
    rand = int.from_bytes(os.urandom(num_bytes), sys.byteorder)
    while rand >= size:
        rand = int.from_bytes(os.urandom(num_bytes), sys.byteorder)
    
    return rand



# Calculate HMO properties of the molecule
def hmo(matrix):
    # matrix = self.matrix
    dict = {}
    p = 5
    beta = 2.5 # eV
    # energies and coefficients calculations
    eigen, vectors = eigh(-matrix)
    vectors = np.around(vectors,p).T
    eigen = np.around(eigen,p)
    # energies = eigen*beta
    
    # gap = -(HOMO - LUMO)
    pielectrons = (len(matrix) + sum(np.diagonal(matrix) == hn) - sum(np.diagonal(matrix) == hb))//2
    gap = -(eigen[pielectrons-1] - eigen[pielectrons])*beta
    gap = np.around(gap,p)

    # IPNs calculations
    ipn = [sum(vectors[i]**4) * sum(vectors[i]**2)**-2 for i in range(len(vectors))]
    ipn = np.around(ipn[pielectrons-1],p)


    total = sum(abs(eigen))
    dict['energies'] = eigen
    # dict['energies'] = energies
    dict['eigenvectors'] = vectors
    dict['\u03C0-electrons, Total \u03C0-Energy'] = (pielectrons,np.round(total,p))
    dict['gap'] = gap
    dict['homo'] = pielectrons-1
    dict['ipn'] = ipn
    # dict['fit'] = np.around(.7*gap + 0.3*ip,p)
    dict['counter'] = Counter(np.diagonal(matrix))
    

    return dict

# Generate a tight-binding matrix from a unitcell string
def tb_matrix_unit(genome):
        
    x = Symbol('x')
    size = [int(j) for j in solve(gmg.TotalNumber(x) - len(genome)*2, x) if j > 0][0]   
    
    skt = gmg.HexagonalLattice(size)
    
    unitcells = gmg.UnitCells(size)    
    
    # Coulomb parameters correction
    diag = [0]*len(genome)*2
    for i in range(len(genome)):
        if genome[i] != 0:
            diag[unitcells[i][0]] = hn
            diag[unitcells[i][1]] = hb
        else:
            diag[unitcells[i][0]] = hc
            diag[unitcells[i][1]] = hc

        
    # Resonance parameters correction    
    for i in gmg.NitrogenPositions(size):
        hood = np.where(skt[i] != 0)[0] # primeiros vizinhos de i : posições
        for j in hood:
            if diag[i] == hn and diag[j] == hc:
                skt[i,j], skt[j,i] = kcn, kcn
            elif diag[i] == hn and diag[j] == hb:
                skt[i,j], skt[j,i] = kbn, kbn
            elif diag[i] == hc and diag[j] == hc:
                skt[i,j], skt[j,i] = kcc, kcc
            elif diag[i] == hc and diag[j] == hb:
                skt[i,j], skt[j,i] = kbc, kbc
                                    
    matrix = np.zeros((len(diag),len(diag)))
    np.fill_diagonal(matrix,diag)
    
    return np.matrix(skt + matrix)



def fitness(bandgaps):
    weights = np.zeros(len(bandgaps))
    for k in range(len(bandgaps)):
        if bandgaps[k] >= 0.5 and bandgaps[k] <= 1.0:
            weights[k] = 1e-16
        elif bandgaps[k] > 1.0:
            weights[k] = np.abs(bandgaps[k] - 1.0)
        else:
            weights[k] = np.abs(bandgaps[k] - 0.5)
    weights = np.power(weights,-1)
    # weights = weights/sum(weights)
    
    return np.round(weights,1)
    

# Proportional selection function
def selection(bandgaps):
    # chooses positions pondered by the inverse of the given weight 'w'
    w = fitness(bandgaps)
    w = w/sum(w)
    options = choice(a = list(range(len(bandgaps))),size =2,replace= False,p=w).tolist()
    return options[0],options[1]


# Reproduction function
def crossover(population,bandgaps):
    # single-point crossover 
    # input: population and band gaps
    # output: 2 new chromossomes
    x, y = selection(bandgaps)
    gene1 = population[x]
    gene2 = population[y]
    k = rand_int(len(gene1))
    while k == 0:
        k = rand_int(len(gene1))
        
    return [*gene1[:k], *gene2[k:]],[*gene2[:k], *gene1[k:]]


# Evolution of population
def evolution(population):    
    dict = {}
    neon = {} 
    generation = 0

    # first evaluation
    matrices = [tb_matrix_unit(i) for i in population]
    diags = [np.diagonal(i).tolist() for i in matrices]
    hmodict = [hmo(i) for i in matrices]
    gaps = [i['gap'] for i in hmodict]
    fit = fitness(gaps).tolist()
    
    # avaliation & append results
    max_ = np.max(fit)
    if max_ == 1e16:
        mx = np.where(fit == max_)[0]
        solutions = []
        pops = []
        pdiags = []
        pfit = []
        pgaps = []
        pipn = []
        phomo = []
        
        for k in mx:
            r = hmodict[k]['counter']
            carbon = r[0.0]
            boron = r[-1.0]
            nitrogen = r[1.5]
            solution = '$C_{'+str(carbon)+'} N_{'+str(nitrogen)+'} B_{'+str(boron)+'}$'
            solutions.append(solution)
            pops.append(population[k])
            pdiags.append(diags[k])
            pfit.append(fit[k])
            pgaps.append(gaps[k])
            pipn.append(hmodict[k]['ipn'])
            phomo.append(hmodict[k]['homo'])
                
        dict['best_solution'] = solutions
        dict['genome'] = pops
        dict['diag'] = pdiags
        dict['fit'] = pfit
        dict['gap'] = pgaps
        dict['ipn'] = pipn
        dict['homo'] = phomo
        dict['population'] = population
        neon[generation] = dict
    else:
        mx = fit.index(max_)
        r = hmodict[mx]['counter']
        carbon = r[0.0]
        boron = r[-1.0]
        nitrogen = r[1.5]
        solution = '$C_{'+str(carbon)+'} N_{'+str(nitrogen)+'} B_{'+str(boron)+'}$'
        solutions = solution
        dict['best_solution'] = solutions
        dict['genome'] = population[mx]
        dict['diag'] = diags[mx]
        dict['fit'] = fit[mx]
        dict['gap'] = gaps[mx]
        dict['ipn'] = hmodict[mx]['ipn']
        dict['homo'] = hmodict[mx]['homo']
        dict['population'] = population
        neon[generation] = dict


    # covergence?
    convergence = all([0.5 <= x <= 1.0 for x in gaps])

    while not convergence:
        dict = {}
        generation += 1
        mx = fit.index(max_)
        new = [population[mx]]

        
        # reproduction
        print('reproducing...')
        for k in range(len(population)//2):
            f,g = crossover(population, gaps)
            new.append(f),new.append(g)
        
        
        
        # Homogeneity
        h = [new[0] == k for k in new]
        # Elitism and mutation
        if all(h):
            print('Applying elitism and mutation')
            
            ### Mutation ###
            # Fast convergence
            new = [randint(2,size=len(population[0])).tolist() for i in range(len(population))]
            
            # ## Slow convergence
            # for i in new:
            #     i[rand_int(len(i))] = 1
                
            
            # ## Random convergence
            # for i in range(len(new)):
            #     k = rand_int(len(new[0])//10)
            #     try:
            #         for j in range(k):
            #             new[i][rand_int(len(new[0]))] = 1
            #     except:
            #         pass
                

            ### Elitism ###
            try:
                for k in range(len(mx)):
                    new[k] = population[mx[k]]
            except:
                new[0] = population[mx]
                
        else:
            pass

        population = new.copy()
                
        # children evaluation
        matrices = [tb_matrix_unit(i) for i in population]
        diags = [np.diagonal(i).tolist() for i in matrices]
        hmodict = [hmo(i) for i in matrices]
        gaps = [i['gap'] for i in hmodict]
        fit = fitness(gaps).tolist()
        
        # avaliation & append results
        max_ = np.max(fit)
        if max_ == 1e16:
            mx = np.where(fit == max_)[0]
            solutions = []
            pops = []
            pdiags = []
            pfit = []
            pgaps = []
            pipn = []
            phomo = []
            
            for k in mx:
                r = hmodict[k]['counter']
                carbon = r[0.0]
                boron = r[-1.0]
                nitrogen = r[1.5]
                solution = '$C_{'+str(carbon)+'} N_{'+str(nitrogen)+'} B_{'+str(boron)+'}$'
                solutions.append(solution)
                pops.append(population[k])
                pdiags.append(diags[k])
                pfit.append(fit[k])
                pgaps.append(gaps[k])
                pipn.append(hmodict[k]['ipn'])
                phomo.append(hmodict[k]['homo'])
                    
            dict['best_solution'] = solutions
            dict['genome'] = pops
            dict['diag'] = pdiags
            dict['fit'] = pfit
            dict['gap'] = pgaps
            dict['ipn'] = pipn
            dict['homo'] = phomo
            dict['population'] = population
            neon[generation] = dict
        else:
            mx = fit.index(max_)
            r = hmodict[mx]['counter']
            carbon = r[0.0]
            boron = r[-1.0]
            nitrogen = r[1.5]
            solution = '$C_{'+str(carbon)+'} N_{'+str(nitrogen)+'} B_{'+str(boron)+'}$'
            solutions = solution
            dict['best_solution'] = solutions
            dict['genome'] = population[mx]
            dict['diag'] = diags[mx]
            dict['fit'] = fit[mx]
            dict['gap'] = gaps[mx]
            dict['ipn'] = hmodict[mx]['ipn']
            dict['homo'] = hmodict[mx]['homo']
            dict['population'] = population
            neon[generation] = dict

        # covergence?
        convergence = all([0.5 <= x <= 1.0 for x in gaps])
        if convergence:
            print('Convergence!') 

        
        if generation == 100:
            print('Max generation reached!')
            break
      
        
    return neon


path = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    size = int(input('ZigZag Hexagonal Lattice Dimension, e.g.: 2x2 -> 2, 3x3 -> 3, etc..., must be greater than 1 (one): \'int\'= '))
    popsize = 5
    tn = gmg.TotalNumber(size)//2
    population = [[0]*tn for i in range(popsize)]
    # max = int(input('limit of generations per trial: \'int\'= '))
    t = int(input('number of trials \'int\'= '))

    results = {}
    for i in range(t):
        print('\n'+10*'##')
        ga = evolution(population)
        results[i] = ga
        print('\n'+10*'##')


    data = pd.DataFrame.from_dict({(i,j): results[i][j]
                                for i in results.keys()
                                for j in results[i].keys()},
                                orient='index')

    data.to_excel(path + '\\results two-dimensional '+str(size)+'x'+str(size)+' system, unit evolution.xlsx')
    print('saved in ', path)






