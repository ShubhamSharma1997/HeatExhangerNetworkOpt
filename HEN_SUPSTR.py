#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 20:31:54 2020

@author: shubham
"""


import numpy as np
import pandas as pd
import random
from deap import base
from deap import creator
from deap import tools



filepath = "/home/shubham/Desktop/Intern/Stream_data.xlsx"
hsdf = pd.read_excel(filepath, sheet_name = "Hot stream")
csdf = pd.read_excel(filepath, sheet_name = "Cold stream")
dTmin = 10
inf = 10000000000
eps = 0.00001


CPh = np.ravel(hsdf['CP'])
CPc = np.ravel(csdf['CP'])

Tih = np.ravel(hsdf['Inlet Temperature'])
Tic = np.ravel(csdf['Inlet Temperature'])

Toh = np.ravel(hsdf['Outlet Temperature'])
Toc = np.ravel(csdf['Outlet Temperature'])

Hh = (Tih-Toh)*CPh
Hc = (Toc-Tic)*CPc

Nh = len(CPh)
Nc = len(CPc)
Ns = max(Nh, Nc)

num_var = Ns*Nc*Nh*(10+Nc+Nh) + (Ns+1)*(Nc+Nh)


def convertToMatrices(individual):
    
    #Extracting alpha variables
    start = 0
    num = Nh*Nc*Ns
    alphah = np.array(individual[start:start+num])
    alphah = alphah.reshape(Nh, Ns, Nc)
    for i in range(Nh):
        alphah[i] = CPh[i]*alphah[i]
    alphah = alphah.transpose(1, 0, 2)
    start = start+num
    alphac = np.array(individual[start:start+num])
    alphac = alphac.reshape(Nc, Ns, Nh)
    for i in range(Nc):
        alphac[i] = CPc[i]*alphac[i]
    alphac = alphac.transpose(1, 0, 2)
    
    #Extracting gamma variables
    start = start+num
    num = Nh*Ns*Nc*Nc
    gammah = np.array(individual[start:start+num])
    gammah = gammah.reshape(Nh, Ns, Nc, Nc)
    for i in range(Nh):
        gammah[i] = gammah[i]*CPh[i]
    gammah = gammah.transpose(1, 0, 2, 3)
    start = start+num
    num = Nh*Ns*Nc*Nh
    gammac = np.array(individual[start:start+num])
    gammac = gammac.reshape(Nc, Ns, Nh, Nh)
    for i in range(Nc):
        gammac[i] = gammac[i]*CPc[i]
    gammac = gammac.transpose(1, 0, 2, 3)
    
    #Extracting T variables
    start = start+num
    num = Nh*Nc*Ns
    Th = np.array(individual[start:start+num])
    Th = Th.reshape(Nh, Ns, Nc)
    for i in range(Nh):
        Th[i] = Th[i]*(Tih[i]-Toh[i]) + Toh[i]
    Th = Th.transpose(1, 0, 2)
    start = start+num
    Tc = np.array(individual[start:start+num])
    Tc = Tc.reshape(Nc, Ns, Nh)
    for i in range(Nc):
        Tc[i] = Tc[i]*(Toc[i]-Tic[i]) + Tic[i]
    Tc = Tc.transpose(1, 0, 2)
    
    #Extracting beta variables
    start = start+num
    betah = np.array(individual[start:start+num])
    betah = betah.reshape(Nh, Ns, Nc)
    for i in range(Nh):
        betah[i] = CPh[i]*betah[i]
    betah = betah.transpose(1, 0, 2)
    start = start+num
    betac = np.array(individual[start:start+num])
    betac = betac.reshape(Nc, Ns, Nh)
    for i in range(Nc):
        betac[i] = CPc[i]*betac[i]
    betac = betac.transpose(1, 0, 2)
    
    #Extracting t variables
    start = start+num
    th = np.array(individual[start:start+num])
    th = th.reshape(Nh, Ns, Nc)
    for i in range(Nh):
        th[i] = th[i]*(Tih[i]-Toh[i]) + Toh[i]
    th = th.transpose(1, 0, 2)
    start = start+num
    tc = np.array(individual[start:start+num])
    tc = tc.reshape(Nc, Ns, Nh)
    for i in range(Nc):
        tc[i] = tc[i]*(Tih[i]-Toh[i]) + Toh[i]
    tc = tc.transpose(1, 0, 2)
    
    #Extracting q variables
    start = start+num
    qh = np.array(individual[start:start+num])
    qh = qh.reshape(Nh, Ns, Nc)
    for i in range(Nh):
        qh[i] = qh[i]*Hh[i]
    qh = qh.transpose(1, 0, 2)
    start = start+num
    qc = np.array(individual[start:start+num])
    qc = qc.reshape(Nc, Ns, Nh)
    for i in range(Nc):
        qc[i] = qc[i]*Hc[i]
    qc = qc.transpose(1, 0, 2)
    
    #Extracting Tbar variables
    start = start+num
    num = Ns*Nh
    Tbarh = np.array(individual[start:start+num])
    Tbarh = Tbarh.reshape(Nh, Ns)
    for i in range(Nh):
        Tbarh[i] = Tbarh[i]*(Tih[i]-Toh[i]) + Toh[i]
    Tbarh = Tbarh.transpose(1, 0)
    start = start+num
    num = Ns*Nc
    Tbarc = np.array(individual[start:start+num])
    Tbarc = Tbarc.reshape(Nc, Ns)
    for i in range(Nc):
        Tbarc[i] = Tbarc[i]*(Toc[i]-Tic[i]) + Tic[i]
    Tbarc = Tbarc.transpose(1, 0)
    
    #Extracting Q variables
    start = start+num
    num = Nh
    Qh = np.array(individual[start:start+num])
    Qh = Qh*Hh
    start = start+num
    num = Nc
    Qc = np.array(individual[start:start+num])
    Qc = Qc*Hc
    
    return  [alphah, alphac, gammah, gammac, Th, Tc, betah, betac, Tbarh, Tbarc, th, tc, qh, qc, Qh, Qc]




def Beta_Expression(alpha, gamma, beta):
    
    alphah, alphac = alpha
    gammah, gammac = gamma
    betah, betac = beta
    
    mh = np.matmul(gammah.transpose(0, 1, 3, 2), np.ones((Nc, 1))) + alphah
    mh = betah - mh
    
    mc = np.matmul(gammac.transpose(0, 1, 3, 2), np.ones((Nh, 1))) + alphac
    mc = betac - mc
    
    return np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)


def Tbar_Expression(T, gamma, Tbar):
    
    Th, Tc = T
    gammah, gammac = gamma
    Tbarh, Tbarc = Tbar
    
    val = np.sum(abs(Tbarh[0]-Tih)>eps) + np.sum(abs(Tbarc[Ns-1]-Tic)>eps)
    
    for k in range(Ns-1):
        mh = np.matmul(np.ones((1, Nc)), (np.diagonal(gammah[k], axis1 = 1, axis2 = 2).reshape(Nh, Nc, 1)))
        mh = mh*Tbarh[k+1] - Th[k]*np.diagonal(gammah[k], axis1 = 1, axis2 = 2)
        
        mc = np.matmul(np.ones((1, Nh)), (np.diagonal(gammac[Ns-1-k], axis1 = 1, axis2 = 2).reshape(Nc, Nh, 1)))
        mc = mc*Tbarc[Ns-2-k] - Tc[Ns-1-k]*np.diagonal(gammac[Ns-1-k], axis1 = 1, axis2 = 2)
        
        val = val + np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)
        
    return val


def t_Expression(T, gamma, Tbar, alpha, beta, t):
    
    Th, Tc = T
    gammah, gammac = gamma
    Tbarh, Tbarc = Tbar
    alphah, alphac = alpha
    betah, betac = beta
    th, tc = t
    val = 0
    
    for k in range(Ns):
        for i in range(Nh):
            dgG = gammah[k][i]
            np.fill_diagonal(dgG, 0)
            mh = (betah[k][i]*th[k][i]) - (Tbarh[k][i]*alphah[k][i])
            mh = mh - (np.matmul(Th[k][i].reshape(1, Nc), dgG)).reshape(Nc)
            val = val + np.sum(abs(mh)>eps)
    for k in range(Ns):
        for j in range(Nc):
            dgG = gammac[k][j]
            np.fill_diagonal(dgG, 0)
            mc = (betac[k][j]*tc[k][j]) - (alphac[k][j]*Tbarc[k][j])
            mc = mc - (np.matmul(Tc[k][j].reshape(1, Nh), dgG)).reshape(Nh)
            val = val + np.sum(abs(mc)>eps)
    return val




def q_Expression(q, beta, t, T):
    qh, qc = q
    betah, betac = beta
    th, tc = t
    Th, Tc = T
    
    mh = qh - betah*(th-Th)
    mc = qc - betac*(Tc-tc)
    
    return np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)



def heat_Conservation(q):
    qh, qc = q
    m = qh - qc.transpose(0, 2, 1)
    
    return np.sum(abs(m)>eps)


def CP_ConserveHE(beta, gamma):
    betah, betac = beta
    gammah, gammac = gamma
    
    mh = betah - np.matmul(gammah, np.ones((Nc, 1))).reshape(Ns, Nh, Nc)
    mc = betac - np.matmul(gammac, np.ones((Nh, 1))).reshape(Ns, Nc, Nh)
    
    return np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)


def CP_ConserveStreams(alpha):
    alphah, alphac = alpha
    val = 0
    for k in range(Ns):
        mh = np.matmul(alphah[k], np.ones((Nc, 1))).reshape(Nh) - CPh
        mc = np.matmul(alphac[k], np.ones((Nh, 1))).reshape(Nc) - CPc
        val = val + np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)
        
    return val


def NetHeat_Reqd(q, Q):
    qh, qc = q
    Qh, Qc = Q
    
    mh = np.matmul(qh.transpose(1, 0, 2), np.ones((Nc, 1)))
    mh = np.matmul(np.ones((1, Ns)), mh).reshape(Nh)
    mh = mh + Qh - Hh
    
    mc = np.matmul(qc.transpose(1, 0, 2), np.ones((Nh, 1)))
    mc = np.matmul(np.ones((1, Ns)), mc).reshape(Nc)
    mc = mc + Qc - Hc
    
    return np.sum(abs(mh)>eps) + np.sum(abs(mc)>eps)





def HeatInequality(q, Q):
    qh, qc = q
    Qh, Qc = Q
    
    return np.sum(qh<0) + np.sum(qc<0) + np.sum(Qh<0) + np.sum(Qc<0)


def SplitInequality(alpha, gamma):
    alphah, alphac = alpha
    gammah, gammac = gamma
    
    return np.sum(gammah<0) + np.sum(gammac<0) + np.sum(alphah<0) + np.sum(alphac<0)


def dTminInequality(T, t, q):
    Th, Tc = T
    th, tc = t
    qh, qc = q
    
    return np.sum((th-Tc.transpose(0, 2, 1)-dTmin)*qh < 0) + np.sum((Th-tc.transpose(0, 2, 1)-dTmin)*qh < 0)



def evalOneMax(individual):
    
    mats = convertToMatrices(individual)
    alphah, alphac, gammah, gammac, Th, Tc, betah, betac, Tbarh, Tbarc, th, tc, qh, qc, Qh, Qc = mats
    
    #Optimization function
    return (Beta_Expression([alphah, alphac], [gammah, gammac], [betah, betac])
          + Tbar_Expression([Th, Tc], [gammah, gammac], [Tbarh, Tbarc])
          + t_Expression([Th, Tc], [gammah, gammac], [Tbarh, Tbarc], [alphah, alphac], [betah, betac], [th, tc])
          + q_Expression([qh, qc], [betah, betac], [th, tc], [Th, Tc])
          + heat_Conservation([qh, qc])
          + CP_ConserveHE([betah, betac], [gammah, gammac])
          + CP_ConserveStreams([alphah, alphac])
          + NetHeat_Reqd([qh, qc], [Qh, Qc])
          + HeatInequality([qh, qc], [Qh, Qc])
          + SplitInequality([alphah, alphac], [gammah, gammac])
          + dTminInequality([Th, Tc], [th, tc], [qh, qc])),



creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
# Attribute generator 
toolbox.register("variables", random.random)
# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.variables, num_var)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)



def main():
    pop = toolbox.population(n=5000)
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    CXPB, MUTPB = 0.5, 0.2
    fits = [ind.fitness.values[0] for ind in pop]

    # Variable keeping track of the number of generations
    g = 0
    
    # Begin the evolution
    while g < 150:
        # A new generation
        g = g + 1
        print("-- Generation %i --" % g)
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        pop[:] = offspring
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
    
    if(g == 150):
        mats = convertToMatrices(tools.selBest(pop, 1)[0])
        alphah, alphac, gammah, gammac, Th, Tc, betah, betac, Tbarh, Tbarc, th, tc, qh, qc, Qh, Qc = mats
        
        print(alphah)
        print(alphac)
        
        print(gammah)
        print(gammac)
        
        print(qh)
        print(qc)
        


if __name__ == '__main__':
    main()
