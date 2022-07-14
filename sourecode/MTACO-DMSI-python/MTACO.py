import numpy as np
import random as rd
import math
import copy
from ScoreFuns import K2score, JensenShannon

# Functions of MTACO_DMSI
# ------------MTACO : the search program of MTACO-DMSI in the 1st stage.
# ------------Serach : the population  serach program based on ACO
# ------------Transfer : The Knowledge transfer program .
# ------------updatePheromones : The pheromones update program

#   The score functions for evaluating the association in ScoreFuns.py
#   ------------K2_score  : Bayesian network-based K2-score .
#   ------------JS_score  : Jensen–Shannon divergence-based JS-score.
#   ------------Gtest_score : the G-test function


# Populations class
class AntPop(object):
    def __init__(self, SNPComDim, antNum, SNPNum):
        self.SNPComDim = SNPComDim
        self.num_ant = antNum
        self.SNPNum = SNPNum
        self.FEs = 0
        self.SNPs = np.zeros((antNum, SNPComDim), dtype='int') # current population's track
        self.Fitness = np.zeros((antNum, 1), dtype='float')  # current population's fitness
        self.Tau = np.ones((SNPNum, 1), dtype='float')  # current population's pheromenons


# Elite class
class EliteSet(object):
    def __init__(self, SNPComDim, Esize):
        self.SNPComDim = SNPComDim
        self.Esize = Esize
        self.SNPs = np.zeros((Esize, SNPComDim), dtype='int')
        self.Fitness = np.ones((Esize, 1), dtype='float')* 1e10
        self.wrost = Esize-1


def search(SNPNum, Task_num, antNum,  Pop):
    T0 = 0.8
    for taskDIm in range(Task_num):
        SNPComDim = taskDIm+2
        SNPs = np.zeros((antNum, SNPComDim), dtype='int')
        # randomly initialing ant position
        if antNum <= SNPNum:
            Randpos = rd.sample(range(0, SNPNum), antNum)
            SNPs[:, 0] = Randpos
        else:
            Randpos = rd.sample(range(0, SNPNum), SNPNum)
            for i in range(SNPNum, antNum):
                Randpos.append(rd.randint(0, SNPNum-1))
            SNPs[:, 0] = Randpos
        Tau = Pop[taskDIm].Tau
        for anti in range(antNum):
            Temp_tau = copy.deepcopy(Tau)
            Temp_tau[SNPs[anti, 0], 0] = 0
            for snpi in range(1, SNPComDim):
                p = Temp_tau
                p = p / (sum(p))
                pcum = np.cumsum(p)
                mRate = rd.random()
                # Roulette selection
                if mRate > T0:
                    sel = math.ceil(rd.random() * SNPNum)
                    while(sel in SNPs[anti, :]) or (sel>=SNPNum):
                        sel = math.ceil(rd.random() * SNPNum)
                else:
                    r = rd.random()
                    sel = np.where(pcum >= r)
                    sel = sel[0][0]
                    # while (sel >= SNPNum):
                    #     r = rd.random()
                    #     sel = np.where(pcum >= r)
                    #     sel = sel[0][0]
                SNPs[anti, snpi] = sel
                Temp_tau[sel, 0] = 0
        SNPTemp = np.sort(SNPs, axis=1)
        _, indicates = np.unique(SNPTemp, axis=0, return_index=True)
        Pop[taskDIm].SNPs = SNPs[indicates, :]
    return Pop


def Array_BcontiansA(ArrayA, ArrayB):
    ArrayA = np.sort(ArrayA)
    ArrayB = np.sort(ArrayB)
    a = np.size(np.unique(ArrayB, axis=0), axis=0)
    b = np.size(np.unique(np.append([ArrayA], ArrayB, axis=0), axis=0), axis=0)
    if a == b:
        return True
    else:
        return False


def Transfer(Pop1, Pop2, taskDIm, Esize, numSNP):
    Dim_SNPs= taskDIm +2
    # (k-1)-order transfer to k-order
    SP1 = np.zeros((Esize, Dim_SNPs))
    if taskDIm > 0:
        for SNPi in range(Esize):
            Parent1 = Pop1[taskDIm-1].Elite.SNPs[SNPi, :]
            Parent2 = Pop1[taskDIm].Elite.SNPs[rd.randint(0, Esize-1), :]
            Spring = np.zeros(Parent2.shape[0])
            if np.unique(np.concatenate((Parent1, Parent2))).shape[0] == Dim_SNPs:
                inew = rd.randint(0, numSNP-1)
                while inew in Parent2:
                    inew = rd.randint(0, numSNP - 1)
                Spring[0] = inew
                Spring[1:] = Parent1
            else:
                newi = rd.randint(0, Dim_SNPs-1)
                while Parent2[newi] in Parent1:
                    newi = rd.randint(0, Dim_SNPs-1)
                Spring[0] = Parent2[newi]
                Spring[1:] = Parent1
            SP1[SNPi, :] = Spring
    # (k+1)-order transfer to k-order
    SP2 = np.zeros((Esize, Dim_SNPs))
    if taskDIm < len(Pop1)-1:
        for SNPi in range(Esize):
            Parent1 = Pop1[taskDIm+1].Elite.SNPs[SNPi, :]
            Parent2 = Pop1[taskDIm].Elite.SNPs[rd.randint(0, Esize-1), :]
            if np.unique(np.concatenate((Parent1, Parent2))).shape[0] == Dim_SNPs+1:
                Spring = np.delete(Parent1, rd.randint(0, Dim_SNPs))
            else:
                inew = rd.randint(0, Dim_SNPs-1)
                while Parent1[inew] in Parent2:
                    inew = rd.randint(0, Dim_SNPs-1)
                Spring = np.delete(Parent1, inew)
            SP2[SNPi, :] = Spring
    if taskDIm > 0 and taskDIm >= len(Pop1)-1:
        return np.concatenate((SP1, Pop2[taskDIm].Elite.SNPs))
    elif taskDIm <= 0 and taskDIm < len(Pop1)-1:
        return np.concatenate((SP2, Pop2[taskDIm].Elite.SNPs))
    elif taskDIm > 0 and taskDIm < len(Pop1)-1:
        return np.concatenate((SP1, SP2, Pop2[taskDIm].Elite.SNPs))
    return Pop2[taskDIm].Elite.SNPs


def updatePheromones(Pop, taskDIm):
    rou = 0.8  # Pheromone Evaporation Factor
    Delta1 = 0.2  # Pheromone increment from elite solutions
    Delta2 = 0.1  # Pheromone increment from Current population
    Tau = Pop[taskDIm].Tau
    max_CSscore = np.max(Pop[taskDIm].Elite.Fitness)
    min_CSscore = np.min(Pop[taskDIm].Elite.Fitness)
    Tau = Tau * (1-rou)
    max_score = np.max(Pop[taskDIm].Fitness)
    min_score = np.min(Pop[taskDIm].Fitness)
    # Pheromone increment from Current population
    Ed1 = max_score - min_score + 1e-10
    for i in range(Pop[taskDIm].SNPs.shape[0]):
        lambda1 = (max_score-Pop[taskDIm].Fitness[i, 0]) / Ed1
        Tau[Pop[taskDIm].SNPs[i, :], 0] += lambda1*Delta2
    # % Pheromone increment from elite solutions
    Ed2 = max_CSscore - min_CSscore + 1e-10
    for i in range(Pop[taskDIm].Elite.SNPs.shape[0]):
        lambda2 = (max_CSscore - Pop[taskDIm].Elite.Fitness[i, 0]) / Ed2
        Tau[Pop[taskDIm].Elite.SNPs[i, :], 0] += lambda2 * Delta1
    Pop[taskDIm].Tau = Tau
    return Pop


def MTACO(data, state, Task_num, antNum, SNPNum, Esize, Max_FEs, dim_epi, aim_snp):
    numSNP = data.shape[1]-1
    PopK2 = [AntPop(taskDIm+2, antNum, SNPNum) for taskDIm in range(Task_num)]
    NumIteration = 0
    for taskDIm in range(0, Task_num):
        setattr(PopK2[taskDIm], 'Elite', EliteSet(taskDIm+2, Esize))
    PopJS = [AntPop(taskDIm+2, antNum, SNPNum) for taskDIm in range(Task_num)]
    for taskDIm in range(0, Task_num):
        setattr(PopJS[taskDIm], 'Elite', EliteSet(taskDIm+2, Esize))
    Total_FEs = 0
    while Total_FEs <= Max_FEs:
        PopK2 = search(SNPNum, Task_num, antNum, PopK2)
        PopJS = search(SNPNum, Task_num, antNum, PopJS)
        # Evalue
        for taskDIm in range(Task_num):
            # pop1 used K2_score as Fitness function.
            for SNP in range(PopK2[taskDIm].SNPs.shape[0]):
                if not Array_BcontiansA(PopK2[taskDIm].SNPs[SNP, :], np.sort(PopK2[taskDIm].Elite.SNPs)):
                    SNPcom = data[:, PopK2[taskDIm].SNPs[SNP, :]]
                    PopK2[taskDIm].Fitness[SNP, 0] = K2score(SNPcom, state)
                    PopK2[taskDIm].FEs += 1
                    if PopK2[taskDIm].Fitness[SNP, 0] < PopK2[taskDIm].Elite.Fitness[PopK2[taskDIm].Elite.wrost, 0]:
                        PopK2[taskDIm].Elite.SNPs[PopK2[taskDIm].Elite.wrost, :] = PopK2[taskDIm].SNPs[SNP, :]
                        PopK2[taskDIm].Elite.Fitness[PopK2[taskDIm].Elite.wrost, :] = PopK2[taskDIm].Fitness[SNP, 0]
                        PopK2[taskDIm].Elite.wrost = np.argmax(PopK2[taskDIm].Elite.Fitness)
            # pop2 used JS_score as Fitness function.
            for SNP in range(PopJS[taskDIm].SNPs.shape[0]):
                if not Array_BcontiansA(PopJS[taskDIm].SNPs[SNP, :], np.sort(PopJS[taskDIm].Elite.SNPs)):
                    SNPcom = data[:, PopJS[taskDIm].SNPs[SNP, :]]
                    PopJS[taskDIm].Fitness[SNP, 0] = JensenShannon(SNPcom, state)
                    PopJS[taskDIm].FEs += 1
                    if PopJS[taskDIm].Fitness[SNP, 0] < PopJS[taskDIm].Elite.Fitness[PopJS[taskDIm].Elite.wrost, 0]:
                        PopJS[taskDIm].Elite.SNPs[PopJS[taskDIm].Elite.wrost, :] = PopJS[taskDIm].SNPs[SNP, :]
                        PopJS[taskDIm].Elite.Fitness[PopJS[taskDIm].Elite.wrost, :] = PopJS[taskDIm].Fitness[SNP, 0]
                        PopJS[taskDIm].Elite.wrost = np.argmax(PopJS[taskDIm].Elite.Fitness)
            # knowledge transfer
            if NumIteration > 1:
                new_SNPs_k2 = Transfer(PopK2, PopJS, taskDIm, Esize, numSNP)
                new_SNPs_k2 = new_SNPs_k2.astype(int)
                sortPopK2 = np.sort(PopK2[taskDIm].Elite.SNPs)
                for SNPi in range(new_SNPs_k2.shape[0]):
                    if not Array_BcontiansA(new_SNPs_k2[SNPi, :], sortPopK2):
                        SNPcom = data[:, new_SNPs_k2[SNPi, :]]
                        Fitness = K2score(SNPcom, state)
                        PopK2[taskDIm].FEs += 1
                        wrost = PopK2[taskDIm].Elite.wrost
                        if Fitness < PopK2[taskDIm].Elite.Fitness[wrost, 0]:
                            PopK2[taskDIm].Elite.SNPs[wrost, :] = new_SNPs_k2[SNPi, :]
                            PopK2[taskDIm].Elite.Fitness[wrost, :] = Fitness
                            # update wrost
                            PopK2[taskDIm].Elite.wrost = np.argmax(PopK2[taskDIm].Elite.Fitness)
                new_SNPs_JS = Transfer(PopJS, PopK2, taskDIm, Esize, numSNP)
                new_SNPs_JS = new_SNPs_JS.astype(int)
                sortPopJS = np.sort(PopJS[taskDIm].Elite.SNPs)
                for SNPi in range(new_SNPs_JS.shape[0]):
                    if not Array_BcontiansA(new_SNPs_JS[SNPi, :], sortPopJS):
                        SNPcom = data[:, new_SNPs_JS[SNPi, :]]
                        Fitness = K2score(SNPcom, state)
                        PopJS[taskDIm].FEs += 1
                        wrost = PopJS[taskDIm].Elite.wrost
                        if Fitness < PopJS[taskDIm].Elite.Fitness[wrost, 0]:
                            PopJS[taskDIm].Elite.SNPs[wrost, :] = new_SNPs_JS[SNPi, :]
                            PopJS[taskDIm].Elite.Fitness[wrost, :] = Fitness
                            PopJS[taskDIm].Elite.wrost = np.argmax(PopJS[taskDIm].Elite.Fitness)
            # update pheromenon
            PopK2 = updatePheromones(PopK2, taskDIm)
            PopJS = updatePheromones(PopJS, taskDIm)
        NumIteration += 1
        solution = np.concatenate((PopK2[dim_epi-2].Elite.SNPs, PopJS[dim_epi-2].Elite.SNPs))
        solution = np.unique(np.sort(solution, axis=1), axis=0)
        Total_FEs = 0
        for d in range(Task_num):
            Total_FEs += (PopK2[d].FEs + PopJS[d].FEs)
        print("Loading:------ %d / %d -------" % (Total_FEs, Max_FEs))
        # if detect the functional SNP combination, stop the search progress.
        if Array_BcontiansA(aim_snp, solution):
            break
    Epi_Dim_FEs = (PopK2[dim_epi-2].FEs + PopJS[dim_epi-2].FEs)
    return solution, Total_FEs, NumIteration, Epi_Dim_FEs


# #  Parameter setting
# antNum = 50
# Esize = 10
#
# filepath = './NDME/NDME-5/003.txt'
# # filepath = 'E:/新建文件夹/threshold/sample=4000/4_0.20_0.10/001.txt'
# data = np.loadtxt(filepath, delimiter='\t', skiprows=1).astype(int)
# SNPNum = data.shape[1]-1
# aim_snp = np.array([96, 97, 98, 99])
# dim_epi = aim_snp.shape[0]
# Task_num = 2*dim_epi - 1
# Max_FEs = 100000*pow(dim_epi-1,2)
# state = data[:, -1]
#
# # Main function
# solution, Total_FEs, NC, Epi_Dim_FEs = MTACO(data, state, Task_num, antNum, SNPNum, Esize, Max_FEs, dim_epi, aim_snp)
# print(solution)
# print(Total_FEs)