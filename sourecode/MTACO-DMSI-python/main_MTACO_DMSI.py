from scipy.special import comb
import time
from MTACO import *
from ScoreFuns import *



clsThreshold = np.array([0.95, 0.85, 0.75, 0.7, 0.65, 0.6, 0.55, 0.53, 0.51, 0.5])
M = np.size(clsThreshold)
Ac = np.zeros(M, dtype=int)
TP = np.zeros(M, dtype=int)
TN = np.zeros(M, dtype=int)
FN = np.zeros(M, dtype=int)
FP = np.zeros(M, dtype=int)
power = 0
power2 = 0
power3 = np.zeros(M, dtype=int)
mdrPower = np.zeros(M, dtype=int)
treePower = np.zeros(M, dtype=int)
startIndex = 1
endIndex = 10
Evaluation_Times = np.zeros(endIndex)
TIME = np.zeros(endIndex)
Dim_FEs = np.zeros(endIndex)
succEvalutation_Times = np.zeros(endIndex)
succTIME = np.zeros(endIndex)
succ = 0
dataIndex = 1
filepath = "./DME/multiplicative/sample=4000/4_0.05_0.10/"
# filepath = "./NDME/NDME-1/"

for dataSet in range(0, endIndex):
    if dataSet < 9:
        ss = '00'+str(dataSet+1)
    elif dataSet < 99:
        ss = '0' + str(dataSet+1)
    else:
        ss = str(dataSet)
    # data = np.loadtxt(filepath +ss+'.txt', delimiter='\t')
    data = np.loadtxt(filepath+ss+'.txt', skiprows=1, delimiter='\t')
    epi_dim = 4  # the dimension of functional snp combination
    CX = [96, 97, 98, 99] # functional SNP combination
    # parameter setting
    m, n = np.shape(data)
    Dim = n - 1
    p_value0 = max(1e-10, 0.05/comb(Dim, epi_dim))
    #    Search k-snp loci using Harmony search algorithm
    s = 2
    time_start = time.time()
    #  Parameter setting
    antNum = 50
    Esize = 10
    SNPNum = data.shape[1] - 1
    Task_num = 2 * epi_dim - 1
    Max_FEs = 100000 * pow(epi_dim - 1, 2) # Maximum number of evaluations of associations between SNP combinations and disease status
    state = data[:, -1]

    # Main function
    solution, Total_FEs, NC, Epi_Dim_FEs = MTACO(data, state, Task_num, antNum, SNPNum, Esize, Max_FEs, epi_dim, CX)
    time_end = time.time()
    runtime = time_end - time_start
    TIME[dataSet] = runtime
    Evaluation_Times[dataSet] = NC
    Dim_FEs[dataSet] = Epi_Dim_FEs
    #     1st stage power
    flag = 0
    if Array_BcontiansA(CX, solution):
        flag = 1
    if flag > 0:
        power = power + 1
        succ = succ + 1
        succEvalutation_Times[dataSet] = NC
        succTIME[dataSet] = runtime
        print('Index:{:d},dataSet {:3d}: success search time={:f},   successTimes  {:2d}/{:2d}, FEs = {:d}'.format(dataIndex, dataSet+1,
               runtime, succ, dataSet - startIndex + 2, Epi_Dim_FEs))
    else:
        print('Index:{:d},dataSet  {:3d}:  search time={:f} **** fail! ***'.format(dataIndex, dataSet, runtime))
    #     2nd and 3rd stage power
    if flag > 0 and GTest_score(data[:, :-1], data[:,-1]) < p_value0:
        power2 = power2 + 1
        CE, CA, PSI, F_score = MDR_2020_2(data[:, :-1], data[:,-1])
        for m in range(np.size(clsThreshold)):
            if CA > clsThreshold[m]:
                mdrPower[m] = mdrPower[m] + 1
                power3[m] = power3[m] + 1
                TP[m] = TP[m] + 1
            else:
                FN[m] = FN[m] + 1
            TN[m] = TN[m] + 1
    elif flag == 0:
        FN[:] = FN[:] + 1
        for k in range(4):
            X = solution
            L = np.size(X)
            for j in range(L):
                x = X[j, :]
                if GTest_score((data[:,x]), data[:,-1]) > p_value0:
                    TN[:] = TN[:] + 1
                    CE, CA, PSI, F_score = MDR_2020_2((data[:, x]),data[:,-1])
                    for m in range(np.size(clsThreshold)):
                        if CA < clsThreshold[m]:
                            TN[m] = TN[m] + 1
                        else:
                            FN[m] = FN[m] + 1
                        FP[m] = FP[m] + 1



meanFEs = np.mean(Evaluation_Times)
meanTime = np.mean(TIME)
Dim_FEs = np.mean(Dim_FEs)
succFEs = np.mean(succEvalutation_Times)
succTime = np.mean(succTIME)
precision = TP/ (TP + FP + 1e-10)
recall = TP / (TP + FN + 1e-10)
F_score = (2*precision * recall) / (precision + recall+1e-10)
sensibility = TP / (TP + FN+1e-10)
specificity = TN / (TN + FP+1e-10)
print(meanFEs, meanTime, Dim_FEs, succFEs, precision, recall, F_score, sensibility, specificity)