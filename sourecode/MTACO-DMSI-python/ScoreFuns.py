import numpy as np
import scipy.stats as states


#  Calculated the frequency of genotype combinations in case-control sample
def GeneTypeCombinationCount(snp_com, state):
    ua = np.unique(state)
    if ua.size > 2:
        print("Class state is not equal to 2 !!!")
    elif min(ua) > 0:
        state = state - min(ua)
    geneType = np.unique(snp_com, axis=0, return_inverse=True)
    lrow = len(geneType[0])  # Types of genotype combinations
    # The first line holds the number of cases sample
    # the second line holds the number of controls sample
    # the third line holds the total number of cases and control sample
    F = np.zeros([3, lrow + 1])
    for j in range(ua.size):
        D = geneType[1]
        D = D[np.where(state == j)]
        for i in range(lrow):
            GeneType_num = np.where(D == i)
            F[j, i] = len(GeneType_num[0])
    F[-1, :] = sum(F[:, :])  # control
    F[0, -1] = sum(F[0, :])  # case
    F[1, -1] = sum(F[1, :])
    F[2, -1] = sum(F[2, :])
    return lrow, F


# G-test score
def GTest_score(snp_com, state):
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    xrow, xcol = snp_com.shape
    G = 0
    Degree = (2 - 1) * (lrow - 1)
    for i in range(2):
        for j in range(lrow):
            O = F[i, j] / xrow
            E = F[i, -1] / xrow * F[-1, j] / xrow
            if F[-1, j] > 0:
                if O > 0:
                    G = G + (xrow * O) * np.log(O / E)
            elif Degree > 1:
                Degree = Degree - 0.5
    G = 2 * G
    P_value = states.chi2.pdf(G, Degree)
    return P_value


# LR-score
def LRScore(snp_com, state):
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    xrow, xcol = snp_com.shape
    G = 0
    Degree = (2 - 1) * (lrow - 1)
    for i in range(2):
        for j in range(lrow):
            O = F[i, j] / xrow
            E = F[i, -1] / xrow * F[-1, j] / xrow
            if F[-1, j] > 0:
                if O > 0:
                    G = G + (xrow * O) * np.log(O / E)
            elif Degree > 1:
                Degree = Degree - 0.5
    G = 2 * G
    LR = 1 / G
    return LR


def My_factorial(e):
    f = 0
    if e > 0:
        e = np.int(e)
        for o in range(e):
            f = f + np.log(o+1)
    return f


# K2-score
def K2score(snp_com, state):
    K2score = 0
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    for i in range(lrow):
        if F[-1, i] > 0:
            y = My_factorial(F[-1, i] + 1)
            r = My_factorial(F[0, i]) + My_factorial(F[1, i])
            K2score = K2score + (r - y)
    K2score = np.abs(K2score)
    return K2score


# Mutual entropy
def MutualEntropy(snp_com, state):
    xrow, xcol = snp_com.shape
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    sCase = F[1, -1]
    P = F[2, 0:-2] / xrow
    Psample = P
    PCC = np.array([F[0, 0:-2] / xrow, F[1, 0:-2] / xrow])
    PCC = PCC[np.where(PCC > 0)]

    s1 = sCase / xrow
    s2 = 1 - s1

    MeY = -(s1 * np.log2(s1) + s2 * np.log2(s2))
    MeX = -sum(np.multiply(Psample, np.log2(Psample)))
    MeXY = -sum(np.multiply(PCC, np.log2(PCC)))
    ME = 1 / (MeX + MeY - MeXY)
    return ME


# Calculated the frequency of genotypes in case-control sample
def GeneTypeCount(snp_com, state):
    xrow, xcol = snp_com.shape
    Case = snp_com[np.where(state == 1)]
    Control = snp_com[np.where(state == 0)]
    Hcase = np.zeros([3, xcol])
    Hcontrol = np.zeros([3, xcol])
    for j in range(xcol):
        for i in range(3):
            geneType_case = np.where(Case[:, j] == i)
            Hcase[i, j] = len(geneType_case[0])

            geneType_control = np.where(Control[:, j] == i)
            Hcontrol[i, j] = len(geneType_control[0])
    return Hcase, Hcontrol


# Joint Entropy of disease genotype combination
def JointEntropy(snp_com, state):
    Hcase, Hcontrol = GeneTypeCount(snp_com, state)
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    A = sum((Hcase - Hcontrol) ** 2)
    SDC = np.sqrt(sum(A)) / (sum(np.abs(F[0, :-1] - F[1, :-1])) + 1e-10)
    Pdisease = F[1, 0:-1] / F[1, -1]
    Pcontrol = F[0, 0:-1] / F[0, -1]
    # Pdisease = Pdisease[np.where(Pdisease > 0)]
    # Pcontrol = Pcontrol[np.where(Pcontrol > 0)]
    # JE = sum(np.multiply(Pdisease, np.log2(Pdisease)))
    # JC = sum(np.multiply(Pcontrol, np.log2(Pcontrol)))
    JC = 0
    JE = 0
    for i in range(lrow):
        if Pdisease[i] > 0:
            JE = JE + Pdisease[i] * np.log2(Pdisease[i])
        if Pcontrol[i] > 0:
            JC = JC + Pcontrol[i] * np.log2(Pcontrol[i])
    JE_Score = (SDC / JC ** 2)
    return JE_Score


# Jensen-Shannon divergence
def JensenShannon(snp_com, state):
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    Pdisease = F[1, 0:-1] / F[1, -1]
    Pcontrol = F[0, 0:-1] / F[0, -1]
    JS1 = 0
    JS2 = 0
    for i in range(lrow):
        if Pdisease[i] > 0:
            JS1 = JS1 + Pdisease[i] * np.log2(2 * Pdisease[i] / (Pdisease[i] + Pcontrol[i]))
        if Pcontrol[i] > 0:
            JS2 = JS2 + Pcontrol[i] * np.log2(2 * Pcontrol[i] / (Pdisease[i] + Pcontrol[i]))
    JS = 1 / (JS1 + JS2+ 1e-10)
    return JS


# Multi-factor dimensionality reduction
# Reference: William S et al. (2008) Alternative contingency table measures improve
# the power and detection of multi-factor dimensionality reduction
def MDR_2020_2(snp_com, state):
    lrow, F = GeneTypeCombinationCount(snp_com, state)
    xrow, xcol = snp_com.shape
    disease = F[1, 0:-1]
    control = F[0, 0:-1]
    ProCaseToControl = (disease + 0.1) / (control + 0.1)
    sample_num = xrow
    caseNum = sum(state)
    controlNum = sample_num - caseNum
    threshold = caseNum / controlNum
    Hcase = 0
    Hcontrol = 0
    Lcase = 0
    Lcontrol = 0
    for i in range(lrow):
        if ProCaseToControl[i] > threshold:
            Hcase = Hcase + disease[i]
            Hcontrol = Hcontrol + control[i]
        else:
            Lcase = Lcase + disease[i]
            Lcontrol = Lcontrol + control[i]
    TP = Hcase
    FP = Hcontrol
    H = TP + FP
    FN = Lcase
    TN = Lcontrol
    L = TN + FN
    T = TP + TN
    F = FP + FN
    ALL = H + L
    O = np.array([[TP, FP, H], [TN, FN, L], [T, F, ALL]])
    E = np.array([[H * T / ALL, F * H / ALL], [L * T / ALL, F * L / ALL]])
    # PSI (Predictive Summary Index) (Linn and Grunau 2006)
    PSI = TP / (TP + FP) + TN / (TN + FN) - 1
    # Classification Error (CE)
    CE = (FP + FN) / (FP + FN + TP + TN)
    # F-Measure
    beta = 0.5
    F = (beta ** 2 + 1) * (TP / (TP + FP) * (TP / (TP + FN))) / (beta ** 2 * TN / (TN + FN) * TN / (TN + FP))
    # Classify accuracy
    CA = 1 - CE
    return CA, CE, PSI, F


def multi_criteriaEvaluationFuns(snp_com, state):
    K2_score = K2score(snp_com, state)
    ME = MutualEntropy(snp_com, state)
    JE = JointEntropy(snp_com, state)
    LR = LRScore(snp_com, state)
    return K2_score, LR, ME, JE