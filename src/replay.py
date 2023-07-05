# Load data
from Bio import SeqIO
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import ShuffleSplit

#import os
#os.chdir("..")

splitter = ShuffleSplit(n_splits=3, test_size=0.33)

def parse_fasta(path):
    return [(r.id, r.seq[:20]) for r in  SeqIO.parse(path, "fasta")]

positives = parse_fasta("tmp/re_positives.faa")
negatives = parse_fasta("tmp/re_negatives.faa")

N = len(positives) + len(negatives)
X = np.zeros((N, 0))
Y = np.concatenate([np.ones((len(positives), 1)), np.zeros((len(negatives), 1))])

Z_SCORES = {
    "A": ( 0.07, -1.73,  0.09),
    "V": (-2.69, -2.53, -1.29),
    "L": (-4.19, -1.03, -0.98),
    "I": (-4.44, -1.68, -1.03),
    "P": (-1.22,  0.88,  2.23),
    "F": (-4.92,  1.30,  0.45),
    "W": (-4.74,  3.65,  0.85),
    "M": (-2.49, -0.27, -0.41),
    "K": ( 2.84,  1.41, -3.14),
    "R": ( 2.88,  2.52, -3.44),
    "H": ( 2.41,  1.74,  1.11),
    "G": ( 2.23, -5.36,  0.30),
    "S": ( 1.96, -1.63,  0.57),
    "T": ( 0.92, -2.09, -1.40),
    "C": ( 0.71, -0.97,  4.13),
    "Y": (-1.39,  2.32,  0.01),
    "N": ( 3.22,  1.45,  0.84),
    "Q": ( 2.18,  0.53, -1.14),
    "D": ( 3.64,  1.13,  2.36),
    "E": ( 3.08,  0.39, -0.07)
}

def test_model(model, X, y, nruns):
    total = []
    acc = []
    for i in range(nruns):
        acc.clear()
        for (traini, testi) in splitter.split(X):
            model.fit(X[traini], y[traini])
            acc.append(model.score(X[testi], y[testi]))
        total.append(sum(acc) / len(acc))
    return total

# Models:
#%% Random model
acc_random = []
for i in range(350):
    for (traini, testi) in splitter.split(X):
        model = DummyClassifier(strategy='uniform').fit(X[traini], Y[traini])
        acc_random.append(model.score(X[testi], Y.ravel()[testi]))
print(sum(acc_random) / len(acc_random))

# Baseline model: Most common class
acc_baseline = []
for i in range(350):
    for (traini, testi) in splitter.split(X):
        model = DummyClassifier(strategy='most_frequent').fit(X[traini], Y[traini])
        acc_baseline.append(model.score(X[testi], Y.ravel()[testi]))
print(sum(acc_baseline) / len(acc_baseline))

#%% Logistic regression of the amino acid regardless of position
# Load data
AA_TO_INDEX = {char: i for (i, char) in enumerate("ACDEFGHIKLMNPQRSTVWY")}
INDEX_TO_AA = {v:k for (k,v) in AA_TO_INDEX.items()}

def count_aas(string):
    y = np.zeros(20)
    for i in string:
        y[AA_TO_INDEX[i]] += 1
    return y

X = np.concatenate([
    np.array([count_aas(str(i[1])) for i in positives]),
    np.array([count_aas(str(i[1])) for i in negatives])
])

#%%
acc_logisitic = []
for i in range(350):
    for (traini, testi) in splitter.split(X):
        model_log = LogisticRegression().fit(X[traini], Y.ravel()[traini])
        acc_logisitic.append(model_log.score(X[testi], Y.ravel()[testi]))
print(sum(acc_logisitic) / len(acc_logisitic))

#%% ACC models
def calc_autocovariance(s, j1, j2, L):
    return sum(
        (Z_SCORES[i][j1] * Z_SCORES[j][j2])
        for (i, j) in zip(s, s[L:])
    ) / (len(s) - L)
# %%
acc_022 = []
acc_X = np.array([[calc_autocovariance(i[1], 0, 2, 2)] for i in positives + negatives])
for i in range(350):
    for (traini, testi) in splitter.split(X):
        model = LogisticRegression().fit(acc_X[traini], Y.ravel()[traini])
        acc_022.append(model.score(acc_X[testi], Y.ravel()[testi]))
print(sum(acc_022) / len(acc_022))
# %%
# More ACCs
acc_accs = []
for j1 in range(3):
    for j2 in range(3):
        for gap in range(1, 4):
            this_acc = []
            acc_accs.append(this_acc)
            acc_X = np.array([[calc_autocovariance(i[1], j1, j2, gap)] for i in positives + negatives])
            for i in range(350):
                for (traini, testi) in splitter.split(X):
                    model = LogisticRegression().fit(acc_X[traini], Y.ravel()[traini])
                    this_acc.append(model.score(acc_X[testi], Y.ravel()[testi]))

max([sum(i)/len(i) for i in acc_accs])

# Ideas:
# Logistic regression per position? Or too obscure?
# %%
def into_aa_matrix(string):
    y = np.zeros(20 * len(string))
    for (i, s) in enumerate(string):
        y[20*i + AA_TO_INDEX[s]] = 1
    return y

X = np.concatenate([
    np.array([into_aa_matrix(str(i[1])) for i in positives]),
    np.array([into_aa_matrix(str(i[1])) for i in negatives])
])

acc_logisitic_2 = []
for i in range(350):
    for (traini, testi) in splitter.split(X):
        model = LogisticRegression().fit(X[traini], Y.ravel()[traini])
        acc_logisitic_2.append(model.score(X[testi], Y.ravel()[testi]))
print(sum(acc_logisitic_2) / len(acc_logisitic_2))
# %%

# Optimal given acc_logistic_2
# Optimal2 MPYYSNSKEKETHSKKNERD
#          ||| | |||||||||  |||
# AfpN18   MPYSSESKEKETHSKETERD

x = model.coef_.copy().reshape((20, -1))

for i in x:
    print("".join(INDEX_TO_AA[i] for i in np.where(i == i.max())[0]))

print("".join([INDEX_TO_AA[np.argmax(i)] for i in x]))

# %%
