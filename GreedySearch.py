from RuntimeDecoratorPattern import ReportRuntime
from copy import deepcopy
import itertools
import random
import os

# read in proteins
with open('proteins.json') as f:
    data = f.read()

data = data[2:-2]
S = data.split('", "')

S = S[:1000]  # subset for testing
k, n = 9, 6

# define functions
def WriteOutputToFile(robustness, motif_set, motif_scores, dt, fname):
    if 'output' not in os.listdir():
        os.mkdir('output')

    with open(f'output/{fname}.log', 'w') as f:
        lines = []
        lines += [f'{round(robustness/n, 1), round(dt, 1)}\n']
        lines += [f"{', '.join(['-'.join(str(x) for x in motif) for motif in motif_set])}\n"]
        lines += [f'{", ".join(str(x) for x in [motif_scores])}\n']
        f.writelines(lines)

def Robustness(S, Motifs):
    recognition_matrix = [[Recognizes(s, m) for m in Motifs] for s in S]  # |S| rows, |Motifs| columns
    robustness = sum([sum([any(s[:j] + s[j + 1:]) for s in recognition_matrix]) for j in range(len(Motifs))])
    return robustness

def Recognizes(s, motif):
    k, P_n2, P_c1 = motif
    for i in range(len(s) - k + 1):
        if (s[i+1] == P_n2) and (s[i + k - 1] == P_c1):
            return 1
    return 0

def PeptidesRecognized(S, motif):
    recognition_vector = [Recognizes(s, motif) for s in S]
    return sum(recognition_vector), recognition_vector

def WeightedPeptidesRecognized(S, motif, recognition_table):
    recognition_vector = [Recognizes(s, motif) for s in S]
    for i in range(len(S)):
        if recognition_table[S[i]] == 0:
            recognition_vector[i] *= n
        elif recognition_table[S[i]] >= 2:
            recognition_vector[i] *= 0
    return sum(recognition_vector), recognition_vector

@ReportRuntime
def SimpleGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set, motif_scores = [], []

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        motif_set += [best_motif]
        motif_scores += [best_score]

        all_motifs.remove(best_motif)
        _, recognition_vector = PeptidesRecognized(S, best_motif)
        S = [S[i] for i in range(len(S)) if not recognition_vector[i]]

    robustness = Robustness(tmp_S, motif_set)
    return robustness, motif_set, motif_scores

@ReportRuntime
def DoubledGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set, motif_scores = [], []

    recognized = set()
    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)
        
        motif_set += [best_motif]
        motif_scores += [best_score]

        all_motifs.remove(best_motif)
        _, recognition_vector = PeptidesRecognized(S, best_motif)

        S2 = []
        for i in range(len(S)):
            if not recognition_vector[i]: # not recognized by current motif
                S2.append(S[i])
            elif (recognition_vector[i]) & (S[i] not in recognized): # first recognition
                recognized.add(S[i])
                S2.append(S[i])
            # (recognition_vector[i]) & (S[i] in recognized) -> second recognition -> remove from S2
        S = S2

    robustness = Robustness(tmp_S, motif_set)
    return robustness, motif_set, motif_scores

@ReportRuntime
def GreedySearchWithOverlap(tmp_S, tmp_all_motifs, p=0.95):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set, motif_scores = [], []

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        motif_set += [best_motif]
        motif_scores += [best_score]

        all_motifs.remove(best_motif)
        _, recognition_vector = PeptidesRecognized(S, best_motif)
        recognition_vector = [x*random.choices([0, 1], [1-p, p])[0] for x in recognition_vector]  # 50% drop-out
        S = [S[i] for i in range(len(S)) if not recognition_vector[i]]

    robustness = Robustness(tmp_S, motif_set)
    return robustness, motif_set, motif_scores

@ReportRuntime
def ScoringGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set, motif_scores = [], []

    recognition_table = {s: 0 for s in tmp_S}

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = WeightedPeptidesRecognized(S, motif, recognition_table)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        motif_set += [best_motif]
        motif_scores += [best_score]

        all_motifs.remove(best_motif)
        _, recognition_vector = WeightedPeptidesRecognized(S, best_motif, recognition_table)

        for i in range(len(S)):
            if recognition_vector[i] >= 1:  # recognition counts will cap out at 2, but won't affect scores
                recognition_table[S[i]] += 1

    robustness = Robustness(tmp_S, motif_set)
    return robustness, motif_set, motif_scores


amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
all_motifs = list(itertools.product(amino_acids, amino_acids))  # Cartesian product to generate all possible motifs
all_motifs = [(9, *m) for m in all_motifs]

for f in [SimpleGreedySearch, GreedySearchWithOverlap, DoubledGreedySearch, ScoringGreedySearch]:
    robustness, motif_set, motif_scores, dt = f(S, all_motifs)
    WriteOutputToFile(robustness, motif_set, motif_scores, dt, f.__name__)
    print(f'{f.__name__} DONE')
