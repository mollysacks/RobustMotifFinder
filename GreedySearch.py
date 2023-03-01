import itertools
import time
from copy import deepcopy
import random

# read in proteins
with open('proteins.json') as f:
    data = f.read()

data = data[2:-2]
S = data.split('", "')

S = S[:10000]  # subset for testing
k, n = 9, 6

# define functions
def Robustness(S, Motifs):
    recognition_matrix = [[Recognizes(s, m) for m in Motifs] for s in S]  # |S| rows, |Motifs| columns
    robustness = sum([sum([any(s[:j] + s[j + 1:]) for s in recognition_matrix]) for j in range(len(Motifs))])
    return robustness

def Recognizes(s, motif):
    k, P_n1, P_c2 = motif
    for i in range(len(s) - k + 1):
        if (s[i] == P_n1) and (s[i + k - 1] == P_c2):
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

def SimpleGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set = []

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        print(f'{best_motif}: {best_score}')
        motif_set += [best_motif]

        all_motifs.remove(best_motif)
        _, recognition_vector = PeptidesRecognized(S, best_motif)
        S = [S[i] for i in range(len(S)) if not recognition_vector[i]]

    robustness = Robustness(tmp_S, motif_set)
    print(f'Robustness: {robustness/n}')

def DoubledGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = 2*deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set = []

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        print(f'{best_motif}: {best_score}')
        motif_set += [best_motif]

        all_motifs.remove(best_motif)

        _, recognition_vector = PeptidesRecognized(S, best_motif)
        deleted = set()

        S2 = []
        for i in range(len(S)):
            if (not recognition_vector[i]) or (S[i] in deleted):
                S2.append(S[i])
            else:
                deleted.add(S[i])
        S = S2

    robustness = Robustness(tmp_S, motif_set)
    print(f'Robustness: {robustness/n}')

def GreedySearchWithOverlap(tmp_S, tmp_all_motifs, p):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set = []

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = PeptidesRecognized(S, motif)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        print(f'{best_motif}: {best_score}')
        motif_set += [best_motif]

        all_motifs.remove(best_motif)
        _, recognition_vector = PeptidesRecognized(S, best_motif)

        recognition_vector = [x*random.choices([0, 1], [1-p, p])[0] for x in recognition_vector]  # 50% drop-out
        S = [S[i] for i in range(len(S)) if not recognition_vector[i]]

    robustness = Robustness(tmp_S, motif_set)
    print(f'Robustness: {robustness/n}')

def ScoringGreedySearch(tmp_S, tmp_all_motifs):
    S, all_motifs = deepcopy(tmp_S), deepcopy(tmp_all_motifs)
    motif_set = []

    recognition_table = {s: 0 for s in tmp_S}

    for _ in range(n):
        best_score, best_motif = 0, None
        for motif in all_motifs:
            score, _ = WeightedPeptidesRecognized(S, motif, recognition_table)
            best_score, best_motif = (score, motif) if score > best_score else (best_score, best_motif)

        print(f'{best_motif}: {best_score}')
        motif_set += [best_motif]

        all_motifs.remove(best_motif)
        _, recognition_vector = WeightedPeptidesRecognized(S, best_motif, recognition_table)
        print(recognition_vector)

        for i in range(len(S)):
            if recognition_vector[i] >= 1:  # recognition counts will cap out at 2, but won't affect scores
                recognition_table[S[i]] += 1

    robustness = Robustness(tmp_S, motif_set)
    print(f'Robustness: {robustness/n}')


amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
all_motifs = list(itertools.product(amino_acids, amino_acids))  # Cartesian product to generate all possible motifs
all_motifs = [(9, *m) for m in all_motifs]

SimpleGreedySearch(S, all_motifs)
GreedySearchWithOverlap(S, all_motifs, p=0.95)  # implement hyper-parameter search and restarts (maybe)
DoubledGreedySearch(S, all_motifs)
ScoringGreedySearch(S, all_motifs)
