from RuntimeDecoratorPattern import ReportRuntime
import itertools

# read in proteins
with open('proteins.json') as f:
    data = f.read()

data = data[2:-2]
S = data.split('", "')

S = S[:100]  # subset for testing
k, n = 9, 6

# define functions
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

@ReportRuntime
def GenerateAllPossibleMotifs():
    # amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    amino_acids = ['A', 'C', 'D']
    all_motifs = list(itertools.product(amino_acids, amino_acids))  # Cartesian product to generate all possible motifs
    all_motifs = [(9, *m) for m in all_motifs]
    return all_motifs

@ReportRuntime
def GenerateMotifSets(motifs, m):
    motif_sets = list(itertools.combinations(motifs, m))
    return motif_sets

@ReportRuntime
def FindBestMotifs(motif_sets):
    best_r, best_motifs = 0, []
    for motifs in motif_sets:
        r = Robustness(S, motifs)
        if r > best_r:
            best_r, best_motifs = r, [motifs]
        elif r == best_r:
            best_motifs.append(motifs)
    return best_r, best_motifs


all_motifs, dt = GenerateAllPossibleMotifs()
print(f'GenerateAllPossibleMotifs: {dt}')

all_motif_sets, dt = GenerateMotifSets(all_motifs, 2)  # subset by using 2 instead of 1
print(f'GenerateMotifSets: {dt}')

(robustness, best_motifs), dt = FindBestMotifs(all_motif_sets)
print(f'FindBestMotifs: {dt}')
