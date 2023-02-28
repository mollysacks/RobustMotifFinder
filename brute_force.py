import itertools

# read in proteins
with open('proteins.json') as f:
    data = f.read()
f.close()

data = data[2:-2]
S = data.split('", "')

# k is 9
k = 9

# n is 6
n = 6

# define functions

def Robustness(S, Motifs):
    # create list to record number of proteins recognized
    numRecognized = []

    # create matrix (|S| rows, |Motifs| columns)
    # RecognitionMatrix[i][j] = 1 if motif j recognizes protein i
    for s in S:
        RecognitionMatrix = []
        for m in Motifs:
            if Recognizes(s, m):
                RecognitionMatrix.append(1)
            else:
                RecognitionMatrix.append(0)
        numRecognized.append(RecognitionMatrix)

    # Calculate Robustness r
    r = 0
    for i in range(len(Motifs)):
        MotifsWithout_i = [s[:i] + s[i+1:] for s in numRecognized]
        # for each row (protein)
        for j in MotifsWithout_i:
            # if at least one motif recognizes protein, add 1 to robustness
            if sum(j) >= 1:
                r += 1
    return r

def Recognizes(s, motif):
    k, P_n1, P_c2 = motif
    for i in range(len(s) - k + 1):
        if (s[i] == P_n1) and (s[i + k -1] == P_c2):
            return True
    return False



AminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# generate all possible motifs
allMotifs = []
for i in AminoAcids:
    for j in AminoAcids:
        allMotifs.append((9, i, j))

#
allMotifSets = list(itertools.combinations(allMotifs, n))

maxRobustness = 0
bestMotifs = []
for Motifs in allMotifSets:
    print(Motifs)
    r = Robustness(S, Motifs)
    if r > maxRobustness:
        maxRobustness = r
        bestMotifs = [Motifs]
    elif r == maxRobustness:
        bestMotifs.append(Motifs)

print(f'{len(bestMotifs)} motif sets have maximum robustness {r}')
