import itertools, json

class RobustMotifFinder():
    def __init__(self,
                 alphabets='ACDEFGHIKLMNPQRSTVWY', # amino acid alphabets
                 num_motifs=6, # number of motifs
                 len_motif=9, # length of motifs
                 anchors=[1, 8] # 2nd and last positions
                ):
        self.alphabets = list(alphabets)
        self.num_motifs = num_motifs
        self.len_motif = len_motif
        self.anchors = anchors
        self.num_anchors = len(anchors)
    

    ###############
    # Motifs
    ###############
    def GenerateAllPossibleMotifs(self, alphabets, num_anchors):
        all_motifs = list(itertools.product(alphabets, repeat=num_anchors))
        return all_motifs
    
    def GenerateMotifSets(self, motifs, num_motifs):
        motif_sets = list(itertools.combinations(motifs, num_motifs))
        return motif_sets
    

    ###############
    # Scoring
    ###############
    def Recognizes(self, seq, motif):
        for i in range(len(seq) - self.len_motif + 1): # for each sub-seq
            sub_seq = seq[i: i+self.len_motif]
            checkpoints = list(zip(self.anchors, motif)) # every alphabet in the motif
            recognize = False
            for anchor, alphabet in checkpoints: # check every alphabet
                if sub_seq[anchor] != alphabet: # not recognize (alphabet level)
                    recognize = False
                    break
                else: # recognize (alphabet level)
                    recognize = True
            if recognize: # recognize (sub-seq level)
                return 1
        return 0
    
    def GenerateRecognitionMatrix(self, seqs, motifs, weight=1):
        return [[self.Recognizes(s, m) * weight for m in motifs] for s in seqs] # |seqs| rows, |motifs| columns
    
    def GenerateRecognitionScore(self, recognition_matrix):
        num_seqs = len(recognition_matrix)
        num_motifs = len(recognition_matrix[0])
        return [sum([recognition_matrix[i][j] for i in range(num_seqs)]) for j in range(num_motifs)]

    def Robustness(self, recognition_matrix):
        num_seqs = len(recognition_matrix)
        num_motifs = len(recognition_matrix[0])
        return sum([sum([any(s[:j] + s[j + 1:]) for s in recognition_matrix])/num_seqs for j in range(num_motifs)])/num_motifs
    

    ###############
    # Algorithms
    ###############
    def BruteForce(self, seqs, all_motifs):
        motif_sets = self.GenerateMotifSets(all_motifs, self.num_motifs)
        best_r, best_motifs = 0, list()
        for motifs in motif_sets:
            recognition_matrix = self.GenerateRecognitionMatrix(seqs, motifs)
            r = self.Robustness(recognition_matrix)
            if r > best_r:
                best_r, best_motifs = r, [motifs]
            elif r == best_r:
                best_motifs.append(motifs)
        return best_r, best_motifs
    
    
    def SimpleGreedySearch(self, seqs, all_motifs):
        num_seqs, num_all_motifs = len(seqs), len(all_motifs)

        # recognition matrix and score for all motifs
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, all_motifs)
        recognition_scores = self.GenerateRecognitionScore(recognition_matrix)
        
        # greedy search
        motif_idx_list = list()
        for _ in range(self.num_motifs):
            # choose top score
            score = max(recognition_scores)
            motif_idx = recognition_scores.index(score)
            motif_idx_list.append(motif_idx)

            # update recognition matrix and score
            recognized_idx_list = [i for i in range(num_seqs) if recognition_matrix[i][motif_idx]==1]
            for idx in recognized_idx_list:
                recognition_matrix[idx] = [0,] * num_all_motifs
            recognition_scores = self.GenerateRecognitionScore(recognition_matrix)

            # remove selected motifs
            for idx in motif_idx_list:
                recognition_scores[idx] = 0
        
        # robustness
        motifs = [all_motifs[i] for i in motif_idx_list]
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, motifs)
        r = self.Robustness(recognition_matrix)
        
        return r, motifs
    
    
    def DoubleGreedySearch(self, seqs, all_motifs):
        num_seqs, num_all_motifs = len(seqs), len(all_motifs)

        # recognition matrix and score for all motifs
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, all_motifs)
        recognition_scores = self.GenerateRecognitionScore(recognition_matrix)

        # greedy search
        motif_idx_list = list()
        recognized_arr = [0,] * num_seqs # recognition times by current motifs
        for _ in range(self.num_motifs):
            # choose top score
            score = max(recognition_scores)
            motif_idx = recognition_scores.index(score)
            motif_idx_list.append(motif_idx)
            recognized_arr = [recognized_arr[i] + recognition_matrix[i][motif_idx] for i in range(num_seqs)]
            
            # update recognition matrix and score
            drop_array = [True if recognized_arr[i] >= 2 else False for i in range(num_seqs)]
            for i in range(num_seqs):
                for j in range(num_all_motifs):
                    recognition_matrix[i][j] = 0 if drop_array[i] else recognition_matrix[i][j]
            recognition_scores = self.GenerateRecognitionScore(recognition_matrix)
            
            # remove selected motifs
            for idx in motif_idx_list:
                recognition_scores[idx] = 0

        # robustness
        motifs = [all_motifs[i] for i in motif_idx_list]
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, motifs)
        r = self.Robustness(recognition_matrix)

        return r, motifs
    

    def ScoringGreedySearch(self, seqs, all_motifs):
        num_seqs, num_all_motifs = len(seqs), len(all_motifs)

        # recognition matrix and score for all motifs
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, all_motifs)
        recognition_scores = self.GenerateRecognitionScore(recognition_matrix)

        # greedy search
        motif_idx_list = list()
        recognized_arr = [0,] * num_seqs # recognition times by current motifs
        for _ in range(self.num_motifs):
            # choose top score
            score = max(recognition_scores)
            motif_idx = recognition_scores.index(score)
            motif_idx_list.append(motif_idx)
            recognized_arr = [recognized_arr[i] + recognition_matrix[i][motif_idx] for i in range(num_seqs)]

            # update recognition score
            recognition_scores = [0,] * num_all_motifs
            for i in range(num_seqs):
                # assign score
                if recognized_arr[i] >= 2:
                    score = 0
                elif recognized_arr[i] == 1:
                    score = 1
                else:
                    score = self.num_motifs - 1
                # for each motif
                for j in range(num_all_motifs):
                    recognition_scores[j] += score if recognition_matrix[i][j] == 1 else 0

            # remove selected motifs
            for idx in motif_idx_list:
                recognition_scores[idx] = 0
        
        # robustness
        motifs = [all_motifs[i] for i in motif_idx_list]
        recognition_matrix = self.GenerateRecognitionMatrix(seqs, motifs)
        r = self.Robustness(recognition_matrix)

        return r, motifs