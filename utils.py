import itertools, time, random

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

    
    def __call__(self, seqs,
                 algos=['BruteForce', 'SimpleGreedySearch', 'DoubledGreedySearch', 'ScoringGreedySearch', 'GibbsSampling']):
        # recognition matrix
        self.seqs = seqs
        self.all_motifs = self.GenerateAllPossibleMotifs(self.alphabets, self.num_anchors)
        self.num_seqs = len(self.seqs)
        self.num_all_motifs = len(self.all_motifs)
        self.recognition_matrix = self.GenerateRecognitionMatrix(self.seqs, self.all_motifs)

        # results
        run_time_dict = dict()
        robustness_dict = dict()
        motifs_dict = dict()

        if 'BruteForce' in algos:
            ts = time.time()
            robustness_dict['BruteForce'], motifs_dict['BruteForce'] = self.BruteForce()
            te = time.time()
            run_time_dict['BruteForce'] = te - ts
        
        if 'SimpleGreedySearch' in algos:
            ts = time.time()
            robustness_dict['SimpleGreedySearch'], motifs_dict['SimpleGreedySearch'] = self.SimpleGreedySearch()
            te = time.time()
            run_time_dict['SimpleGreedySearch'] = te - ts
        
        if 'DoubledGreedySearch' in algos:
            ts = time.time()
            robustness_dict['DoubledGreedySearch'], motifs_dict['DoubledGreedySearch'] = self.DoubledGreedySearch()
            te = time.time()
            run_time_dict['DoubledGreedySearch'] = te - ts
        
        if 'ScoringGreedySearch' in algos:
            ts = time.time()
            robustness_dict['ScoringGreedySearch'], motifs_dict['ScoringGreedySearch'] = self.ScoringGreedySearch()
            te = time.time()
            run_time_dict['ScoringGreedySearch'] = te - ts

        if 'GibbsSampling' in algos:
            ts = time.time()
            robustness_dict['GibbsSampling'], motifs_dict['GibbsSampling'] = self.GibbsSampling()
            te = time.time()
            run_time_dict['GibbsSampling'] = te - ts
        
        return robustness_dict, motifs_dict, run_time_dict


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
    
    def GenerateRecognitionMatrix(self, seqs, motifs):
        return [[self.Recognizes(s, m) for m in motifs] for s in seqs] # |seqs| rows, |motifs| columns
    
    def Robustness(self, recognition_matrix, motif_idx_list=None, return_sub_scores=False):
        num_seqs = len(recognition_matrix)
        if motif_idx_list is None: # all motifs
            num_motifs = len(recognition_matrix[0])
            scores = [sum([any(s[:j] + s[j + 1:]) for s in recognition_matrix])/num_seqs for j in range(num_motifs)]
        else: # selected motifs
            num_motifs = len(motif_idx_list)
            scores = list()
            for j in range(len(motif_idx_list)):
                idx_list = motif_idx_list[:j] + motif_idx_list[j+1:]
                score = sum([any([s[i] for i in idx_list]) for s in recognition_matrix]) / num_seqs
                scores.append(score)
        if return_sub_scores:
            return scores
        else:
            return sum(scores) / num_motifs
    

    ###############
    # Algorithms
    ###############
    def BruteForce(self):
        motif_sets = list(itertools.combinations(list(range(self.num_all_motifs)), self.num_motifs))
        best_score, best_motif_idxs = 0, list()
        for motif_idxs in motif_sets:
            score = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idxs)
            if score > best_score:
                best_score = score
                best_motif_idxs = motif_idxs
        motifs = [self.all_motifs[i] for i in best_motif_idxs]
        return best_score, motifs
    
    
    def SimpleGreedySearch(self):
        # initialization
        seq_candidate_idxs = list(range(self.num_seqs))
        motif_candidate_idxs = list(range(self.num_all_motifs))
        
        # greedy search
        motif_idx_list = list()
        for _ in range(self.num_motifs):
            # choose motif with top score
            best_motif_idx, best_score = 0, -1
            for j in motif_candidate_idxs:
                score = sum([self.recognition_matrix[i][j] for i in seq_candidate_idxs])
                if score > best_score:
                    best_score = score
                    best_motif_idx = j
            motif_idx_list.append(best_motif_idx)

            # update candidates
            motif_candidate_idxs.remove(best_motif_idx)
            seq_candidate_idxs = [i for i in seq_candidate_idxs if self.recognition_matrix[i][best_motif_idx]==0]

        # robustness
        motifs = [self.all_motifs[i] for i in motif_idx_list]
        robustness = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list)
        
        return robustness, motifs

    
    def DoubledGreedySearch(self):
        # initialization
        seq_candidate_idxs = list(range(self.num_seqs))
        motif_candidate_idxs = list(range(self.num_all_motifs))

        # greedy search
        motif_idx_list = list()
        recognized_arr = [0,] * self.num_seqs # recognition times by current motifs
        for _ in range(self.num_motifs):
            # choose motif with top score
            best_motif_idx, best_score = 0, -1
            for j in motif_candidate_idxs:
                score = sum([self.recognition_matrix[i][j] for i in seq_candidate_idxs])
                if score > best_score:
                    best_score = score
                    best_motif_idx = j
            motif_idx_list.append(best_motif_idx)
            recognized_arr = [recognized_arr[i] + self.recognition_matrix[i][best_motif_idx] for i in range(self.num_seqs)]

            # update candidates
            motif_candidate_idxs.remove(best_motif_idx)
            seq_candidate_idxs = [i for i in seq_candidate_idxs if recognized_arr[i]<2]

        # robustness
        motifs = [self.all_motifs[i] for i in motif_idx_list]
        robustness = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list)
        
        return robustness, motifs
    

    def ScoringGreedySearch(self):
        # initialization
        seq_candidate_idxs = list(range(self.num_seqs))
        motif_candidate_idxs = list(range(self.num_all_motifs))

        # greedy search
        motif_idx_list = list()
        recognized_arr = [0,] * self.num_seqs # recognition times by current motifs
        for _ in range(self.num_motifs):
            # choose motif with top score
            best_motif_idx, best_score = 0, -1
            for j in motif_candidate_idxs:
                score = 0
                for i in seq_candidate_idxs:
                    if self.recognition_matrix[i][j] == 0:
                        continue
                    if recognized_arr[i] == 1: # recognized 1 time
                        score += 1
                    elif recognized_arr[i] == 0: # recognized 0 time
                        score += self.num_motifs
                if score > best_score:
                    best_score = score
                    best_motif_idx = j
            motif_idx_list.append(best_motif_idx)
            recognized_arr = [recognized_arr[i] + self.recognition_matrix[i][best_motif_idx] for i in range(self.num_seqs)]

            # update candidates
            motif_candidate_idxs.remove(best_motif_idx)
            seq_candidate_idxs = [i for i in seq_candidate_idxs if recognized_arr[i]<2]
        
        # robustness
        motifs = [self.all_motifs[i] for i in motif_idx_list]
        robustness = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list)
        
        return robustness, motifs
    

    def GibbsSampling(self):
        # randomly sample motifs
        motif_idx_list = random.sample(list(range(self.num_all_motifs)), self.num_motifs)
        scores = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list, return_sub_scores=True)

        num_pass = 0
        while num_pass < self.num_motifs:
            # find the max score: larger score mean the motif is less important
            max_score = max(scores)
            max_idx = scores.index(max_score)

            # find the motif with best robustness
            best_score, best_motif_idx = 0, 0
            used_motif_idx_list = motif_idx_list[:max_idx] + motif_idx_list[max_idx+1:]
            for j in range(self.num_all_motifs): # for each motif
                if j in used_motif_idx_list: # skip used motifs
                    continue
                tmp_motif_idx_list = motif_idx_list[:max_idx] + [j,] + motif_idx_list[max_idx+1:]
                tmp_score = self.Robustness(self.recognition_matrix, motif_idx_list=tmp_motif_idx_list)
                if tmp_score > best_score:
                    best_score = tmp_score
                    best_motif_idx = j
            
            # update
            if best_motif_idx == motif_idx_list[max_idx]: # no replacement
                scores[max_idx] = 0 # for moving to the next max score
                num_pass += 1
            else: # replacement
                motif_idx_list = motif_idx_list[:max_idx] + [best_motif_idx,] + motif_idx_list[max_idx+1:]
                scores = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list, return_sub_scores=True)
                num_pass = 0
            
        # robustness
        motifs = [self.all_motifs[i] for i in motif_idx_list]
        robustness = self.Robustness(self.recognition_matrix, motif_idx_list=motif_idx_list)
        
        return robustness, motifs