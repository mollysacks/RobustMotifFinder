import os, json
import pandas as pd
from collections import OrderedDict, defaultdict
from tqdm import tqdm
from utils import *

seq_file = '../../data/proteins.json' # input sequence file
outdir = '../../result/' # output directory

# create output directory
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# load sequences
seqs = json.load(open(seq_file, 'r'))
len_seqs = len(seqs)

# arguments = (alphabets, num_motifs, num_seqs, num_trials, run_brute_force)
args = OrderedDict({
    'exp1': ('ALV', 6, 10000, 100, True),
    'exp2': ('ACDEFGHIKLMNPQRSTVWY', 4, 5000, 100, False),
    'exp3': ('ACDEFGHIKLMNPQRSTVWY', 6, 5000, 100, False),
    'exp4': ('ACDEFGHIKLMNPQRSTVWY', 8, 5000, 100, False),
    'exp5': ('ALV', 6, 5000, 100, False),
    'exp6': ('ALV', 6, 10000, 100, False),
    'exp7': ('ALV', 6, 20000, 100, False),
    'exp8': ('ACDEFGHIKLMNPQRSTVWY', 6, len_seqs, 1, False) # all samples
})

# experiments
for exp, arg in args.items():
    # mkdir
    if not os.path.isdir('{}/{}'.format(outdir, exp)):
        os.mkdir('{}/{}'.format(outdir, exp))

    # arguments
    alphabets, num_motifs, num_seqs, num_trials, run_brute_force = arg
    if run_brute_force:
        algos = ['BruteForce', 'SimpleGreedySearch', 'DoubledGreedySearch', 'GreedySearchWithOverlap', 'ScoringGreedySearch', 'GibbsSampling']
    else:
        algos = ['SimpleGreedySearch', 'DoubledGreedySearch', 'GreedySearchWithOverlap', 'ScoringGreedySearch', 'GibbsSampling']
    
    # RobustMotifFinder
    Main = RobustMotifFinder(alphabets=alphabets, num_motifs=num_motifs)

    # recording
    seq_idx_list = list()
    r_df = pd.DataFrame(columns=algos)
    t_df = pd.DataFrame(columns=['BuildMatrix',] + algos)
    motif_dict = defaultdict(list)
    
    # trials
    for _ in tqdm(range(num_trials)):
        # main function
        tmp_seq_idx = random.sample(list(range(len_seqs)), num_seqs)
        tmp_seqs = [seqs[i] for i in tmp_seq_idx]
        r, m, t = Main(tmp_seqs, algos=algos)

        # recording
        seq_idx_list.append(tmp_seq_idx)
        r_df.loc[r_df.shape[0]] = r
        t_df.loc[t_df.shape[0]] = t
        for algo in algos:
            motif_dict[algo].append(m[algo])

    # saving
    r_df.to_csv('{}/{}/robustness.csv'.format(outdir, exp), index=False)
    t_df.to_csv('{}/{}/time.csv'.format(outdir, exp), index=False)
    json.dump(seq_idx_list, open('{}/{}/seq_idx.json'.format(outdir, exp), 'w'))
    json.dump(motif_dict, open('{}/{}/motif.json'.format(outdir, exp), 'w'))

    # log result
    with open('{}/{}/result.txt'.format(outdir, exp), 'w') as f:
        # run time
        f.write('Run Time\n')
        build_matrix_time = t_df['BuildMatrix'].mean()
        for algo in algos:
            run_time = build_matrix_time + t_df[algo].mean()
            f.write('{}: {}\n'.format(algo, run_time))
        
        # robustness
        f.write('\nRobustness\n')
        for algo in algos:
            robustness = r_df[algo].mean()
            f.write('{}: {}\n'.format(algo, robustness))

    print('{} completed'.format(exp))