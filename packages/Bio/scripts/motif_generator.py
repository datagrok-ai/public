#!/usr/bin/env python3

import random
from math import sqrt
import argparse
import sys

from typing import List, Tuple

letter_choice_type = List[str]
motif_template_type = List[letter_choice_type]

default_alphabet = 'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y'

def meanrange(mean:int,disp:int) -> int:
    return random.randint(mean - disp, mean + disp)

def generate_modif_template(motif_length:int, alphabet:List[str], max_variants_cluster:int, prob_any:float=0.2) -> motif_template_type:  # Making a template to generate from it some random motifs
    motif_template = []
    for position in range(motif_length):
        # Selecting letters for position i
        if (0 < position < motif_length-1) and (random.random() <= prob_any):
                letters = ['?'] # this stands for any symbol
        else:
            n_variants = random.randrange(max_variants_cluster) + 1
            letters = [ random.choice(alphabet) for i in range(n_variants)]
        motif_template.append(letters)
    return motif_template
    
def generate_motif(template: motif_template_type, alphabet:List[str]) -> str:
    # Sunbtituting the ? in template for any letter
    template_with_any = [ (letters if not '?' in letters else alphabet) for letters in template ]
    return ''.join([ random.choice(letters) for letters in template_with_any ])

def motif_notation(motif_template: motif_template_type) -> str:
    def motif_notation_code(letter_choice:letter_choice_type) -> str:
        if len(letter_choice) == 1:
            return(letter_choice[0])
        else:
            return f"[{''.join(letter_choice)}]"
    
    return ''.join([ motif_notation_code(letter_choice) for letter_choice in motif_template])

def generate_random(n:int, alphabet:List[str]) -> str:
    return ''.join([ random.choice(alphabet) for i in range(n) ])

def make_cliff(motif_template:motif_template_type, alphabet:List[str] , motif:str) -> str:
    # Selecting conservative letter in motif
    pos = random.randrange(len(motif_template))
    while '?' in motif_template[pos]:
        pos = (pos + 1) % len(motif_template) # always will find letters since ends of motif can't be any symbol
    outlier_letters = list(set(alphabet) - set (motif_template[pos]))
    return motif[:pos] + random.choice(outlier_letters) + motif[pos+1:]
    
# ====================================================================================

parser = argparse.ArgumentParser(prog='MotifSequencesGenerator',
                    description='The program generates set of sequences containing sequence motifs for SAR fucntionality testing',
                    epilog='Unitity support: Gennadii Zakharov ')

parser.add_argument("-a", "--alphabet", type=str, default=default_alphabet, help="Alphabet to generate sequences, separated by comma",)
parser.add_argument("-c", "--clusters", type=int, default=1, help="Number of clusters")
parser.add_argument("-s", "--sequences", type=int, default=500, help="Number of sequences in each cluster",)
parser.add_argument("-m,", "--motif", type=int, default=12, help="Average length of motif",)
parser.add_argument("-r,", "--random", type=int, default=4, help="Average length of random sequence parts before and after motif",)
parser.add_argument("-d,", "--dispersion", type=int, default=2, help="Variation of total sequence lengths",)

parser.add_argument("--max-variants-position", type=int, default=3, help="maximum number of different letters in motif position",)
parser.add_argument("--cliff-probability", type=float, default=0.01, help="Probabaility to make activity cliff of a sequence",)
parser.add_argument("--cliff-strength", type=float, default=4.0, help="Strength of cliff",)

args = parser.parse_args()

alphabet:List[str] = args.alphabet.split(',')

print('cluster\tsequence_id\tsequence\tactivity\tis_cliff')

line_number = 0

for n_cluster in range(args.clusters):
    activity_average = random.random() * 10 
    activity_dispersion = random.random()
    
    # Generatin motif template for cluster
    motif_length = meanrange(args.motif, args.dispersion)
    motif_template = generate_modif_template(motif_length, alphabet, args.max_variants_position)
    sys.stderr.write(f"Cluster {n_cluster:2} motif template: {motif_notation(motif_template)}\n")
    
    total_length  = meanrange(args.random * 2, args.dispersion) + motif_length
    prefix_length = meanrange(args.random, args.dispersion//2)
    suffix_length = total_length - motif_length - prefix_length
    
    cliff_made = False
    for n_seq in range(args.sequences):
        line_number +=1
        activity = random.gauss(activity_average, activity_dispersion)
        
        motif  = generate_motif(motif_template, alphabet)
        prefix = generate_random(prefix_length, alphabet)
        suffix = generate_random(suffix_length, alphabet)
        seq = prefix + motif + suffix

        is_cliff = random.random() <= args.cliff_probability
        if is_cliff:
            # Making activity cliff
            cliff_motif = make_cliff(motif_template, alphabet, motif)
            cliff_seq = prefix + cliff_motif + suffix
            # Recalculating activity
            cliff_disp = activity_dispersion * args.cliff_strength * (0.5 + random.random())
            activity = activity_average - cliff_disp
            cliff_activity = activity_average + cliff_disp
            
            sys.stderr.write(f"Cliff for sequence #{line_number:4}, cluster {n_cluster} \n")
            sys.stderr.write(f"{activity_average}\t{motif}\t{activity}\n")
            sys.stderr.write(f"{activity_average}\t{cliff_motif}\t{cliff_activity}\n")
            print(f"{n_cluster}\tc{n_cluster}_seq{line_number}\t{cliff_seq}\t{cliff_activity:5.2f}\t{is_cliff}")
            line_number +=1
        print(f"{n_cluster}\tc{n_cluster}_seq{line_number}\t{seq}\t{activity:5.2f}\t{is_cliff}")
            
