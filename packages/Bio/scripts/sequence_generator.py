#!/usr/bin/env python3
# name: Sequence generator
# description: Create the model peptides/DNA sequences with peptides data
# language: python
# tags: template, demo
# input: int clusters = 1 [Number of superclusters]
# input: int sequences = 500 [Number of sequences in each supercluster]
# input: int motif_length = 12 [Average length of motif]
# input: int max_variants_position = 3 [Maximum number of different letters in conservative position in motif]
# input: int random_length = 3 [Average length of random sequence parts before and after motif]
# input: int dispersion = 2 [Variation of total sequence length]
# input: string alphabet_key = 'PT' [Sequence alphabet: PT/DNA/RNA/custom. Custom alphabet is a list of values separated by comma]
# input: bool disable_cliffs = False [Disable generation of cliffs]
# input: float cliff_probability = 0.01 [Probability to make activity cliff of a sequence]
# input: float cliff_strength = 4.0 [Strength of cliff]
# output: dataframe sequences

import random
import argparse
import sys

from typing import List, Tuple, Dict, Iterator, Any

alphabet_type = List[str]

letter_choice_type = List[str]
motif_template_type = List[letter_choice_type]

sequence_record_type = Tuple[int, str, float, bool]
sequence_record_cluster_type = Tuple[int, str, str, float, bool]

alphabets: Dict[str, str] = {
    "PT": "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y",
    "DNA": "A,T,G,C",
    "RNA": "A,U,G,C",
}


def mean_range(mean: int, disp: int) -> int:
    return random.randint(max(mean - disp, 0), mean + disp)


def generate_motif_template(
    motif_length: int,
    alphabet: alphabet_type,
    max_variants_cluster: int,
    prob_any: float = 0.2,
) -> motif_template_type:
    motif_template = []
    for position in range(motif_length):
        # Selecting letters for position i
        if (0 < position < motif_length - 1) and (random.random() <= prob_any):
            letters = ["?"]  # this stands for any symbol
        else:
            n_variants = random.randrange(max_variants_cluster) + 1
            letters = [random.choice(alphabet) for i in range(n_variants)]
        motif_template.append(letters)
    return motif_template


def generate_motif(template: motif_template_type, alphabet: alphabet_type) -> str:
    template_with_any = [(letters if not "?" in letters else alphabet) for letters in template]
    return "".join([random.choice(letters) for letters in template_with_any])


def motif_notation(motif_template: motif_template_type) -> str:
    def motif_notation_code(letter_choice: letter_choice_type) -> str:
        if len(letter_choice) == 1:
            return letter_choice[0]
        else:
            return f"[{''.join(letter_choice)}]"

    return "".join([motif_notation_code(letter_choice) for letter_choice in motif_template])


def generate_random(n: int, alphabet: alphabet_type) -> str:
    return "".join([random.choice(alphabet) for i in range(n)])


def make_cliff(motif_template: motif_template_type, alphabet: alphabet_type, motif: str) -> str:
    # Mutate conservative letter in motif
    pos = random.randrange(len(motif_template))
    while "?" in motif_template[pos]:
        pos = (pos + 1) % len(motif_template)  # always will find letters since ends of motif can't be any symbol
    outlier_letters = list(set(alphabet) - set(motif_template[pos]))
    return motif[:pos] + random.choice(outlier_letters) + motif[pos + 1 :]


def generate_cluster(
    n_sequences: int,
    motif_length: int,
    prefix_length: int,
    suffix_length: int,
    max_variants_position: int,
    make_cliffs: bool,
    alphabet: alphabet_type,
    cliff_probability: float,
    cliff_strength: float,
) -> Iterator[sequence_record_type]:
    motif_template = generate_motif_template(motif_length, alphabet, max_variants_position)

    activity_average = random.random() * 10
    activity_dispersion = random.random()
    sys.stderr.write(f"Motif template: {motif_notation(motif_template)}\n")

    for n_seq in range(n_sequences):
        activity = random.gauss(activity_average, activity_dispersion)

        motif = generate_motif(motif_template, alphabet)
        prefix = generate_random(prefix_length, alphabet)
        suffix = generate_random(suffix_length, alphabet)
        seq = prefix + motif + suffix

        is_cliff = make_cliffs and (random.random() <= cliff_probability)
        sequence_record: sequence_record_type = (n_seq, seq, activity, is_cliff)
        yield sequence_record

        if is_cliff:
            # Making activity cliff
            cliff_motif = make_cliff(motif_template, alphabet, motif)
            cliff_seq = prefix + cliff_motif + suffix
            # Recalculating activity
            cliff_disp = activity_dispersion * cliff_strength * (0.5 + random.random())
            activity = activity_average - cliff_disp
            cliff_activity = activity_average + cliff_disp

            # sys.stderr.write(f"Cliff for sequence #{line_number:4}, cluster {n_cluster} \n")
            # sys.stderr.write(f"{activity_average}\t{motif}\t{activity}\n")
            # sys.stderr.write(f"{activity_average}\t{cliff_motif}\t{cliff_activity}\n")
            n_seq += 1
            sequence_record = (n_seq, cliff_seq, cliff_activity, is_cliff)
            yield sequence_record


def generate_sequences(
    n_clusters: int,
    n_sequences: int,
    average_motif_length: int,
    max_variants_position: int,
    average_random_length: int,
    dispersion: int,
    alphabet: alphabet_type,
    make_cliffs: bool,
    cliff_probability: float,
    cliff_strength: float,
) -> Tuple[List[str], List[sequence_record_cluster_type]]:
    headers: List[str] = ["cluster", "sequence_id", "sequence", "activity", "is_cliff"]
    sequences: List[sequence_record_cluster_type] = []

    for n_cluster in range(n_clusters):
        motif_length = mean_range(average_motif_length, dispersion)

        # sys.stderr.write(f"Cluster {n_cluster:2} motif template: {motif_notation(motif_template)}\n")
        total_length = mean_range(average_random_length * 2, args.dispersion) + motif_length
        prefix_length = mean_range(average_random_length, args.dispersion // 2)
        suffix_length = total_length - motif_length - prefix_length
        sys.stderr.write(f"Generating sequences for cluster {n_cluster}\n")
        for n_seq, seq, activity, is_cliff in generate_cluster(
            n_sequences,
            motif_length,
            prefix_length,
            suffix_length,
            max_variants_position,
            make_cliffs,
            alphabet,
            cliff_probability,
            cliff_strength,
        ):
            sequences.append((n_cluster, f"c{n_cluster}_s{n_seq}", seq, activity, is_cliff))
    return headers, sequences


def parse_command_line_args() -> Any:
    parser = argparse.ArgumentParser(
        prog="MotifSequencesGenerator",
        description="The program generates set of sequences containing sequence motifs "
        "for SAR fucntionality testing",
        epilog="Utility support: Gennadii Zakharov",
    )

    parser.add_argument("-c", "--clusters", type=int, default=1, help="Number of superclusters")
    parser.add_argument(
        "-s",
        "--sequences",
        type=int,
        default=500,
        help="Number of sequences in each supercluster",
    )
    parser.add_argument("-m,", "--motif-length", type=int, default=12, help="Average length of motif")

    parser.add_argument(
        "-r,",
        "--random-length",
        type=int,
        default=3,
        help="Average length of random sequence parts before and after motif",
    )
    parser.add_argument(
        "-d,",
        "--dispersion",
        type=int,
        default=2,
        help="Variation of total sequence length",
    )

    available_alphabets = ",".join(list(alphabets.keys()) + ["custom"])
    parser.add_argument(
        "--alphabet",
        type=str,
        default=list(alphabets.keys())[0],
        help=f"Sequence alphabet: {available_alphabets}. Custom alphabet is a list of values separated " f"by comma",
    )
    parser.add_argument(
        "--max-variants-position",
        type=int,
        default=3,
        help="Maximum number of different letters in conservative position in motif",
    )
    parser.add_argument(
        "--cliff-probability",
        type=float,
        default=0.01,
        help="Probability to make activity cliff of a sequence",
    )
    parser.add_argument(
        "--cliff-strength",
        type=float,
        default=4.0,
        help="Strength of cliff",
    )
    parser.add_argument(
        "--disable-cliffs",
        type=bool,
        default=False,
        help="Disable generation of cliffs",
    )

    command_line_args = parser.parse_args()

    return command_line_args


# ====================================================================================

grok = "clusters" in globals()

if not grok:
    # We are not in Datagrok - need to parse command line arguments
    args = parse_command_line_args()
    clusters = args.clusters
    sequences = args.sequences
    motif_length = args.motif_length
    max_variants_position = args.max_variants_position
    random_length = args.random_length
    dispersion = args.dispersion
    alphabet_key = args.alphabet
    disable_cliffs = args.disable_cliffs
    cliff_probability = args.cliff_probability
    cliff_strength = args.cliff_strength

alphabet: alphabet_type = alphabets[alphabet_key].split(",") if alphabet_key in alphabets else alphabet_key.split(",")

# Running sequence generator
header, data = generate_sequences(
    clusters,
    sequences,
    motif_length,
    max_variants_position,
    random_length,
    dispersion,
    alphabet,
    not disable_cliffs,
    cliff_probability,
    cliff_strength,
)

if grok:
    # Exporting data to Datagrok as a pandas dataframe
    import pandas as pd

    sequences = pd.DataFrame.from_records(data, columns=header)
else:
    # Writing results to stdout - no need to work with big and heavy Pandas
    import csv

    csv_writer = csv.writer(sys.stdout, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
    csv_writer.writerow(header)
    for line in data:
        csv_writer.writerow(line)
