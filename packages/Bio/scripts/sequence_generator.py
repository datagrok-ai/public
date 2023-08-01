#!/usr/bin/env python3
# name: Sequence generator
# description: Create the model peptides/DNA sequences with peptides data
# language: python
# tags: template, demo
# input: int clusters = 5 [Number of superclusters]
# input: int num_sequences = 50 [Number of sequences in each supercluster]
# input: int motif_length = 12 [Average length of motif]
# input: int max_variants_position = 3 [Maximum number of different letters in conservative position in motif]
# input: int random_length = 3 [Average length of random sequence parts before and after motif]
# input: int dispersion = 2 [Variation of total sequence length]
# input: string alphabet_key = 'PT' [Sequence alphabet: PT/DNA/RNA/custom. Custom alphabet is a list of values separated by comma]
# input: bool disable_cliffs = False [Disable generation of cliffs]
# input: double cliff_probability = 0.01 [Probability to make activity cliff of a sequence]
# input: double cliff_strength = 4.0 [Strength of cliff]
# input: string fasta_separator = '' {nullable: true}
# output: dataframe sequences

"""
The most simple options set running from command line
  python sequence_generator.py -c 4 -s 50 > output_file.tsv
Basic options:
  -с number of clusters
  -s cluster size (number of sequences per cluster)
"""

import random
import argparse
import sys
from enum import Enum

from typing import List, Tuple, Dict, Iterator, Any


# --- Type definitions ---

Letter = str
Alphabet = List[str]

LetterChoice = List[Letter]
MotifTemplate = List[LetterChoice]

Sequence = List[Letter]  # The sequence in a form of list
SequenceSquashed = str  # Sequence, joined together in string form

SequenceRecord = Tuple[int, Sequence, float, bool]
ClusterSequenceRecord = Tuple[int, str, Sequence, float, bool]

# --- constants ---

HelmConnectionMode = Enum("HelmConnectionMode", ["linear", "cyclic", "mixed"])

alphabets: Dict[str, str] = {
    "PT": "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y",
    "DNA": "A,T,G,C",
    "RNA": "A,U,G,C",
    "PT_HELM": "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,dA,dC,dD,dE,dF,dH,dI,dK,dL,dM,dN,dP,dQ,dR,dS,dT,dV,dW,dY,meA,meD,meS,meT,meV,meY,meE,meG,meI,meK,meM,meN,meQ,meC,meR,meW,meF,meH,meL,Nle,Nva,Orn,Iva,aIle,gGlu,Hcy,Hse,Hyp,D-gGlu,D-Nle,D-hPhe,D-Hyp,D-Nva,D-Orn,Pyr,Phe_3Cl,Phe_4Cl,Phe_4NH2,Phg,Ser_tBu,Tyr_Bn,Tza,1Nal,Cha,Lys_Boc,aThr,D-2Nal,D-2Thi,D-aHyp,D-aIle,D-Phg,D-Ser_tBu,Cya,Lys_Me3,Pen,Phe_4Me,Ser_Bn,Tyr_tBu,2Nal,Thi,aHyp,Ala_tBu,hPhe,D-1Nal,D-aThr,D-Cha,D-Pen,D-Phe_4Cl,D-Ser_Bn,Wil,Oic_3aS-7aS,Pip,3Pal,4Pal,Abu,Apm,Chg,Dab,Dap,D-3Pal,D-aMeAbu,D-Chg,D-Cit,D-Dab,D-Pip,D-Tic,Aca,Tic,Aad,Cit,Aze,Ac5c,Aib,D-2Pal,D-Abu,D-Dap,Asu,D-Thz,D-Trp_For,D-Tyr_Et,Lys_Ac,Asp_OMe,Phe_ab-dehydro,Sta_3xi4xi,Tyr_ab-dehydroMe,App,Cap,Cys_SEt,Dsu,pnC,pnG,Pqa,Pro_4Me3OH,Met_O2,Phe_2Me,Phe_34diCl,Phe_4Br,Phe_4I,Phe_4Sdihydroorotamido,Pyl,Ser_PO3H2,Thr_PO3H2,Thz,Trp_Me,Tyr_26diMe,Tyr_3I,Tyr_3NO2,Tyr_Ph4OH,Tyr_SO3H,Val_3OH,xiIle,NMe2Abz,NMebAla,aMePhe,aMePro,aMeTyr_3OH,Bmt,Bmt_E,Cys_Bn,Gla,hHis,His_1Me,Gly_allyl,Gly_cPr,Asp_Ph2NH2,Azi,2Abz,3Abz,4Abz,Ac3c,Ac6c,bAla,D-Bmt,D-Bmt_E,D-hArg,D-Phe_4F,D-Trp_2Me,D-Tyr_Me,D-xiIle,Lys_iPr,Phe_ab-dehydro_3NO2,Sta_3S4S,Bux,Dpm,pnA,pnT,seC,Met_O,nTyr,Oic_3aR-7aS,Oic_3axi-7axi,Phe_2F,Phe_3F,Phe_4F,Phe_4NO2,Phe_bbdiMe,Trp_5OH,Trp_Ome,Tyr_35diI,Tyr_3OH,Tyr_Me,Tyr_PO3H2,xiHyp,xiThr,NMe4Abz,aMeTyr,Aoda,Bpa,Cys_Me,Dip,hArg,His_1Bn,His_3Me,Hyl_5xi,Bip,Abu_23dehydro,D-Dip,Dha,D-hArg_Et2,D-Met_S-O,D-His_1Bn,D-nTyr,D-Phe_4ureido",
}


def mean_range(mean: int, disp: int) -> int:
    return random.randint(max(mean - disp, 0), mean + disp)


def generate_motif_template(
    motif_length: int,
    alphabet: Alphabet,
    max_variants_cluster: int,
    prob_any: float = 0.2,
) -> MotifTemplate:
    motif_template = []
    for position in range(motif_length):
        # Selecting letters for position i
        if (0 < position < motif_length - 1) and (random.random() <= prob_any):
            letters = ["?"]  # this stands for any symbol
        else:
            n_variants = random.randrange(max_variants_cluster) + 1
            letters = list(set((random.choice(alphabet) for i in range(n_variants))))
        motif_template.append(letters)
    return motif_template


def generate_motif(template: MotifTemplate, alphabet: Alphabet) -> Sequence:
    template_with_any = [
        (letters if not "?" in letters else alphabet) for letters in template
    ]
    return [random.choice(letters) for letters in template_with_any]


def motif_notation(motif_template: MotifTemplate) -> str:
    def motif_notation_code(letter_choice: LetterChoice) -> str:
        if len(letter_choice) == 1:
            return letter_choice[0]
        else:
            return f"[{''.join(letter_choice)}]"

    return "".join(
        [motif_notation_code(letter_choice) for letter_choice in motif_template]
    )


def generate_random(n: int, alphabet: Alphabet) -> Sequence:
    return [random.choice(alphabet) for i in range(n)]


def make_cliff(
    motif_template: MotifTemplate, alphabet: Alphabet, motif: Sequence
) -> Sequence:
    # Mutate conservative letter in motif
    motif_len = len(motif_template)
    pos = random.randrange(motif_len)
    while "?" in motif_template[pos]:
        pos = (
            pos + 1
        ) % motif_len  # always will find letters since ends of motif can't be any symbol
    outlier_letters = list(set(alphabet) - set(motif_template[pos]))
    new_letter = random.choice(outlier_letters)
    return (
        motif[:pos]
        + [
            new_letter,
        ]
        + motif[pos + 1 :]
    )


def sequence_to_fasta(sequence: Sequence, separator: str) -> SequenceSquashed:
    return separator.join(sequence)


def sequence_to_helm(
    sequence: Sequence, helm_connection_mode: str = HelmConnectionMode.linear.name
) -> SequenceSquashed:
    def is_cyclic(helm_connection_mode: str) -> bool:
        return helm_connection_mode == HelmConnectionMode.cyclic.name or (
            helm_connection_mode == HelmConnectionMode.mixed.name
            and random.random() < 0.5
        )

    sequence_escaped: Sequence = [
        f"[{letter}]" if len(letter) > 1 else letter for letter in sequence
    ]
    connection_format = ""
    if is_cyclic(helm_connection_mode):
        connection_format = f"PEPTIDE1,PEPTIDE1,{len(sequence_escaped)}:R2-1:R1"
    return f"PEPTIDE1{{{sequence_to_fasta(sequence_escaped,'.')}}}${connection_format}$$$V2.0"


def generate_cluster(
    n_sequences: int,
    motif_length: int,
    prefix_length: int,
    suffix_length: int,
    max_variants_per_position: int,
    make_cliffs: bool,
    alphabet: Alphabet,
    cliff_probability: float,
    cliff_strength: float,
) -> Iterator[SequenceRecord]:
    # Making a motif template
    motif_template = generate_motif_template(
        motif_length, alphabet, max_variants_per_position
    )
    # Setting average and dispersion for activity
    activity_average = random.random() * 10
    activity_dispersion = random.random()
    sys.stderr.write(f"Motif template: {motif_notation(motif_template)}\n")

    for n_seq in range(n_sequences):
        activity = random.gauss(activity_average, activity_dispersion)

        motif = generate_motif(motif_template, alphabet)
        prefix = generate_random(prefix_length, alphabet)
        suffix = generate_random(suffix_length, alphabet)
        seq = prefix + motif + suffix
        sequence_record: SequenceRecord = (n_seq, seq, activity, False)
        yield sequence_record

        is_cliff = make_cliffs and (random.random() <= cliff_probability)
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
    max_variants_per_position: int,
    average_random_length: int,
    dispersion: int,
    alphabet: Alphabet,
    make_cliffs: bool,
    cliff_probability: float,
    cliff_strength: float,
) -> Tuple[List[str], List[ClusterSequenceRecord]]:
    headers: List[str] = ["cluster", "sequence_id", "sequence", "activity", "is_cliff"]
    sequences: List[ClusterSequenceRecord] = []

    for n_cluster in range(n_clusters):
        motif_length = mean_range(average_motif_length, dispersion)

        # sys.stderr.write(f"Cluster {n_cluster:2} motif template: {motif_notation(motif_template)}\n")
        total_length = mean_range(average_random_length * 2, dispersion) + motif_length
        prefix_length = mean_range(average_random_length, dispersion // 2)
        suffix_length = total_length - motif_length - prefix_length
        sys.stderr.write(f"Generating sequences for cluster {n_cluster}\n")
        for n_seq, seq, activity, is_cliff in generate_cluster(
            n_sequences,
            motif_length,
            prefix_length,
            suffix_length,
            max_variants_per_position,
            make_cliffs,
            alphabet,
            cliff_probability,
            cliff_strength,
        ):
            sequences.append(
                (n_cluster, f"c{n_cluster}_s{n_seq:03d}", seq, activity, is_cliff)
            )
    return headers, sequences


def convert_to_fasta(
    cluster_sequence_records: List[ClusterSequenceRecord], separator: str
) -> List[Tuple[int, str, str, float, bool]]:
    return [
        (n_cluster, name_cluster, sequence_to_fasta(seq, separator), activity, is_cliff)
        for n_cluster, name_cluster, seq, activity, is_cliff in cluster_sequence_records
    ]


def convert_to_helm(
    cluster_sequence_records: List[ClusterSequenceRecord], helm_connection_mode: str
) -> List[Tuple[int, str, str, float, bool]]:
    return [
        (
            n_cluster,
            name_cluster,
            sequence_to_helm(seq, helm_connection_mode),
            activity,
            is_cliff,
        )
        for n_cluster, name_cluster, seq, activity, is_cliff in cluster_sequence_records
    ]


def is_monomer_suitable(monomer: Any) -> bool:
    return (
        monomer["polymerType"] == "PEPTIDE"
        and monomer["monomerType"] == "Backbone"
        and len(monomer["rgroups"]) == 2
    )


def alphabet_from_helm(helm_library_file: str) -> Alphabet:
    import json

    alphabet: Alphabet = []
    with open(helm_library_file) as helm_library:
        for monomer in json.load(helm_library):
            if is_monomer_suitable(monomer):
                alphabet.append(monomer["symbol"])
    return alphabet


def parse_command_line_args() -> Any:
    parser = argparse.ArgumentParser(
        prog="MotifSequencesGenerator",
        description="The program generates set of sequences containing sequence motifs "
        "for SAR functionality testing",
        epilog="Utility author and support: Gennadii Zakharov <Gennadiy.Zakharov@gmail.com>",
    )

    parser.add_argument(
        "-c", "--clusters", type=int, default=5, help="Number of clusters"
    )
    parser.add_argument(
        "-s",
        "--sequences",
        type=int,
        default=50,
        help="Number of sequences in each supercluster",
    )
    parser.add_argument(
        "-m,", "--motif-length", type=int, default=12, help="Average length of motif"
    )

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

    parser.add_argument(
        "-h,",
        "--helm-library-file",
        type=str,
        help="JSON file containing the HELM monomer library in the same format as used for Datagrok. "
        + "The alphabet property is ignored when helm library is specified.",
    )

    parser.add_argument(
        "--helm-connection-mode",
        type=str,
        default=HelmConnectionMode.linear.value,
        help=f"HELM peptide generation mode: {'/'.join([mode.name for mode in HelmConnectionMode])}",
    )

    available_alphabets = ",".join(list(alphabets.keys()) + ["custom"])
    parser.add_argument(
        "--alphabet",
        type=str,
        default=list(alphabets.keys())[0],
        help=f"Sequence alphabet: {available_alphabets}. Custom alphabet is a list of values separated "
        f"by comma",
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
    parser.add_argument(
        "--fasta-separator",
        type=str,
        default="",
        help="Separator symbol for FASTA sequence",
    )
    command_line_args = parser.parse_args()

    return command_line_args


# ====================================================================================

grok = "clusters" in globals()

if not grok:
    # We are not in Datagrok - need to parse command line arguments
    args = parse_command_line_args()
    clusters = args.clusters
    num_sequences = args.sequences
    motif_length = args.motif_length
    max_variants_position = args.max_variants_position
    random_length = args.random_length
    dispersion = args.dispersion
    alphabet_key = args.alphabet
    disable_cliffs = args.disable_cliffs
    cliff_probability = args.cliff_probability
    cliff_strength = args.cliff_strength
    fasta_separator = args.fasta_separator
    helm_library_file = args.helm_library_file
    helm_connection_mode = args.helm_connection_mode

helm_init = "helm_library_file" in globals() and helm_library_file is not None

if not helm_init:
    alphabet: Alphabet = (
        alphabets[alphabet_key].split(",")
        if alphabet_key in alphabets
        else alphabet_key.split(",")
    )
else:
    alphabet = alphabet_from_helm(helm_library_file)

# Running sequence generator
header, data = generate_sequences(
    clusters,
    num_sequences,
    motif_length,
    max_variants_position,
    random_length,
    dispersion,
    alphabet,
    not disable_cliffs,
    cliff_probability,
    cliff_strength,
)
if not helm_init:
    data_formatted = convert_to_fasta(data, fasta_separator)
else:
    data_formatted = convert_to_helm(data, helm_connection_mode)

if grok:
    # Exporting data to Datagrok as a Pandas dataframe
    import pandas as pd

    sequences = pd.DataFrame.from_records(data_formatted, columns=header)
else:
    # Writing results to stdout - no need to work with big and heavy Pandas
    import csv

    csv_writer = csv.writer(sys.stdout, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
    csv_writer.writerow(header)
    for line in data_formatted:
        csv_writer.writerow(line)
