#!/usr/bin/env python3
# name: Sequence generator
# description: Create the model peptides/DNA sequences with peptides data
# language: python
# tags: template, demo
# input: int clusters = 5 { caption: Clusters; category: Clusters } [Number of clusters]
# input: int num_sequences = 50 { caption: Sequences; category: Clusters } [Number of sequences in each cluster]
# input: string alphabet_key = "Protein" {  choices: ["Protein", "DNA", "RNA", "Protein_EXT"]; caption: Alphabet; category: Clusters;} [Sequence alphabet. Ignored if the HELM library is specified.]
# input: int motif_length = 12 { caption: Motif length; category: Motif } [Average length of motif]
# input: int max_variants_position = 3 { caption: Position variants; category: Motif } [Maximum number of different letters in a conservative position of the motif]
# input: int random_length = 3 { caption: Randon length; category: Motif } [Average length of random sequence parts before and after motif]
# input: int dispersion = 2 { caption: Length variation; category: Motif } [Variation of total sequence length]
# input: double activity_range = 0.2 { caption: Activity range; category: Activity parameters; format: 0.000} [Range of the mean activity value difference between clusters]
# input: double cliff_probability = 0.05 { caption: Cliff probability; category: Activity parameters; format: 0.000} [Probability to make activity cliff of a sequence]
# input: double cliff_strength = 5.0 { caption: Cliff strength; category: Activity parameters } [The size of the cliff comparing to the dispersion of the initial activity]
# input: double cliff_strength_dispersion = 1.0 { caption: Cliff dispersion; category: Activity parameters } [Dispersion of cliff strength]
# input: string assay_noise_levels = "0.4, 0.85" { caption: Noise levels; category: Assay settings } [List of assay noise levels, separated by comma]
# input: string assay_scales = "(0|10), (0|150.0)"  { caption: Assay scales; category: Assay settings } [Typical scale size for each assay. Assays are separated by comma. Minimum and maximum values are separated by pipe. Brackets are optional]
# input: bool disable_negatives = true { caption: Crop negatives; category: Assay settings } [Set negative measurements for assay to zero]
# input: string fasta_separator = "" { caption: Fasta separator; nullable: true; category: Output format} [Monomers separator for FASTA format]
# input: file helm_library_file { caption: HELM library; nullable: true; category: Output format} [HELM library to load alphabet. Output format is set to HELM if the HELM library is specified]
# input: string helm_connection_mode = "linear" { choices: ["linear", "cyclic", "mixed"]; caption: Connection mode; category: Output format} [Peptides connection mode for HELM output)]
# output: dataframe sequences_data

description = """The utility generates clusters of macromolecule sequences to test SAR functionality. 
Each cluster contains a randomly generated sequence motif.
Each sequence has activity - a Gauss-distributed random value. 
The utility can simulate activity cliffs - random changes in the conservative motif letters,
leading to the significant change in the activity.
Utility can simulate multiple experimental assays measuring activity, with different scales and noise levels."""

import random
import argparse
import sys
from collections import namedtuple
from enum import Enum

from typing import List, Tuple, NamedTuple, Dict, Set, Any

# ===== Type definitions =====
Letter = str
Alphabet = List[Letter]
LetterChoice = List[Letter]
MotifTemplate = List[LetterChoice]

# The sequence in a list of a monomers from the alphabet.
# We can't use string because monomers can have several letters
Sequence = List[Letter]
SequenceList = List[Sequence]
SequenceSquashed = str  # Sequence, joined together in string form

CliffPair = Tuple[int, int]
CliffList = List[CliffPair]

Activity = float
ActivityList = List[Activity]

ClusterParameters = NamedTuple(
    "ClusterParameters",
    [
        ("motif_length", int),
        ("max_variants_per_position", int),
        ("random_length", int),
        ("dispersion", int),
    ],
)
CliffParameters = namedtuple(
    "CliffParameters",
    ["cliff_probability", "cliff_strength", "cliff_strength_dispersion"],
)
AssayParameters = NamedTuple(
    "AssayParameters", [("noise_level", float), ("min", float), ("max", float)]
)

DataLine = Tuple[
    Any, ...
]  # Contains strings and 1+ number of floats - can't type more exactly

# ===== Constants =====
OutputFormat = Enum("OutputFormat", ["Fasta", "Helm"])
HelmConnectionMode = Enum("HelmConnectionMode", ["linear", "cyclic", "mixed"])

alphabets: Dict[str, Alphabet] = {
    "Protein": "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y".split(","),
    "DNA": "A,T,G,C".split(","),
    "RNA": "A,U,G,C".split(","),
    "Protein_EXT": "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,dA,dC,dD,dE,dF,dH,dI,dK,dL,dM,dN,dP,dQ,dR,dS,dT,dV,dW,dY,meA,meD,meS,meT,meV,meY,meE,meG,meI,meK,meM,meN,meQ,meC,meR,meW,meF,meH,meL,Nle,Nva,Orn,Iva,aIle,gGlu,Hcy,Hse,Hyp,D-gGlu,D-Nle,D-hPhe,D-Hyp,D-Nva,D-Orn,Pyr,Phe_3Cl,Phe_4Cl,Phe_4NH2,Phg,Ser_tBu,Tyr_Bn,Tza,1Nal,Cha,Lys_Boc,aThr,D-2Nal,D-2Thi,D-aHyp,D-aIle,D-Phg,D-Ser_tBu,Cya,Lys_Me3,Pen,Phe_4Me,Ser_Bn,Tyr_tBu,2Nal,Thi,aHyp,Ala_tBu,hPhe,D-1Nal,D-aThr,D-Cha,D-Pen,D-Phe_4Cl,D-Ser_Bn,Wil,Oic_3aS-7aS,Pip,3Pal,4Pal,Abu,Apm,Chg,Dab,Dap,D-3Pal,D-aMeAbu,D-Chg,D-Cit,D-Dab,D-Pip,D-Tic,Aca,Tic,Aad,Cit,Aze,Ac5c,Aib,D-2Pal,D-Abu,D-Dap,Asu,D-Thz,D-Trp_For,D-Tyr_Et,Lys_Ac,Asp_OMe,Phe_ab-dehydro,Sta_3xi4xi,Tyr_ab-dehydroMe,App,Cap,Cys_SEt,Dsu,pnC,pnG,Pqa,Pro_4Me3OH,Met_O2,Phe_2Me,Phe_34diCl,Phe_4Br,Phe_4I,Phe_4Sdihydroorotamido,Pyl,Ser_PO3H2,Thr_PO3H2,Thz,Trp_Me,Tyr_26diMe,Tyr_3I,Tyr_3NO2,Tyr_Ph4OH,Tyr_SO3H,Val_3OH,xiIle,NMe2Abz,NMebAla,aMePhe,aMePro,aMeTyr_3OH,Bmt,Bmt_E,Cys_Bn,Gla,hHis,His_1Me,Gly_allyl,Gly_cPr,Asp_Ph2NH2,Azi,2Abz,3Abz,4Abz,Ac3c,Ac6c,bAla,D-Bmt,D-Bmt_E,D-hArg,D-Phe_4F,D-Trp_2Me,D-Tyr_Me,D-xiIle,Lys_iPr,Phe_ab-dehydro_3NO2,Sta_3S4S,Bux,Dpm,pnA,pnT,seC,Met_O,nTyr,Oic_3aR-7aS,Oic_3axi-7axi,Phe_2F,Phe_3F,Phe_4F,Phe_4NO2,Phe_bbdiMe,Trp_5OH,Trp_Ome,Tyr_35diI,Tyr_3OH,Tyr_Me,Tyr_PO3H2,xiHyp,xiThr,NMe4Abz,aMeTyr,Aoda,Bpa,Cys_Me,Dip,hArg,His_1Bn,His_3Me,Hyl_5xi,Bip,Abu_23dehydro,D-Dip,Dha,D-hArg_Et2,D-Met_S-O,D-His_1Bn,D-nTyr,D-Phe_4ureido".split(
        ","
    ),
}

# ===== Motif and sequence generation functions =====


def alphabet_from_helm(helm_library_file: str) -> Alphabet:
    """
    Reads the HELM library from a JSON file and extracts only backbone monomers suitable for sequence generation
    """
    import json

    def is_monomer_suitable(monomer: Any) -> bool:
        return (
            monomer["polymerType"] == "PEPTIDE"
            and monomer["monomerType"] == "Backbone"
            and len(monomer["rgroups"]) == 2
        )

    alphabet: Alphabet = []
    with open(helm_library_file) as helm_library:
        for monomer in json.load(helm_library):
            if is_monomer_suitable(monomer):
                alphabet.append(monomer["symbol"])
    return alphabet


def mean_range(mean: int, disp: int) -> int:
    """
    Returns random positive value around some mean with selected dispersion
    """
    return random.randint(max(mean - disp, 0), mean + disp)


def generate_motif_template(
    motif_length: int,
    alphabet: Alphabet,
    max_variants_cluster: int,
    prob_any: float = 0.2,  # The probability to have a non-conservative letter (the `?` sign in notation) inside motif
) -> MotifTemplate:
    """
    Generated random template from the alphabet
    """
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
    """
    Generate sequence motif by motif template
    """
    template_with_any = [
        (letters if not "?" in letters else alphabet) for letters in template
    ]
    return [random.choice(letters) for letters in template_with_any]


def motif_notation(motif_template: MotifTemplate, fasta_separator: str = "") -> str:
    """
    Returns string representation of motif template
    """

    def motif_notation_code(letter_choice: LetterChoice) -> str:
        if len(letter_choice) == 1:
            return letter_choice[0] + fasta_separator
        else:
            return f"[{fasta_separator.join(letter_choice)}]"

    return "".join(
        [motif_notation_code(letter_choice) for letter_choice in motif_template]
    )


def generate_random_sequence(n: int, alphabet: Alphabet) -> Sequence:
    """
    Generate a sequence containing n random letters from the alphabet
    """
    return [random.choice(alphabet) for i in range(n)]


def make_motif_cliff(
    motif_template: MotifTemplate, alphabet: Alphabet, motif: Sequence
) -> Sequence:
    """
    Mutates a random conservative letter in the motif
    """
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


def generate_cluster_sequences(
    n_sequences: int,
    motif_template: MotifTemplate,
    prefix_length: int,
    suffix_length: int,
    alphabet: Alphabet,
    cliff_probability: float,
) -> Tuple[SequenceList, CliffList]:
    """
    Returns set of sequences for one cluster and introduces sequence cliffs
    Also makes activity cliffs
    """
    n_seq = 0
    sequences: SequenceList = []
    cliffs: CliffList = []

    while n_seq < n_sequences:
        motif = generate_motif(motif_template, alphabet)
        prefix = generate_random_sequence(prefix_length, alphabet)
        suffix = generate_random_sequence(suffix_length, alphabet)
        seq = prefix + motif + suffix
        sequences.append(seq)
        n_seq += 1
        if n_seq >= n_sequences:
            break  # This is the last sequence - can't do cliff
        is_cliff = random.random() <= cliff_probability
        if is_cliff:
            # Making activity cliff
            cliff_motif = make_motif_cliff(motif_template, alphabet, motif)
            cliff_seq = prefix + cliff_motif + suffix
            sequences.append(cliff_seq)
            cliffs.append((n_seq - 1, n_seq))
            n_seq += 1
            # sys.stderr.write(f"Cliff for sequence #{line_number:4}, cluster {n_cluster} \n")
    return sequences, cliffs


def sequence_to_fasta(sequence: Sequence, separator: str) -> SequenceSquashed:
    """
    Converts the sequence to FASTA format
    """
    return separator.join(sequence)


def sequence_to_helm(
    sequence: Sequence,
    helm_connection_mode: HelmConnectionMode = HelmConnectionMode.linear,
) -> SequenceSquashed:
    """
    Converts the sequence to HELM format
    """

    def is_cyclic(helm_connection_mode: HelmConnectionMode) -> bool:
        return helm_connection_mode == HelmConnectionMode.cyclic or (
            helm_connection_mode == HelmConnectionMode.mixed and random.random() < 0.5
        )

    sequence_escaped: Sequence = [
        f"[{letter}]" if len(letter) > 1 else letter for letter in sequence
    ]
    connection_format = ""
    if is_cyclic(helm_connection_mode):
        connection_format = f"PEPTIDE1,PEPTIDE1,{len(sequence_escaped)}:R2-1:R1"
    return f"PEPTIDE1{{{sequence_to_fasta(sequence_escaped,'.')}}}${connection_format}$$$V2.0"


# ===== Activity generation functions =====
def generate_ideal_activities(n: int, activity_range: float = 0) -> ActivityList:
    """
    Generate ideal activities with Gauss distribution
    The distribution center is chosen randomly with some dispersion
    """
    mean = random.uniform(-activity_range, activity_range) if activity_range > 0 else 0
    return [random.gauss(mean, 1) for _ in range(n)]


def make_activity_cliff(
    activities: ActivityList,
    cliffs: List[CliffPair],
    cliff_strength: float,
    cliff_strength_dispersion: float,
) -> ActivityList:
    """
    Introduce activity cliffs -
    make a pair of activities differ for random gauss-distributed value defined by cliff_strength and cliff_strength_dispersion
    """
    cliff_activities = activities[:]
    for first, second in cliffs:
        activity1 = activities[first]
        activity2 = activities[second]
        average = (activity1 + activity2) / 2
        scale = random.gauss(cliff_strength, cliff_strength_dispersion) / abs(
            activity1 - activity2
        )
        cliff_activities[first] = average + (activity1 - average) * scale
        cliff_activities[second] = average + (activity2 - average) * scale
    return cliff_activities


def generate_assay_activities(
    activities: ActivityList,
    assay: AssayParameters,
    disable_negatives: bool = True,
) -> ActivityList:
    """
    Generates activities measured in assay from some "ideal" activities.
    Adds noise and scales the values to emulate some assay measurement scale
    """
    assay_activities = []
    scale_factor = 3 * (
        1 + assay.noise_level
    )  # real activity 3-sigma in the interval [-scale_factor,+scale_factor]
    for activity in activities:
        noise = random.uniform(
            -3, 3
        )  # some random noize in [-3,3] - 3 sigma for ideal activity
        # Adding noize and normalizing
        noised_activity = activity + noise * assay.noise_level
        rescaled_activity = (
            noised_activity / (scale_factor * 2)
        ) + 0.5  # rescaling activity to the interval [0;1]

        assay_result = assay.min + (rescaled_activity * (assay.max - assay.min))

        if disable_negatives and assay_result < 0:
            assay_result = 0

        assay_activities.append(assay_result)
    return assay_activities


def generate_data(
    n_clusters: int,
    n_sequences: int,
    cluster_parameters: ClusterParameters,
    assays: List[AssayParameters],
    disable_negatives: bool,
    alphabet: Alphabet,
    output_format: OutputFormat,
    fasta_separator: str,
    helm_connection_mode: HelmConnectionMode,
    activity_range: float,
    cliff_probability: float = 0.05,
    cliff_strength: float = 5.0,
    cliff_dispersion: float = 1.0,
) -> Tuple[List[str], List[DataLine]]:
    """
    Main function generating all data set - sequences, activities, etc
    """
    headers: List[str] = ["cluster", "sequence_id", "sequence", "is_cliff"]
    headers += [f"Assay_{i+1}" for i in range(len(assays))]
    data: List[DataLine] = []

    def cliffs_to_positions(cliffs: CliffList) -> Set[int]:
        """
        Convert CliffList to a set containing positions of cliffs
        """
        unique_pos = {pos for cliff in cliffs for pos in cliff}
        return unique_pos

    for n_cluster in range(n_clusters):
        motif_length = mean_range(
            cluster_parameters.motif_length, cluster_parameters.dispersion
        )
        total_length = (
            mean_range(
                cluster_parameters.random_length * 2, cluster_parameters.dispersion
            )
            + motif_length
        )
        prefix_length = mean_range(
            cluster_parameters.random_length, cluster_parameters.dispersion // 2
        )
        suffix_length = total_length - motif_length - prefix_length

        # Making a motif template
        motif_template = generate_motif_template(
            motif_length, alphabet, cluster_parameters.max_variants_per_position
        )
        sys.stderr.write(
            f"Motif template for cluster {n_cluster}: {motif_notation(motif_template, fasta_separator)}\n"
        )

        sequences, cliffs = generate_cluster_sequences(
            n_sequences,
            motif_template,
            prefix_length,
            suffix_length,
            alphabet,
            cliff_probability,
        )

        if output_format == OutputFormat.Fasta:
            squashed_sequences = [
                sequence_to_fasta(seq, fasta_separator) for seq in sequences
            ]
        elif output_format == OutputFormat.Helm:
            squashed_sequences = [
                sequence_to_helm(seq, helm_connection_mode) for seq in sequences
            ]
        else:
            print("Unsupported output format")
            exit(-1)

        ideal_activities = generate_ideal_activities(n_sequences, activity_range)
        cliffed_activities = make_activity_cliff(
            ideal_activities, cliffs, cliff_strength, cliff_dispersion
        )

        assay_activities = [
            generate_assay_activities(cliffed_activities, assay, disable_negatives)
            for assay in assays
        ]

        cliffs_positions = cliffs_to_positions(cliffs)
        is_cliffs = [pos in cliffs_positions for pos in range(n_sequences)]
        sequence_IDs = [f"c{n_cluster}_s{n:03d}" for n in range(n_sequences)]

        cluster_data = zip(
            [n_cluster] * n_sequences,
            sequence_IDs,
            squashed_sequences,
            is_cliffs,
            *assay_activities,
        )

        data.extend(cluster_data)

    return headers, data


def repack_assays(noise_levels_str: str, scales_str: str) -> List[AssayParameters]:
    """
    Converts strings passed from the input data to the list of AssayParameters namedtuples
    """
    noise_levels = [float(s) for s in noise_levels_str.split(",")]
    scales = [s.strip().split("|") for s in scales_str.split(",")]
    minmaxes = [(float(x[0].strip("() ")), float(x[1].strip("()"))) for x in scales]
    if not (len(noise_levels) == len(minmaxes)):
        print("Not equal range of parameters for assay definition")
        exit(-1)
    assays = [
        AssayParameters(noise, min, max)
        for noise, (min, max) in zip(noise_levels, minmaxes)
    ]
    return assays


# ===== Tests =====


def test_activities_correlation() -> None:
    import numpy as np

    ideal_activities = generate_ideal_activities(25, 0.1)
    cliff_activities = make_activity_cliff(
        ideal_activities, [(0, 1)], cliff_strength=5.0, cliff_strength_dispersion=1.0
    )
    assay_parameters = AssayParameters(0.3, 0, 10)
    x = generate_assay_activities(cliff_activities, assay_parameters)
    assay_parameters = AssayParameters(0.5, 0, 250)
    y = generate_assay_activities(cliff_activities, assay_parameters)

    print("Assay1: " + ",".join([str(a) for a in x]))
    print("Assay2: " + ",".join([str(a) for a in y]))
    corr = np.corrcoef(x, y)
    print("Correlation: ", corr[1, 0])
    assert corr[1, 0] >= 0.5


# ===== Command-line arguments parsing =====


def parse_command_line_args() -> Any:
    parser = argparse.ArgumentParser(
        prog="MotifSequencesGenerator",
        description=description,
        epilog="Utility author and support: Gennadii Zakharov <Gennadiy.Zakharov@gmail.com>",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    cluster_group = parser.add_argument_group("Cluster parameters")

    cluster_group.add_argument(
        "-c", "--clusters", type=int, default=5, help="Number of clusters"
    )
    cluster_group.add_argument(
        "-s",
        "--sequences",
        type=int,
        default=50,
        help="Number of sequences in each supercluster",
    )

    available_alphabets = ",".join(list(alphabets.keys()))
    cluster_group.add_argument(
        "--alphabet",
        type=str,
        default=list(alphabets.keys())[0],
        help=f"Sequence alphabet: {available_alphabets}.\n"
        + "Ignored if the HELM library is specified",
    )

    motif_group = parser.add_argument_group("Motif parameters")

    motif_group.add_argument(
        "-m,", "--motif-length", type=int, default=12, help="Average length of motif"
    )

    motif_group.add_argument(
        "--max-variants-position",
        type=int,
        default=3,
        help="Maximum number of different letters in a conservative position of the motif",
    )

    motif_group.add_argument(
        "-r,",
        "--random-length",
        type=int,
        default=3,
        help="Average length of random sequence parts before and after motif",
    )
    motif_group.add_argument(
        "-d,",
        "--dispersion",
        type=int,
        default=2,
        help="Variation of total sequence length",
    )

    cliffs_group = parser.add_argument_group("Activity parameters")

    cliffs_group.add_argument(
        "--activity-range",
        type=float,
        default=0.5,
        help="Range of the mean activity value difference between clusters",
    )

    cliffs_group.add_argument(
        "--cliff-probability",
        type=float,
        default=0.05,
        help="Probability to make activity cliff of a sequence",
    )
    cliffs_group.add_argument(
        "--cliff-strength",
        type=float,
        default=5.0,
        help="Average strength of cliff",
    )

    cliffs_group.add_argument(
        "--cliff-strength-dispersion",
        type=float,
        default=1.0,
        help="Cliff strength dispersion",
    )

    assay_group = parser.add_argument_group("Assay parameters")

    assay_group.add_argument(
        "--assay-noise-levels",
        type=str,
        default="0.4, 0.85",
        help="Noise level(s) for assays. A list of values separated by comma.",
    )

    assay_group.add_argument(
        "--assay-scales",
        type=str,
        default="(0|10), (0|150.0)",
        help="Typical scale size for each assay. Assays are separated by comma. Minimum and maximum values are separated by pipe. Brackets are optional."
        + "Activity outliers may be located outside the specified scale",
    )

    assay_group.add_argument(
        "--enable-negatives",
        type=bool,
        help="Enable negative values for assays results",
    )

    output_group = parser.add_argument_group("Output parameters")

    output_group.add_argument(
        "--custom-alphabet",
        type=str,
        default="",
        help=f"Custom sequence alphabet: list of letters separated by comma. Used only if the --alphabet=custom",
    )

    output_group.add_argument(
        "--fasta-separator",
        type=str,
        default="",
        help="Monomers separator for FASTA format",
    )

    output_group.add_argument(
        "-H,",
        "--helm-library-file",
        type=str,
        help="JSON file containing the HELM monomer library. "
        + "The alphabet property is ignored when helm library is specified.",
    )

    output_group.add_argument(
        "--helm-connection-mode",
        type=str,
        default=HelmConnectionMode.linear.name,
        help=f"HELM peptide generation mode: {'/'.join([mode.name for mode in HelmConnectionMode])}",
    )

    command_line_args = parser.parse_args()

    return command_line_args


# ===== Main part of script =====

if __name__ == "__main__":
    grok = "clusters" in globals()

    if not grok:
        # We are not in Datagrok - need to parse command line arguments
        args = parse_command_line_args()
        #
        clusters = args.clusters
        num_sequences = args.sequences
        alphabet_key = args.alphabet
        #
        motif_length = args.motif_length
        max_variants_position = args.max_variants_position
        random_length = args.random_length
        dispersion = args.dispersion
        #
        activity_range = args.activity_range
        cliff_probability = args.cliff_probability
        cliff_strength = args.cliff_strength
        cliff_strength_dispersion = args.cliff_strength_dispersion
        #
        assay_noise_levels = args.assay_noise_levels
        assay_scales = args.assay_scales
        disable_negatives = not args.enable_negatives
        #
        custom_alphabet = args.custom_alphabet
        fasta_separator = args.fasta_separator
        helm_library_file = args.helm_library_file
        helm_connection_mode = args.helm_connection_mode

    helm_init = helm_library_file is not None and helm_library_file != ""

    if helm_init:
        alphabet = alphabet_from_helm(helm_library_file)
        output_format = OutputFormat.Helm
        fasta_separator = "|"
    else:
        output_format = OutputFormat.Fasta
        if not alphabet_key in alphabets:
            pass  # TBD: custom alphabet
        alphabet = alphabets[alphabet_key]

    # Packing parameters to structures to simplify function signatures
    cluster_parameters = ClusterParameters(
        motif_length, max_variants_position, random_length, dispersion
    )
    assays = repack_assays(assay_noise_levels, assay_scales)

    # Running sequence generator
    header, data = generate_data(
        clusters,
        num_sequences,
        cluster_parameters,
        assays,
        disable_negatives,
        alphabet,
        output_format,
        fasta_separator,
        HelmConnectionMode[helm_connection_mode],
        activity_range,
        cliff_probability,
        cliff_strength,
        cliff_strength_dispersion,
    )

    if grok:
        # Exporting data to Datagrok as a Pandas dataframe
        import pandas as pd

        sequences_data = pd.DataFrame.from_records(data, columns=header)
    else:
        # Writing results to stdout - no need to work with big and heavy Pandas
        import csv

        csv_writer = csv.writer(sys.stdout, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(header)
        for line in data:
            csv_writer.writerow(line)
