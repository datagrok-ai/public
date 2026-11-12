/* The datasets Bio's features open: the package's own test and sample files on the stand
   (`System:AppData/Bio/...`, published with the package). Platform names (toolbox, grid, ...)
   are reserved; a Bio view's own names go on a context here when a feature needs one. */
import {dataset, element} from '@datagrok-libraries/bdd';

dataset('filter_FASTA', {path: 'System:AppData/Bio/tests/filter_FASTA.csv', aliases: ['filter-fasta'],
  description: 'one column "fasta" (peptides, fasta notation, 9 sequences + 4 empty cells, 14 rows)'});
dataset('filter_HELM', {path: 'System:AppData/Bio/tests/filter_HELM.csv', aliases: ['filter-helm'],
  description: 'one column "HELM string" (helm notation, 4 rows, one bracketed monomer [dV])'});
dataset('filter_MSA', {path: 'System:AppData/Bio/tests/filter_MSA.csv', aliases: ['filter-msa'],
  description: 'columns "MSA" (aligned, separator "/", multichar monomers, width 17) and "Activity"'});
dataset('antibodies', {path: 'System:AppData/Bio/samples/antibodies.csv',
  description: '493 antibodies: AntibodyHC and AntibodyLC (fasta peptides ~130–215 aa), Antigen, Y'});
dataset('FASTA_sample', {path: 'System:AppData/Bio/samples/FASTA.csv', aliases: ['fasta-sample'],
  description: '64 UniProt peptides: Entry, Length, UniProtKB, Sequence (fasta), Activity, Cluster'});
dataset('HELM_sample', {path: 'System:AppData/Bio/samples/HELM.csv', aliases: ['helm-sample'],
  description: '540 peptides: HELM (helm notation, multichar monomers), Activity'});
dataset('MSA_sample', {path: 'System:AppData/Bio/samples/MSA.csv', aliases: ['msa-sample'],
  description: '540 aligned peptides: MSA (separator "/", multichar monomers), Activity'});
dataset('cyclized', {path: 'System:AppData/SequenceTranslator/samples/cyclized.csv',
  description: 'SequenceTranslator\'s sample: "seqs", separator sequences with cyclized monomers ("C(1)"), which its notation refiner turns into custom units'});
dataset('BILN', {path: 'System:AppData/Bio/samples/BILN.csv', aliases: ['biln-sample'],
  description: 'column "biln" (BILN notation: "-" between monomers, "(n,m)" ring bonds, "." between chains) and Activity'});
dataset('BILN_complex', {path: 'System:AppData/Bio/samples/BILN_complex.csv', aliases: ['biln-complex'],
  description: '120 peptides in BILN with their connections: name, category (disulfide, multi-chain, head-to-tail, side-chain ring, branched, scaffold, linear, long names), biln, notes'});
dataset('HELM_cyclic', {path: 'System:AppData/Bio/samples/HELM_cyclic.tsv', aliases: ['helm-cyclic'],
  description: '41 head-to-tail cyclic HELM peptides of non-natural monomers: cluster, sequence_id, sequence, activity, is_cliff'});
dataset('FASTA_DNA', {path: 'System:AppData/Bio/samples/FASTA_DNA.csv', aliases: ['fasta-dna'],
  description: '100 DNA sequences: id, sequence (fasta)'});
dataset('FASTA_RNA', {path: 'System:AppData/Bio/samples/FASTA_RNA.csv', aliases: ['fasta-rna'],
  description: '100 RNA sequences: id, sequence (fasta)'});
dataset('SEPARATOR_PT', {path: 'System:AppData/Bio/samples/SEPARATOR_PT.csv', aliases: ['separator-pt'],
  description: '99 peptides: id, sequence (separator "-")'});

element('monomer libraries link', {selector: '.d4-link-action[data-source="Bio:Bioinformatics | Manage Monomer Libraries"]',
  description: 'the "Manage monomer libraries" link of a Macromolecule column\'s Bioinformatics > Manage Monomer Libraries pane: the dialog, or the Manage Monomer Libraries view when one is open'});
