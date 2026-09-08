/* The datasets Bio's features open: the package's own test and sample files on the stand
   (`System:AppData/Bio/...`, published with the package). Platform names (toolbox, grid, ...)
   are reserved; a Bio view's own names go on a context here when a feature needs one. */
import {dataset} from '@datagrok-libraries/bdd';

dataset('filter_FASTA', {path: 'System:AppData/Bio/tests/filter_FASTA.csv', aliases: ['filter-fasta'],
  description: 'one column "fasta" (peptides, fasta notation, 9 sequences + 4 empty cells, 14 rows)'});
dataset('filter_HELM', {path: 'System:AppData/Bio/tests/filter_HELM.csv', aliases: ['filter-helm'],
  description: 'one column "HELM string" (helm notation, 3 rows, one bracketed monomer [dV])'});
dataset('filter_MSA', {path: 'System:AppData/Bio/tests/filter_MSA.csv', aliases: ['filter-msa'],
  description: 'columns "MSA" (aligned, separator "/", multichar monomers, width 17) and "Activity"'});
dataset('antibodies', {path: 'System:AppData/Bio/samples/antibodies.csv',
  description: '493 antibodies: AntibodyHC and AntibodyLC (fasta peptides ~130–215 aa), Antigen, Y'});
dataset('FASTA_PT_activity', {path: 'System:AppData/Bio/samples/FASTA_PT_activity.csv', aliases: ['fasta-pt-activity', 'peptides with activity'],
  description: '99 peptides: cluster, sequence_id, sequence (16-mers), activity, is_cliff'});
dataset('helm_cyclic_cliffs', {path: 'System:AppData/Bio/tests/helm_cyclic_cliffs.csv', aliases: ['helm-cyclic-cliffs']});
dataset('FASTA_sample', {path: 'System:AppData/Bio/samples/FASTA.csv', aliases: ['fasta-sample'],
  description: '64 UniProt peptides: Entry, Length, UniProtKB, Sequence (fasta), Activity, Cluster'});
