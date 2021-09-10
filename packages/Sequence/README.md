Support for the biological sequences (DNA, RNA, proteins), including FASTA readers, 
3D visualization of molecules, MSA viewer, and functions for working with sequences.

## Fasta format support

Sequence Package supports import of FASTA files. The file is loaded into a dataframe (table) where all sequence lines belonging to one molecule are concatenated into a text field.

- description line is parsed to extract sequence id, description and metainformation (features) of the defline like “[taxid=9606]”. 
  No feature name or value checks are performed.
- supports both nucleotides and proteins or alphabets
- supports ambiguity codes ‘N’, ‘M’ and alignment gap symbols
