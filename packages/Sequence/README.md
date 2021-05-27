Sequence is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.
It showcases various bioinformatics functionality related to sequences with using 3-rd party libraries.


##FASTA Format support


Sequence Package supports import of FASTA files. The file is loaded into a dataframe (table) where all sequence lines belonging to one molecule are concatenated into a text field.

FASTA description line is parsed to extract sequence id, description and metainformation (features) of the defline like “[taxid=9606]”. No feature name or value checks are performed.

FASTA reader supports both nucleotides and proteins or alphabets. 
FASTA reader supports ambiguity codes ‘N’, ‘M’ and alignment gap symbols.
