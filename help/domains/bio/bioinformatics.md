---
title: "Bioinformatics"
---

Bioinformatics (see also the [Wikipedia](https://en.wikipedia.org/wiki/Bioinformatics) article) is the application of
computer and information science methods for understanding biological data, especially when the data sets are large 
and complex. 

Bioinformatics itself is a
complex discipline, serving multiple areas, such as:

* New generation clinical diagnostics 
* Personalized predictive medicine
* Cancer mutations analysis, retargeting of existing anti-cancer medicines
* Analysis of biological (especially antimicrobial and anti-cancer) activity of natural biomolecules 
* Design of new biomolecules with target activity
* High-throughput analysis of biological images
* Mining data and making insights from scientific publications  
* etc.

The associated _in silico_ techniques be very different in terms of input data, analysis procedures,
and final results.

Datagrok supports large set of bioinformatics capabilities, 
especially for analysis of sequence-activity relations (SAR) for macromolecules.

## Biological data format support

DataGrok automatically determines data format and provides custom loaders for several biological data formats:
FASTA (protein and DNA/RNA), HELM, PDB, Newick (NWK).
 
This functionality is not bound to specific file extension. 
For example DataGrok can find columns containing DNA/peptide sequences in any dataframe.

## Custom data rendering

For each supported biological data type DataGrok provides a custom renderer working directly inside a dataframe: 

* Peptide sequences are colored basing on aminoacid properties. 
* DNA sequences are colored by nucleotide. The user can control coloring by changing the alphabet.
* HELM structures are colored by the monomer. HELM renderer supports circular and branching structures. 
If monomer library for HELM structure is loaded 
(to do this open column properties by left-clicking on the column header 
and select in the column menu “Manage libraries”) 
the macromolecule preview and editing become availiable. HELM renderer and loader supports circular and branching structures.
* PDB structures has graphical preview. By clicking on the cell you can open separate 3D viewer, 
supporting rotating, zooming and changing the color scheme.    

## Data editing

For several data types the data can be modified directly in the dataframe. 

* You can edit FASTA sequences for DNA, RNA and proteins
* If the HELM monomer library is loaded, you can edit HELM sequences: add/remove monomers, change connections. 
The editor supports circular and branching structures.

## Data filtering 

For sequence column you can specify the subsequence and filter only sequences containing this subsequence.

## Specific viewers and actions for macromolecules

Datagrok has set of specifc actions for macromolecules. 
All of them are in the "Bio" section of the main menu.

### MSA (alignment)

Make a multiple sequence alignment. Create a new column with aligned sequences.
Important note: the alignment utility supports alignment of non-natural peptides, 
containing modified and artificial amino acids.
MSA works both for FASTA format and for HELM with linear sequences .

### Compositional analysis

Makes a “weblogo” representation of letter composition at each position. Usually used on aligned data.

### To atomic level

Reconstructs the chemical structure of peptides/DNAs. Requires presence of a monomer library containing all monomers.

### Sequence activity cliffs

Determines the most sharp sequence-activity changes. This action requires two columns of the dataset – the sequence and activity

### Sequence space

Calculates the clusters of peptides/DNAs based on the sequence. The FASTA/HELM sequences are supported as an input. In most cases to obtain valuable data need to be run on aligned data.

### Similarity search

Creates linked table where for each selected sequence form initial table it displays the most similar sequences (10 by default, the number and cutoff can be changed in table properties)

### Diversity search

Very similar to “Similarity search”, just displays the sequences most diverse form the current.

### Peptides SAR

The bid module analysing the peptides sequence-activity relationships and identifying the most promising changes.
See detailed description here: 

## ToDo

DataGrok can open and render Phylogenetic trees in Newick notation. 
DataGrok has its own viewer “Dendrogram”, specially designed for fast rendering of big trees. The tree about 100000 nodes are rendered in several seconds.
Also DataGrok has built-in connectors to PhyloTree package which is not so fast, but provides some additional functionality.
