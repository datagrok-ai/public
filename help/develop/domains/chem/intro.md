<!-- TITLE: Cheminformatics concepts -->
<!-- SUBTITLE: -->

1. [Cheminformatics Introduction](https://www.youtube.com/watch?v=yM0ums_ur78&ab_channel=JeremyYang-Datascience%2Cetc.)

_[Wikipedia](https://en.wikipedia.org/wiki/Cheminformatics): Cheminformatics is the use of 
computer and informational techniques applied to a range of 
problems in the field of chemistry. These in silico techniques are used, for example, in pharmaceutical companies 
and academic settings in the process of drug discovery. These methods can also be used in chemical and allied 
industries in various other forms._

But what does *in silico* actually means? It is conterposed to *in vitro* (undertaken in glass without cells) and *in vivo* (undertaken in cell cultures or organisms) biolchemical experiments. So *in silico* was defined to imply a biolchemical experiment performed via computational approaches (as silicon is the main component of modern CPU). 

Cheminformatics is an *in silico* discipline serving to handle chemical entities for a number of different purposes:
- generate possible chemical structures [all possible](https://gdb.unibe.ch/downloads/) or with [synthetic rules implementation](https://cactus.nci.nih.gov/download/savi_download/);
- storing compounds [collections](https://www.merckgroup.com/en/research/open-innovation/biopharma-open-innovation-portal/open-compound-sourcing.html);
- facilitation of chemists' work via visualization and reaction prediction;
- prediction of physical and chemical properties of the compounds (QSPR) and biological activities of the compounds ([QSAR](https://en.wikipedia.org/wiki/Quantitative_structure%E2%80%93activity_relationship));
- [ADMET](https://en.wikipedia.org/wiki/ADME) prediction;
- virtual screening to evaluate the most potent drug candidates.

2. Storing chemical information

As one may know the most convinient form of molecular representation for human is [molecular graph](https://en.wikipedia.org/wiki/Molecular_graph) and Datagrok indeed fulfills this convinience to its users, check it at https://public.datagrok.ai/f/Demo.TestJobs.Files.DemoFiles/chem/smiles.csv (double click the molecular graph to modify it). But mathematically graphs are sets of vertices and edges and could be processed by machines in different ways.

3. Operating with chemical information


See also:
* [Cheminformatics API](../../cheminformatics.md)