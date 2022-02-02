<!-- TITLE: Cheminformatics concepts -->
<!-- SUBTITLE: -->

# Cheminformatics

Please see the
20-minute [cheminformatics introduction video](https://www.youtube.com/watch?v=yM0ums_ur78&ab_channel=JeremyYang-Datascience%2Cetc)
.

_[Wikipedia](https://en.wikipedia.org/wiki/Cheminformatics): Cheminformatics is the use of computer and informational
techniques applied to a range of problems in the field of chemistry. These in silico techniques are used, for example,
in pharmaceutical companies and academic settings in the process of drug discovery. These methods can also be used in
chemical and allied industries in various other forms._

But what does *in silico* actually means? It is counterposed to *in vitro* (undertaken in glass without cells) and *in
vivo* (undertaken in cell cultures or organisms) biochemical experiments. So, *in silico* was defined to imply a
biochemical experiment performed via computational approaches (as silicon is the main component of modern CPU).

Cheminformatics is an *in silico* discipline serving to handle chemical entities for a number of different purposes:

- generate possible chemical structures [all possible](https://gdb.unibe.ch/downloads/) or
  with [synthetic rules implementation](https://cactus.nci.nih.gov/download/savi_download/);
- storing
  compounds [collections](https://www.merckgroup.com/en/research/open-innovation/biopharma-open-innovation-portal/open-compound-sourcing.html)
  ;
- facilitation of chemists' work via visualization and reaction prediction;
- prediction of physical and chemical properties of the compounds (QSPR) and biological activities of the
  compounds ([QSAR](https://en.wikipedia.org/wiki/Quantitative_structure%E2%80%93activity_relationship))
  ;
- [ADMET](https://en.wikipedia.org/wiki/ADME) prediction;
- virtual screening to evaluate the most potent drug candidates.

## Storing chemical information

The most convenient form of molecular representation for human
is [molecular graph](https://en.wikipedia.org/wiki/Molecular_graph), and Datagrok indeed fulfills this convenience to
its users. Check it at [smiles.csv](https://public.datagrok.ai/f/Demo.TestJobs.Files.DemoFiles/chem/smiles.csv),
double-click the molecular graph to modify it. Mathematically, graphs are sets of vertices and edges, and could be
processed by machines in different ways.

To store the whole information about the molecule, the [MOL](http://c4.cabrillo.edu/404/ctfile.pdf)
format is widely used. It contains atoms (vertices) and bonds (edges) of the molecule with all associated information as
atom coordinates, charges, isotopes etc. Multiple MOL files are stored as an SDF where any other additional information
could be added (e.g. experimental activity values). Being very convenient, MOL is also commonly used in the overwhelming
majority of computational cheminformatics software.

The tightest popular format for the storage is [SMILES](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
string. Special rules are used for such string generation. Datagrok uses SMILES to restore MOL files for subsequent
processing. For chemical reactions, the modified strings
of [SMIRKS](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html) format are used.

MOL and SMILES are possibly the best formats for storing a single molecule. However, some applications require the
notion of uncertain chemical structures (e.g. fragments that do not correspond to any real molecule). For such needs,
logical expressions are essential, (e.g. carbon or oxygen atom, double or aromatic bond)
and [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) format provides this option. One may say
that SMARTS is the same thing as regular expressions in cheminformatics. Datagrok uses SMARTS for
finding [structural alerts](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/info-panels/structural-alerts.md)
, performing substructure
search, [R-group analysis](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/r-group-analysis.md) (*
revision required*).

## Operating with chemical information

Though molecular graph is an appropriate form of data for storage it is useless for the wide spectrum of cheminformatics
applications especially machine learning procedures. For such applications, molecules are represented as a set of
molecular [descriptors](https://www.youtube.com/watch?v=0j1Eeexd1ig&ab_channel=ChemAxon) or molecular fingerprints.
Datagrok provides generation of different sets
of [descriptors](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/descriptors.md)
and [fingerprints](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/fingerprints.md)
. The objective of these representations is to satisfy linear algebra requirements applied in the majority of ML methods
and provide a vector of values describing the molecule. Vector (linear) spaces based on molecular descriptors are called
chemical spaces.

Descriptors are frequently used for proceeding [similar](https://en.wikipedia.org/wiki/Chemical_similarity) chemical
structures. These principles yield
in [similarity search](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/similarity-search.md)
and [diversity search](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/diversity-search.md)
. In combination with clustering and self-organizing maps, methods such
as [stochastic proximity embedding](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/similarity-spe.md)
allow to reduce dimensionality of used vectors and to separate the most significant features of the molecule. It helps
us to visualize the chemical space in [2D maps](https://public.datagrok.ai/script/c7e74227-ca12-5337-b7d0-8c4d1bacfbb9).

Pharmaceutical needs demand wide use of cheminformatics methods for chemical datasets exploration analysis and following
modelling studies. These datasets are always accompanied with experimental values (e.g. measured biological activity of
a compound). One of the most common task is evaluation of structure-activity relationships, which are essential in drug
development as they contribute in [hit](https://drugdevelopment.fi/drug-development/hit-to-lead/) compound
identification and lead compound optimization. (Q)SAR studies are performed to find possible leads in the screening
datasets. (*revision required - possible field of Datagrok interests, prediction of solubility, activity cliffs analysis
are currently implemented*)

All described methods are implemented in different analysis pipelines, and assume that descriptors describe a real
molecule perfectly. Data-associated errors lead to biases in descriptors, wrong interpretations of modeling outputs, and
irrelevance of the whole work. The most sensitive cases are duplicated vectors in the training set, and errors derived
from the incorrect structure representation.
Thus, [curation](https://github.com/datagrok-ai/public/blob/master/help/domains/chem/chem-curate.md)
of chemical data is usually integrated to analysis pipeline.

## Freely available tools and public databases

There is a number of freely-available cheminformatics toolkits available for investigators and
developers. [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html)
, [CDK](https://cdk.github.io/), [Open Babel](http://openbabel.org/wiki/Main_Page) are among them. RDKit is widely used
at Datagrok in form of Python, C++ and JS code.

To provide an opportunity of having up-to-date data for its users, Datagrok is connected
to [Pubchem](https://pubchem.ncbi.nlm.nih.gov/), [ChEMBL](https://www.ebi.ac.uk/chembl/) and other databases that
contain continuously updated chemical structures data associated with biological experimental values.

## Cheminformatics limitations

Though it might seem that cheminformatics covers all molecular needs, it has its limits. The applicability of this
discipline is limited to small molecules and some types of peptides. It shows real power when it is combined with
other *in silico* approaches including but not limited to:

- docking and molecular dynamics
- systems biology and pharmacology
- bioinformatics
- proteomics.

See also:

- [Cheminformatics in Datagrok](https://datagrok.ai/cheminformatics)
- [Cheminformatics API](cheminformatics.md)
