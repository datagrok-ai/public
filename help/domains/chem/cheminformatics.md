<!-- TITLE: Cheminformatics -->
<!-- SUBTITLE: -->

# Cheminformatics

_[Wikipedia](https://en.wikipedia.org/wiki/Cheminformatics): Cheminformatics is the use of 
computer and informational techniques applied to a range of 
problems in the field of chemistry. These in silico techniques are used, for example, in pharmaceutical companies 
and academic settings in the process of drug discovery. These methods can also be used in chemical and allied 
industries in various other forms._

Datagrok provides first-class support for small molecules, as well as most popular building blocks for 
cheminformatics. It understand several popular notations for representing chemical (sub)structures, 
such as SMILES and SMARTS. Molecules can be rendered in either 2D or 3D with different visualization options.
They can be sketched as well. Chemical properties, descriptors, and fingerprints can be extracted.
Predictive models that accept molecules as an input can be easily trained, assessed, executed, deployed, reused by other scientists, and used in pipelines or in [info panels](../../discover/info-panels.md).
 
Several toxicity and drug-likeness prediction models are supported. Substructure and similarity search
works out-of-the box for imported data, and can be efficiently utilized for querying databases using 
Postgres chemical cartridge. To further explore collections
of molecules, use advanced tools such as [diversity search](diversity-search.md) and 
[similarity search](similarity-search.md).  

## Importing molecular data

Simply [import the dataset](../../access/importing-data.md) as you normally would - by opening a file, 
querying a database, connecting to a webservice, or by any other method. The platform is smart enough 
to automatically recognize chemical structures. 

## Molecule sketcher

Sketch a molecule using the built-in editor, or retrieve one by entering a compound
identifier. The following compound identifiers are natively understood since they
have a prefix that uniquely identifies source system: SMILES, InChI, InChIKey,
CHEMBL, MCULE, comptox, and zinc. The rest of the 30+ identifier systems can be 
referenced by prefixing source name followed by colon to the identifier, i.e. 'pubchem:11122'.

![Sketcher](../../uploads/chem/sketcher.png "Sketcher")

## Chemically-aware viewers

Many viewers, such as 
[grid](chemically-aware-viewers.md#grid),
[scatter plot](chemically-aware-viewers.md#scatter-plot),
[network diagram](chemically-aware-viewers.md#network-diagram),
[tile viewer](chemically-aware-viewers.md#tile-viewer),
[bar chart](chemically-aware-viewers.md#bar-chart), 
[form viewer](chemically-aware-viewers.md#form-viewer), and
[trellis plot](chemically-aware-viewers.md#trellis-plot)
 will recognize and render chemical structures.   

## Accessing cheminformatics tools

Chemical intelligence tools are natively integrated into the platform, so in most cases the 
appropriate functionality is automatically presented based on the user actions and context.
For instance, when user clicks on a molecule, it becomes a [current object](../../overview/navigation.md#current-object), 
and its properties are shown in the [property panel](../../overview/navigation.md#properties). To see chemically-related actions
applicable for the specified column, right-click on the column, and look under
`Current column | Chem` and `Current column | Extract`. Alternatively, click on the column of 
interest, and expand the 'Actions' section in the property panel. 

Check out 'Tools | Chemistry' to see additional functionality.

As always, it is a good idea to search for functionality using the smart search (Alt+Q), or
by opening the registry of available functions `Help | Functions`.   

## Chemical properties and descriptors

Use 'Extract' popup menu to calculate the following properties: 
formula, drug likeness, acceptor count, donor count, logP, logD, polar surface area, 
rotatable bond count, stereo center count. 

Chemical descriptors are numerical features extracted from chemical structures for molecular 
data mining, compound diversity analysis and compound activity prediction. In addition to properties, 
the platform also makes it easy to compute different sets of [molecular descriptors](descriptors.md). 
Supported descriptor sets are: Lipinski, Crippen, EState, EState VSA, Fragments, Graph, MolSurf, QED. 

## Fingerprints

[Fingerprints](fingerprints.md) are a very abstract representation of certain structural features 
of a molecule. Similarity measures, calculations that quantify the similarity of two molecules, 
and screening, a way of rapidly eliminating molecules as candidates in a substructure search, 
are both processes that use fingerprints. Grok supports the following fingerprints: RDKFingerprint,
MACCSKeys, AtomPair, TopologicalTorsion, Morgan/Circular.

## Similarity / diversity analyses

We have implemented few tools that help scientists analyze a collection of molecules in 
terms of molecular similarity. Both tools are based on applying different distance metrics 
(such as Tanimoto) to fingerprints.

* [Similarity Search](similarity-search.md) - finds structures similar to the specified one 
* [Diversity Search](diversity-search.md) - finds 10 most distinct molecules

These tools can be used together as a collection browser. 'Diverse structures' window shows different classes
of compounds present in the dataset; when you click on a molecule representing a class, similar
molecules will be shown in the 'Similar structures' window.

## Clustering

* [Similarity analysis using Stochastic Proximity Embedding](similarity-spe.md)

## Physical predictive models

* Toxicity - predicts the following toxicity properties: mutagenicity, tumorigenicity, 
irritating effects, reproductive effects.
* Drug likeness - a score that shows how likely this molecule is to be a drug. The score comes
with an interpretation of how different sub-structure fragments contribute to the score.

## Machine learning predictive models

In contrast to the physical predictive models, machine learning predictive models do not have
any intrinsic knowledge about the physical and biological processes. Instead, they use techniques
such as random forests or deep learning to discern mathematical relationships between empirical 
observations of small molecules and extrapolate them to predict chemical, biological and physical 
properties of novel compounds.

Datagrok enables machine learning predictive models by using chemical 
[properties, descriptors](#chemical-properties-and-descriptors), 
and [fingerprints](#fingerprints) as features, and the observed properties as results when 
[building predictive models](../../learn/data-science.md#predictive-modeling). This lets scientists
build predictive models that can be trained, assessed, executed, reused by other scientists, and used in pipelines.   

See [Cheminformatics predictive modeling](chem-predictive-modeling.md) for
more details and a demo of building and applying a model.

References:
* [Machine learning in chemoinformatics and drug discovery](https://www.sciencedirect.com/science/article/pii/S1359644617304695)

## Molecule identifier conversions

Grok lets users easily and efficiently convert molecule identifiers between different 
source systems, including proprietary company identifiers.

Supported sources are: 
chembl, pdb, drugbank, pubchem_dotf, gtopdb, ibm, kegg_ligand, zinc, nih_ncc,
emolecules, atlas, chebi, fdasrs, surechembl, pubchem_tpharma, pubchem, recon,
molport, bindingdb, nikkaji, comptox, lipidmaps, carotenoiddb, metabolights,
brenda, pharmgkb, hmdb, nmrshiftdb2, lincs, chemicalbook, selleck, mcule, actor,
drugcentral, rhea 

To map the whole column containing identifiers, use #{x.ChemMapIdentifiers} function.

IUPAC name is located in the "Properties" panel.

In order to retrieve a single structure by an identifier, it might be handy to use
[Sketcher](sketcher.md)

## Info panels  

Click on a molecule to select it as a current object. This will bring up
this molecule's properties to the property panel. The following panels are part of the 
'chem' plugin: 

* Structure - 2D structure
* Properties - all above-mentioned properties
* SDF - molfile
* 3D - interactive 3D rendering
* [Toxicity](info-panels/toxicity-risks.md) - results of the toxicity prediction
* [Drug likeness](info-panels/drug-likeness.md) - a score that shows how likely this molecule is to be a drug. The score comes
with an interpretation of how different sub-structure fragments contribute to the score.
* Identifiers - all known identifiers for the specified structure ([UniChem](https://www.ebi.ac.uk/unichem/ucquery/listSources)) 
* Patents - patents associated with that structure ([SureChEMBL](https://www.surechembl.org)) 
* ChEMBL similar structures 
* ChEMBL substructure search 
* Gasteiger Partial Charges visualization 
* [Structural alerts](info-panels/structural-alerts.md)

![](../../uploads/gifs/chem-model-augment.gif "Toxicity, Gasteiger Partial Charges, Solubility Prediction")

In addition to these pre-defined info panels, users can develop their own using any
[scripting language](../../develop/scripting.md) supported by the Grok platform. For example, 
#{x.demo:pythonscripts:GasteigerPartialCharges}. 
 
## In-memory substructure search

To search for molecules within the table that contain specified substructure, click on the molecule column,
and press Ctrl+F. To add a substructure filter to [column filters](../../visualize/viewers/filters.md), click on the
'☰' icon on top of the filters, and select the molecular column under the 'Add column filter' submenu.

## Most common substructure

The maximum common substructure (MCS) problem is of great importance in multiple aspects of cheminformatics. 
It has diverse applications ranging from lead prediction to automated reaction mapping and visual alignment of 
similar compounds.

To find MCS for the column with molecules, run `Chem | Find MCS` command from column's context menu. To execute 
it from the [console](), use `chem:findMCS(tableName, columnName)` command. 

## R-Group analysis

[R-Group Analysis](r-group-analysis.md) is a common function in chemistry. Typically, it involves R-group decomposition, 
followed by the visual analysis of the obtained R-groups. Grok's chemically-aware 
[Trellis Plot](../../visualize/viewers/trellis-plot.md) is a natural fit for such an analysis. 

![R-Group Analysis](../../uploads/chem/r-group-analysis.png "R-Group Analysis"){:height="100px" width="60px"}

## Scaffolds

The scaffold concept is widely applied in medicinal chemistry. Scaffolds are mostly used to 
represent core structures of bioactive compounds. Although the scaffold concept has limitations 
and is often viewed differently from a chemical and computational perspective, it has provided 
a basis for systematic investigations of molecular cores and building blocks, going far beyond 
the consideration of individual compound series.

[Murcko Scaffolds](functions/murcko-scaffolds.md)

## Chemical reactions

Applies a specified reaction to two columns containing molecules. The output table contains a row 
for each product produced by applying the reaction to the inputs. Each row contains the product 
molecule, index information, and the reactant molecules that were used. 

'Do Matrix Expansion': If checked, each reactant 1 will be combined with each reactant 2 yielding 
the combinatorial expansion of the reactants. If not checked, reactants 1 and 2 will be combined 
sequentially, with the longer list determining the number of output rows.

Corresponding function: #{x.demo:pythonscripts:TwoComponentReaction} 

See details [here](functions/reactions.md).

## Functions

The following cheminformatics-related [functions](../../overview/functions/function.md) are exposed:

* #{x.ChemSubstructureSearch}
* #{x.ChemFindMCS}
* #{x.ChemDescriptors}
* #{x.ChemGetRGroups}
* #{x.ChemFingerprints}
* #{x.ChemSimilaritySPE}
* #{x.ChemSmilesToInchi}
* #{x.ChemSmilesToCanonical}
* #{x.ChemMapIdentifiers}

Lot of chemical analysis is implemented using [scripting](../../develop/scripting.md) functionality:
* #{x.ChemScripts:ButinaMoleculesClustering}
* #{x.ChemScripts:FilterByCatalogs}
* #{x.ChemScripts:GasteigerPartialCharges}
* #{x.ChemScripts:MurckoScaffolds}
* #{x.ChemScripts:SaltStripper}
* #{x.ChemScripts:SimilarityMapsUsingFingerprints}
* #{x.ChemScripts:ChemicalSpaceUsingtSNE}
* #{x.ChemScripts:TwoComponentReaction}
* #{x.ChemScripts:ChemicalSpaceUsingUMAP}
* #{x.ChemScripts:USRCAT}

# Performance

| Function                              | Molecules | Execution time, s |
|---------------------------------------|-----------|-------------------|
| ChemSubstructureSearch                | 1M        | 70                |
| ChemFindMcs                           | 100k      | 43                |
| ChemDescriptors (201 descriptor)      | 1k        | 81                |
| ChemDescriptors (Lipinski)            | 1M        | 164               |
| ChemGetRGroups                        | 1M        | 233               |
| ChemFingerprints (TopologicalTorsion) | 1M        | 782               |
| ChemFingerprints (MACCSKeys)          | 1M        | 770               |
| ChemFingerprints (Morgan/Circular)    | 1M        | 737               |
| ChemFingerprints (RDKFingerprint)     | 1M        | 2421              |
| ChemFingerprints (AtomPair)           | 1M        | 1574              |
| ChemSmilesToInChI                     | 1M        | 946               |
| ChemSmilesToInChIKey                  | 1M        | 389               |
| ChemSmilesToCanonical                 | 1M        | 331               |


## Database substructure and similarity search

Efficient substructure and similarity searching in a database containing information about molecules is a key requirement for 
any chemical information management system. This is typically done by installing a so-called chemical cartridge
on top of a database server. The cartridge extends server's functionality with the molecule-specific operations,
which are made efficient by using chemically-aware indexes, which are often based on molecular fingerprints. Typically,
these operations are functions that can be used as part of the SQL query.

Datagrok provides mechanisms for the automated translation of queries into SQL statements for several commonly
used chemical cartridges. We support the following ones:     

1) [RDKit Postgres cartridge](https://www.rdkit.org/docs/Cartridge.html)   
2) [JChem cartridge](https://docs.chemaxon.com/display/docs/JChem+Cartridge)  (todo)

![DB Substructure and Similarity Search](../../uploads/gifs/db-substructure-similarity-search.gif "DB Substructure and Similarity Search")

See [DB Substructure and similarity search](db-substructure-similarity-search.md) for details.

## Public datasets deployed on our servers

* [ChEMBL](ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/) (Postgres)
* [UniChem](ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/) (Postgres)
* TODO: cheminformatics training/demo datasets

### Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/k1NVdTRpYOM" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also:

* [Descriptors](descriptors.md)
* [Diversity Search](diversity-search.md)
* [Similarity Search](similarity-search.md)
* [Fingerprints](fingerprints.md)
* [Similarity SPE](similarity-spe.md)
* [Info Panels](../../discover/info-panels.md)
* [GrokCompute](../../develop/admin/grok-compute.md)

References:
* [Machine learning in chemoinformatics and drug discovery](https://www.sciencedirect.com/science/article/pii/S1359644617304695)
