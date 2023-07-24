---
title: "Chemical scripts"
---

| Name                               | Function                  |
|------------------------------------|---------------------------|
| Substructure search                |<br /><pre>`\#{x.ChemSubstructureSearch}`</pre>|
| Find MCS                           |<br /><pre>`\#{x.ChemFindMCS}`</pre> |
| Descriptors                        |<br /><pre>`\#{x.ChemDescriptors}`</pre> |
| R-Groups                           |<br /><pre>`\#{x.ChemGetRGroups}`</pre> |
| Fingerprints                       |<br /><pre>`\#{x.ChemFingerprints}`</pre> |
| Similarity SPE                     |<br /><pre>`\#{x.ChemSimilaritySPE}`</pre> |
| SMILES to InchI                    |<br /><pre>`\#{x.ChemSmilesToInchi}`</pre> |
| SMILES to Canonical                |<br /><pre>`\#{x.ChemSmilesToCanonical}`</pre> |
| Chemical map identifiers           |<br /><pre>`\#{x.ChemMapIdentifiers}`</pre> |
| Butina cluster                     |<br /><pre>`\#{x.ChemScripts:ButinaMoleculesClustering}`</pre> |
| [Filter by catalogs](filter-catalogs.md)                 |<br /><pre>`\#{x.ChemScripts:FilterByCatalogs}`</pre> |
| [Gasteiger partial charges](gasteiger-charges.md)          |<br /><pre>`\#{x.ChemScripts:GasteigerPartialCharges}`</pre> |
| [Murcko scaffolds](murcko-scaffolds.md)                   |<br /><pre>`\#{x.ChemScripts:MurckoScaffolds}`</pre>|
| [Similarity maps using fingerprints](sim-maps.md) |<br /><pre>`\#{x.ChemScripts:SimilarityMapsUsingFingerprints}`</pre> |
| [Chemical space using tSNE](tsne.md)          |<br /><pre>`\#{x.ChemScripts:ChemicalSpaceUsingtSNE}`</pre> |
| [Two component reactions](reactions.md)           |<br /><pre>`Chem:TwoComponentReaction`</pre> |
| [Chemical space using UMAP](umap.md)          |<br /><pre>`\#{x.ChemScripts:ChemicalSpaceUsingUMAP}`</pre> |
| [USRCAT](usrcat.md)                |<br /><pre>`\#{x.ChemScripts:USRCAT}`</pre> |
| [Mutate](mutate.md)                            ||
| [Solubility prediction](solubility-prediction.md)             |<br /><pre>`\#{x.18b704d0-0b50-11e9-b846-1fa94a4da5d1."Predict Solubility"}`</pre>|
| Curate                             ||

The following table gives an indicative data for the performance of certain chemical functions:

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
