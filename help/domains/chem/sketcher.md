<!-- TITLE: Molecule sketcher -->
<!-- SUBTITLE: -->

# Molecule sketcher

Sketch a molecule using the built-in editor, or retrieve one by entering compound
identifier. The following compound identifiers are natively understood since they
have a prefix that uniquely identifies source system:

| Identifier       | Example                                            |
|------------------|----------------------------------------------------|
| SMILES           | CC1=CC(=O)C=CC1=O                                  |
| InChI            | InChI=1S/C7H6O2/c1-5-4-6(8)2-3-7(5)9/h2-4H,1H3     |
| InChIKey         | VTWDKFNVVLAELH-UHFFFAOYSA-N                        |
| CHEMBL           | CHEMBL358225                                       |
| MCULE            | MCULE-9650878933                                   |
| comptox          | DTXSID8060290                                      |
| zinc             | ZINC000000080829                                   |

All of the following sources can be referenced by prefixing source name followed by colon to the 
identifier, i.e. 'pubchem:11122':

chembl, pdb, drugbank, pubchem_dotf, gtopdb, ibm, kegg_ligand, zinc, nih_ncc,
emolecules, atlas, chebi, fdasrs, surechembl, pubchem_tpharma, pubchem, recon,
molport, bindingdb, nikkaji, comptox, lipidmaps, carotenoiddb, metabolights,
brenda, pharmgkb, hmdb, nmrshiftdb2, lincs, chemicalbook, selleck, mcule, actor,
drugcentral, rhea

![Sketcher](../../uploads/chem/sketcher.png "Sketcher")

See also

* [Mapping identifiers](cheminformatics.md)
* [Cheminformatics](cheminformatics.md)
