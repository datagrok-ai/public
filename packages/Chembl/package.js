class ChemblPackage extends DG.Package {
  //name: chemblIdToSmilesTest
  //input: string id
  //output: string smiles { semType: Molecule }
  //meta.role: converter
  //meta.inputRegexp: CHEMBL[0-9]+
  chemblIdToSmiles(id) {
    return 'COCCN=C(S)N(Cc1cccnc1)Cc1ccc(OC)c(OC)c1';
  }
}