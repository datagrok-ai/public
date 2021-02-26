class ChemPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    if ((
      col.name === 'smiles' ||
      col.name === 'rdkit' ||
      col.name === 'molecule' ||
      col.name === 'scaffold' ||
      col.name === 'Scaffold'
    ) && col.type === DG.TYPE.STRING) {
      col.semType = DG.SEMTYPE.MOLECULE;
      return col.semType;
    }

    return null;
  }
}
