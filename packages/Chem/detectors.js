class ChemPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    if ((
      col.name.toLowerCase().endsWith('smiles') ||
      col.name.toLowerCase().endsWith('rdkit') ||
      col.name.toLowerCase().endsWith('molecule') ||
      col.name.toLowerCase().endsWith('scaffold')
    ) && col.type === DG.TYPE.STRING) {
      col.semType = DG.SEMTYPE.MOLECULE;
      return col.semType;
    }

    return null;
  }
}
