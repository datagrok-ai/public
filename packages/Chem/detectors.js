class ChemPackageDetectors extends DG.Package {

  static likelyNames = [
    'structure', 'mol', 'molecule', 'smiles', 'rdkit',
    'canonical_smiles', 'core', 'scaffold',
    'r1', 'r2', 'r3', 'r4', 'r5'];

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;

    let name = col.name.toLowerCase();
    if (ChemPackageDetectors.likelyNames.some((likelyName) => name.endsWith(likelyName))) {
      col.semType = DG.SEMTYPE.MOLECULE;

      // smiles or molblock?
      let str = col.length > 0 ? col.get(0) : null;
      col.tags[DG.TAGS.UNITS] = (str !== null && str.includes('M  END')) ? DG.UNITS.Molecule.MOLBLOCK : DG.UNITS.Molecule.SMILES;

      return col.semType;
    }

    return null;
  }
}
