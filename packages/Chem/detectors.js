class ChemPackageDetectors extends DG.Package {

  static likelyNames = [
    'structure', 'mol', 'molecule', 'smiles', 'rdkit',
    'canonical_smiles', 'core', 'scaffold',
    'r1', 'r2', 'r3', 'r4', 'r5'];

  static validSmilesChars = new Uint32Array([0, 671083320, 1073741823, 134217726, 0, 0, 0, 0]);
  static validFirstSmilesChars = new Uint32Array([0, 256, 134857308, 573448, 0, 0, 0, 0]);
  static numberOfMolsToCheck = 100;
  static numberOfMolsToCheckRdKit = 10;

  static validChar(pos) {
    return (ChemPackageDetectors.validSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  static validFirstChar(pos) {
    return (ChemPackageDetectors.validFirstSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  static likelyValidSmiles(s) {
    if (!ChemPackageDetectors.validFirstChar(s[0].charCodeAt(0)))
      return false;
    for (let j = 1; j < s.length; ++j)
      if (!ChemPackageDetectors.validChar(s[j].charCodeAt(0)))
        return false;
    return true;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMolBlocks(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.includes('M  END'))) {
      col.tags[DG.TAGS.UNITS] = DG.UNITS.Molecule.MOLBLOCK;
      return DG.SEMTYPE.MOLECULE;
    }
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    console.log('detecting');
    let lowerCaseName = col.name.toLowerCase();
    let likelyMolName = ChemPackageDetectors.likelyNames.some((likelyName) => lowerCaseName.endsWith(likelyName));
    let minUnique = likelyMolName ? 1 : 3;
    let longest = grok_Column_Aggregate(col.dart, 'longest');
    if (!likelyMolName && longest.length < 5)
      return null;

    // temporary workaround to make it understand r-groups like "CC[*:2]"
    if (likelyMolName) {
      grok.functions.call('Chem:detectSmiles', { col: col, min: minUnique }).then(() => {});
      return null;
    }

    if (DG.Detector.sampleCategories(col, ChemPackageDetectors.likelyValidSmiles, minUnique, 10, 0.8))
      grok.functions.call('Chem:detectSmiles', { col: col, min: minUnique }).then(() => {});
  }
}
