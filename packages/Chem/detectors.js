class ChemPackageDetectors extends DG.Package {

  static likelyNames = [
    'structure', 'mol', 'molecule', 'smiles', 'rdkit',
    'canonical_smiles', 'core', 'scaffold'];

  // typical r-group names, such as "R1" or "R100-family"
  static likelyRegexps = [/[r|R][0-9]+(\W|$)/];

  static likelyChemicalName(s) {
    return ChemPackageDetectors.likelyNames.some((likelyName) => s.includes(likelyName)) ||
      ChemPackageDetectors.likelyRegexps.some((regexp) => regexp.test(s));
  }

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

  //name: detectMolecules
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //meta.skipTest: GROK-17630
  detectMolecules(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.includes('M  END'), 1)) {
      col.meta.units = DG.UNITS.Molecule.MOLBLOCK;
      return DG.SEMTYPE.MOLECULE;
    }

    let lowerCaseName = col.name.toLowerCase();
    let likelyMolName = ChemPackageDetectors.likelyChemicalName(lowerCaseName);
    let minUnique = likelyMolName ? 1 : 3;
    let longest = '';
    try {
      longest = col.aggregate('longest') ?? '';
    } catch(x) {}
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
