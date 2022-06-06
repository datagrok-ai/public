class ChemPackageDetectors extends DG.Package {
  static sValidSmilesChar = new Uint32Array([0, 671083320, 1073741823, 134217726, 0, 0, 0, 0]);

  static sValidFirstSmilesChar = new Uint32Array([0, 256, 134857308, 573448, 0, 0, 0, 0]);

  static validChar(char) {
    const pos = char.charCodeAt(0);
    return (ChemPackageDetectors.sValidSmilesChar[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  static validFstChar(char) {
    const pos = char.charCodeAt(0);
    return (ChemPackageDetectors.sValidFirstSmilesChar[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    // return DG.SEMTYPE.MOLECULE;
    if (col.toList().every((el) => el && el.length > 1 &&
        ChemPackageDetectors.validFstChar(el[0]) &&
        el.split('').every((sChar) => ChemPackageDetectors.validChar(sChar)))) {
      grok.functions.call('Chem:checkSmilesValidity', {col: col}).then((validPerc) => {
        if (validPerc > 0.79) {
          col.semType = DG.SEMTYPE.MOLECULE;
          // smiles or molblock?
          const str = col.length > 0 ? col.get(0) : null;
          col.tags[DG.TAGS.UNITS] = (str !== null && str.includes('M  END')) ?
            DG.UNITS.Molecule.MOLBLOCK : DG.UNITS.Molecule.SMILES;
        }
      });
    }
  }
}
