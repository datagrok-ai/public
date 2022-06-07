class ChemPackageDetectors extends DG.Package {
  static sValidSmilesChar = new Uint32Array([0, 671083320, 1073741823, 134217726, 0, 0, 0, 0]);

  static sValidFirstSmilesChar = new Uint32Array([0, 256, 134857308, 573448, 0, 0, 0, 0]);

  static numberOfMolsToCheck = 10;

  static validChar(char) {
    const pos = char.charCodeAt(0);
    return (ChemPackageDetectors.sValidSmilesChar[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  static validFirstChar(char) {
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
    for (let i = 0; i < col.length; ++i) {
      const el = col.get(i);
      if (el && el.length > 0) {
        if (!ChemPackageDetectors.validFirstChar(el[0]))
          return null;
        for (let j = 1; j < el.length; ++j)
          if (!ChemPackageDetectors.validChar(el[j]))
            return null;
      }
    }

    let molStrings = [];
    if (col.length <= ChemPackageDetectors.numberOfMolsToCheck)
      molStrings = col.toList().filter(el => el);
    else {
      const inc = Math.floor(col.length / ChemPackageDetectors.numberOfMolsToCheck);
      for (let i = 0; i < col.length; i += inc) {
        let el = col.get(i);
        if (el)
          molStrings.push(el);
        else {
          for (let j = 0; j < Math.min(inc, 5) && i + j < col.length; ++j) {
            el = col.get(i + j);
            if (el) {
              molStrings.push(el);
              break;
            }
          }
        }
      }
    }

    grok.functions.call('Chem:checkSmilesValidity', {molStrings: molStrings}).then((validPerc) => {
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
