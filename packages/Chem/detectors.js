class ChemPackageDetectors extends DG.Package {
  static sValidSmilesChar = Uint32Array.from(
    Buffer.from('AAAAABzX/+T////8f///4AAAAAAAAAAAAAAAAAAAAAA', 'base64'));

  static sValidFirstSmilesChar = Uint32Array.from(
    Buffer.from('AAAAAACAAAA6Q5AQEAMQAAAAAAAAAAAAAAAAAAAAAAA', 'base64'));

  static validChar(char) {
    const pos = char.charCodeAt(0);
    return (this.sValidSmilesChar[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  static validFstChar(char) {
    const pos = char.charCodeAt(0);
    return (this.sValidFirstSmilesChar[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRDSmiles(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    console.log('hmmm');
    if (col.toList().every((el) => el && el.length > 1 &&
      ChemPackageDetectors.validFstChar[el[0]] &&
      el.split('').every((sChar) => ChemPackageDetectors.validChar[sChar]))) {
      col.semType = DG.SEMTYPE.MOLECULE;

      // smiles or molblock?
      const str = col.length > 0 ? col.get(0) : null;
      col.tags[DG.TAGS.UNITS] = (str !== null && str.includes('M  END')) ?
        DG.UNITS.Molecule.MOLBLOCK : DG.UNITS.Molecule.SMILES;

      return col.semType;
    }
  }
}
