const RAW_PNG_SEM_TYPE = 'BinaryImage';

class DemoPackageDetectors extends DG.Package {
  static likelyNames = ['smiles', 'mol'];

  //name: detectMolecules
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectMolecules(col) {
    if (col.type === DG.TYPE.STRING && DemoPackageDetectors.likelyNames.includes(col.name)) {
      col.semType = DG.SEMTYPE.MOLECULE;
      return DG.SEMTYPE.MOLECULE;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //meta.semType: Macromolecule
  //meta.skipTest: GROK-1
  //input: column col
  //output: string semType
  detectSequences(col) {
    return grok.functions.call('Demo:detectMolecules', {col}) ? 'Macromolecule' : null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  static async detectImages(col) {
    if (col.name === 'png') {
      col.semType = RAW_PNG_SEM_TYPE;
      return col.semType;
    }
    if (col.name === 'txt')
      return 'Text';
    return grok.functions.call('Demo:detectMolecules', {col});
  }
}
