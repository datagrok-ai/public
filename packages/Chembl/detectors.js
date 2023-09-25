class ChemblPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMolRegNo(col) {
    if (col.name === 'molregno')
      return 'molregno';
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectChembl(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.startsWith('CHEMBL'), 1))
      return 'chembl';
    return null;
  }
}
