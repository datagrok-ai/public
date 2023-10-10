class ChemblPackageDetectors extends DG.Package {

  //meta.role: semValueExtractor
  //input: string s
  //output: semantic_value result
  detectChemblId(s) {
    if (s.startsWith('CHEMBL'))
      return DG.SemanticValue.fromValueType(s, 'chembl');
  }

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
