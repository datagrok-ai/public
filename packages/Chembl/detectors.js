class ChemblPackageDetectors extends DG.Package {

  //meta.role: semValueExtractor
  //input: string s
  //output: semantic_value result
  detectChemblId(s) {
    if (/^CHEMBL\d+$/.test(s))
      return DG.SemanticValue.fromValueType(s, 'CHEMBL_ID');
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectMolRegNo(col) {
    if (col.name === 'molregno')
      return 'molregno';
    return null;
  }
}
