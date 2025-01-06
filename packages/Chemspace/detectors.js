class ChemspacePackageDetectors extends DG.Package {

  //meta.role: semValueExtractor
  //input: string s
  //output: semantic_value result
  detectCsId(s) {
    if (/^(CSMS|CSMB|CSCS|CSSB|CSSS)\d{11}$/.test(s))
      return DG.SemanticValue.fromValueType(s, 'chemspace-id');
  }
}
