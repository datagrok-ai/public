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
    const chembl = Array.from(col.values()).filter((val) => val.startsWith('CHEMBL'));
    if (chembl.length == col.length)
      return 'chembl';
    return null;
  }
}
