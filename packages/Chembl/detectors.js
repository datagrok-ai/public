const { chem } = require("datagrok-api/dg");

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
    const chembl = col.getRawData().map((value) => value.startsWith('CHEMBL'));
    if (chembl.length == col.length)
      return 'chemblId';
    return null;
  }
}
