class ChemspacePackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectCsId(col) {
    if (['CS-id', 'csId'].includes(col.name)) {
      col.semType = 'chemspace-id';
      return col.semType;
    }

    return null;
  }
}
