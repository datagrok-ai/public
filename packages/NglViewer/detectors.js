class NglViewerPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPdbId(col) {
    let name = col.name.toLowerCase();
    if (col.type === DG.TYPE.STRING && name.contains('pdb') && DG.Detector.sampleCategories(col,(s) => s.length === 4))
      return 'pdb_id';
  }
}
