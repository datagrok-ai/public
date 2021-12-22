class PeptidesPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {
    const regexp = new RegExp(/^([^-^\n]*-){7,49}(\w|\(|\))+$/);
    return DG.Detector.sampleCategories(col, (s) => regexp.test(s.trim())) ? 'alignedSequence' : null;
  }
}
