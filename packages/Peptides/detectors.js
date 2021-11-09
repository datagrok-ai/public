class PeptidesPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {
    const regexp = new RegExp('^((.+)?-){5,49}(\\w|\\(|\\))+$');
    return DG.Detector.sampleCategories(col, (s) => s.length < 1000 && regexp.test(s)) ? 'alignedSequence' : null;
  }
}
