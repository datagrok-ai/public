class PeptidesPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {
    const regexp = new RegExp(/^([^-^\n]*-){2,49}(\w|\(|\))+$/);
    return DG.Detector.sampleCategories(col, (s) => {
      let res = regexp.test(s.trim());
      if (!res) {
        console.error(`Failed on ${s}`);
      }
      return res;
    }) ? 'alignedSequence' : null;
  }
}
