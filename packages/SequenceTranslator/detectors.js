class SequencetranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (DG.Detector.sampleCategories(col, (s) => /^[ATGC]*$/.test(s)))
        return 'nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => /^[*56789ATGC]*$/.test(s)))
        return 'BioSpring / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => /^(moe|5mC|n|ps|A|C|G|T|U)*$/.test(s)))
        return 'GCRS / Gapmers';
    }
  }
}