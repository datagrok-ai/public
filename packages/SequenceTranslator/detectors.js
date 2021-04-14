class SequencetranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (DG.Detector.sampleCategories(col, (s) => /^[ATGCU]{10,}$/.test(s)))
        return 'nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => /^[*56789ATGC]{30,}$/.test(s)))
        return 'BioSpring / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => /^(?=.*moe)(?=.*5mC)(?=.*ps){30,}/.test(s)))
        return 'GCRS / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => /^[*1-8]{30,}$/.test(s)))
        return 'BioSpring / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => /^[fsACGUacgu]{30,}$/.test(s)))
        return 'Axolabs / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => /^[fmpsACGU]{30,}$/.test(s)))
        return 'GCRS / siRNA';
    }
  }
}