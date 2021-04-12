class SequencetranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      const checksNum = Math.max(1, Math.round(0.1 * col.length));
      if (DG.Detector.sampleCategories(col, (s) => /^[ATGC]*$/.test(s), 1, checksNum))
        return 'dna_sequence/classic code';
      if (DG.Detector.sampleCategories(col, (s) => /^[*56789ATGC]*$/.test(s), 1, checksNum))
        return 'BioSpring Code For ASO Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => /^(moe|5mC|n|ps|A|C|G|T|U)*$/.test(s), 1, checksNum))
        return 'Janssen GCRS code For ASO Gapmers';
    }
  }
}