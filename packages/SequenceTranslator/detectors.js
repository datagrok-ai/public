class SequenceTranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (DG.Detector.sampleCategories(col, (s) => isDnaNucleotides(s)))
        return 'DNA nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => isRnaNucleotides(s)))
        return 'RNA nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => isAsoGapmerBioSpring(s)))
        return 'BioSpring / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => isAsoGapmerGcrs(s)))
        return 'GCRS / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => isSiRnaBioSpring(s)))
        return 'BioSpring / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => isSiRnaAxolabs(s)))
        return 'Axolabs / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => isGcrs(s)))
        return 'GCRS';
      if (DG.Detector.sampleCategories(col, (s) => isMermade12(s)))
        return 'Mermade 12 / siRNA';
    }
  }
}
