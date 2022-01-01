class SequenceTranslatorPackageDetectors extends DG.Package {
  isDnaNucleotides(sequence) {return /^[ATGC]{6,}$/.test(sequence);}
  isRnaNucleotides(sequence) {return /^[AUGC]{6,}$/.test(sequence);}
  isAsoGapmerBioSpring(sequence) {return /^[*56789ATGC]{6,}$/.test(sequence);}
  isAsoGapmerGcrs(sequence) {return /^(?=.*moe)(?=.*5mC)(?=.*ps){6,}/.test(sequence);}
  isSiRnaBioSpring(sequence) {return /^[*1-8]{6,}$/.test(sequence);}
  isSiRnaAxolabs(sequence) {return /^[fsACGUacgu]{6,}$/.test(sequence);}
  isSiRnaGcrs(sequence) {return /^[fmpsACGU]{6,}$/.test(sequence);} // TODO: insert into detectNucleotides
  isGcrs(sequence) {return /^[fmpsACGU]{6,}$/.test(sequence);}
  isMermade12(sequence) {return /^[IiJjKkLlEeFfGgHhQq]{6,}$/.test(sequence);}

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