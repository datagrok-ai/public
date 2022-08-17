class SequenceTranslatorPackageDetectors extends DG.Package {
  isDnaNucleotides(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|A|T|G|C){6,}$/.test(sequence);}
  isRnaNucleotides(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C){6,}$/.test(sequence);}
  isAsoGapmerBioSpring(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|5|6|7|8|9|A|T|G|C){6,}$/.test(sequence);}
  isAsoGapmerGcrs(sequence) {return /^(?=.*moe)(?=.*5mC)(?=.*ps){6,}/.test(sequence);}
  isSiRnaBioSpring(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|1|2|3|4|5|6|7|8){6,}$/.test(sequence);}
  isSiRnaAxolabs(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|s|A|C|G|U|a|c|g|u){6,}$/.test(sequence);}
  isSiRnaGcrs(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|m|p|s|A|C|G|U){6,}$/.test(sequence);}
  isGcrs(sequence) {return /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|m|p|s|A|C|G|U){6,}$/.test(sequence);}
  isMermade12(sequence) {
    return /(\(invabasic\)|\(GalNAc-2-JNJ\)|I|i|J|j|K|k|L|l|E|e|F|f|G|g|H|h|Q|q){6,}$/.test(sequence);
  }
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (DG.Detector.sampleCategories(col, (s) => this.isDnaNucleotides(s))) return 'DNA nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => this.isRnaNucleotides(s))) return 'RNA nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => this.isAsoGapmerBioSpring(s))) return 'BioSpring / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => this.isAsoGapmerGcrs(s))) return 'GCRS / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => this.isSiRnaBioSpring(s))) return 'BioSpring / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => this.isSiRnaAxolabs(s))) return 'Axolabs / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => this.isGcrs(s))) return 'GCRS';
      if (DG.Detector.sampleCategories(col, (s) => this.isMermade12(s))) return 'Mermade 12 / siRNA';
    }
  }
}
