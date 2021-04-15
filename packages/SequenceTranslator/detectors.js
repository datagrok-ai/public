import {isClassicCode, isAsoGapmerBioSpringCode, isSiRnaAxolabsCode, isAsoGapmerGCRSCode, isSiRnaBioSpringCode, isSiRnaGCRSCode} from "./src/package.ts";

class SequencetranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (DG.Detector.sampleCategories(col, (s) => isClassicCode(s))) return 'nucleotides';
      if (DG.Detector.sampleCategories(col, (s) => isAsoGapmerBioSpringCode(s))) return 'BioSpring / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => isAsoGapmerGCRSCode(s))) return 'GCRS / Gapmers';
      if (DG.Detector.sampleCategories(col, (s) => isSiRnaBioSpringCode(s))) return 'BioSpring / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => isSiRnaAxolabsCode(s))) return 'Axolabs / siRNA';
      if (DG.Detector.sampleCategories(col, (s) => isSiRnaGCRSCode(s))) return 'GCRS / siRNA';
    }
  }
}