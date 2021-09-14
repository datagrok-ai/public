import * as DG from 'datagrok-api/dg';

class PeptidesPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {
    let regexp = new RegExp("^((\\w+)?-){5,49}\\w+$");
    return Detector.sampleCategories(col, (s) => regexp.test(s)) ? 'alignedSequence' : null;
  }
}
