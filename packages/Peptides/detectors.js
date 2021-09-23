class PeptidesPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {
    let regexp = new RegExp("^((\\w+)?-){5,49}\\w+$");
    return DG.Detector.sampleCategories(col, (s) => regexp.test(s)) ? 'alignedSequence' : null;
  }

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectAminoAcids(col) {
        for (let i = 0; i < col.categories.length; i++) {
            if (col.categories[i]) {
                if (!['_','NH2','COOH','C', 'U', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q']
                    .includes(col.categories[i])) {
                    return null;
                }
            }
        }
        return 'aminoAcids';
    }

}
