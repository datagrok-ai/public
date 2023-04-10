/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
// eslint-disable-next-line no-unused-vars
class DavitSequencePackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectDetectSequence(col) {
    if (col.type === 'string') {
      let isSeq = true;
      Array.from(col.categories).forEach((seq) =>{
        Array.from(seq.replace(/\s/g, '')).forEach((n) =>{
          if (!'ATGCatgc'.includes(n)) {
            isSeq = false;
            return null;
          }
        });
        if (!isSeq)
          return null;
      });
      if (!isSeq)
        return null;

      col.semType = 'dna_sequences';
      return col.semType;
    }
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectENAID(col) {
    // returns semType 'EnaID', if name is of the form [A-Z]{2}[0-9]{6}
    if (col.type == 'string') {
      let isOfType = true;
      const regexPattern = /^[A-Z]{2}[0-9]{6}$/;
      Array.from(col.categories).forEach((seq) =>{
        if (seq.length !== 8 ) {
          isOfType = false;
          return;
        }
        if (!regexPattern.test(seq)) {
          isOfType = false;
          return;
        }
      });
      if (!isOfType)
        return null;


      col.semType = 'EnaID';
      return col.semType;
    }
    return null;
  }
}


