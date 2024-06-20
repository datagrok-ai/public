/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class AndrewSequencePackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    const validNucleotides = new Set(['A', 'C', 'G', 'T', 'R', 'D', 'N']);
    const toSkip = [' ', '\n'];
    for (const sequence of col.categories) {
      for (const char of sequence) {
        if (toSkip.includes(char)) continue;
        if (!validNucleotides.has(char.toUpperCase()))
          return null;
      }
    }
    return 'dna_nucleotide';
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectENAID(col) {
    // returns semType 'EnaID', if name is of the form [A-Z]{2}[0-9]{6}
    const idTest = /^[A-Z]{2}[0-9]{6}(\.\d+)?$/;
    for (const id of col.categories)
      if (!idTest.test(id)) return null;

    return 'EnaID';
  }
}
