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
    const validNucleotides = new Set(['A', 'C', 'G', 'T']);
    for (const sequence of col.categories) {
      for (const char of sequence) {
        if (char === ' ') continue;
        if (!validNucleotides.has(char.toUpperCase()))
          return null;
      }
    }
    return 'dna_nucleotide';
  }
}
