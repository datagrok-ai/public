/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class GuillermoSequencePackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //meta.testData: test-sequences.csv
  //meta.testDataColumnName: dna_nucleotide
  detectNucleotides(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    // Unprefixed values (e.g. the bio demo files) must be long enough not to match ordinary words like "cat".
    const isNucleotides = (c) => /^fasta:\s*[ATGC]+$/.test(c) ||
      (/^[ACGT\s]+$/i.test(c) && c.replace(/\s/g, '').length >= 20);
    return col.categories.every(isNucleotides) ? 'dna_nucleotide' : null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //meta.testData: test-ena-ids.csv
  //meta.testDataColumnName: EnaID
  detectEnaId(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    // The optional version suffix covers the demo files, which store IDs like MT079845.1.
    return col.categories.every((c) => /^[A-Z]{2}\d{6}(\.\d+)?$/.test(c)) ? 'EnaID' : null;
  }
}
