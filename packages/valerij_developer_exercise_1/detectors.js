/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 */
class Valerij_developer_exercise_1PackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type !== 'string')
      return null;

    const categories = col.categories;
    if (!categories || categories.length === 0)
      return null;

    for (const value of categories) {
      if (value == null)
        continue;

      const s = value
        .toString()
        .trim()
        .toUpperCase()
        .replace(/^FASTA:\s*/, '')
        .replace(/\s+/g, '');

      if (s.length === 0)
        return null;

      if (!/^[ACGTNRYSWKMDHBV-]+$/.test(s))
        return null;
    }

    return 'dna_nucleotide';
  }
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectENAID(col) {
    if (col.type !== 'string')
      return null;

    const categories = col.categories;
    if (!categories || categories.length === 0)
      return null;

    for (const value of categories) {
      if (value == null)
        continue;

      if (!/^[A-Z]{2}[0-9]{6}$/.test(value.toString().trim()))
        return null;
    }

    return 'EnaID';
  }
}
