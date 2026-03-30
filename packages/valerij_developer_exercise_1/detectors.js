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
        .replace(/^FASTA:\s*/, '');

      if (!/^[ATGC]+$/.test(s))
        return null;
    }

    return 'dna_nucleotide';
  }
}