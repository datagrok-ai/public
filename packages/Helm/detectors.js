/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class should follow the `<PackageName>Detectors` notation, otherwise it won't be loaded.
 */
class HelmPackageDetectors extends DG.Package {

  /** @param s {String} - string to check
   * @returns {boolean} */
  static isHelm(s) {
    return s.startsWith('PEPTIDE1{') || s.startsWith('RNA1{') || s.startsWith('CHEM1{') || s.startsWith('BLOB1{');
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectHelm(col) {
    if (DG.Detector.sampleCategories(col, (s) => HelmPackageDetectors.isHelm(s)))
      return 'HELM';
  }
}
