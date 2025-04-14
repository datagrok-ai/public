/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name consists of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class PowerPackPackageDetectors extends DG.Package {
  /* @param s {String} - string to check
   * @returns {boolean} */
  static likelyText(s) {
    if (s == null)
      return false;
    return false;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //meta.skipTest: GROK-17759
  detectText(col) {
    if (col.type === DG.COLUMN_TYPE.STRING && col.name === 'Questions')
      col.semType = DG.SEMTYPE.TEXT;

    // if (DG.Detector.sampleCategories(col, (s) => PowerPackGridPackageDetectors.likelyText(s), 1))
    //   return DG.SEMTYPE.TEXT;
  }
}
