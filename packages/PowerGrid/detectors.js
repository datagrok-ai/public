/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class PowerGridPackageDetectors extends DG.Package {
  /* @param s {String} - string to check
   * @returns {boolean} */
  static likelyImageUrl(s) {
    if (s == null)
      return false;
    s = s.toLowerCase();
    return (s.startsWith('http') || s.startsWith('system:')) && ['png', 'jpg', 'jpeg'].some((x) => s.endsWith(x));
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectImageColumn(col) {
    if (DG.Detector.sampleCategories(col, (s) => PowerGridPackageDetectors.likelyImageUrl(s), 1))
      return 'ImageUrl';
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectByteArrayColumn(col) {
    if (col.type === DG.COLUMN_TYPE.BYTE_ARRAY)
      return 'BinaryImage';
  }
}
