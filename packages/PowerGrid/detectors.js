/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
const RAW_PNG_SEM_TYPE = 'rawPng';
class PowerGridPackageDetectors extends DG.Package {
  /* @param s {String} - string to check
   * @returns {boolean} */
  static likelyImageUrl(s) {
    if (s == null)
      return false;
    s = s.toLowerCase();
    return (s.startsWith('http') || s.startsWith('system:')) && ['png', 'jpg', 'jpeg'].some((x) => s.endsWith(x));
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectImageColumn(col) {
    if (DG.Detector.sampleCategories(col, (s) => PowerGridPackageDetectors.likelyImageUrl(s), 1))
      return 'ImageUrl';
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectByteArrayColumn(col) {
    if (col.type === DG.COLUMN_TYPE.BYTE_ARRAY)
      return 'BinaryImage';
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectRawPng(col) {
    if (col.type === DG.COLUMN_TYPE.STRING && DG.Detector.sampleCategories(col, (s) => {
      return s && s.length > 1000 && s.startsWith('iVBORw0KGgo');
    }, 1))
      return RAW_PNG_SEM_TYPE;
    return null;
  }
}
