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
  static likelySvg(s) {
    if (s == null || s.length < 50)
      return false;
    const t = s.trimStart();
    return t.startsWith('<svg') || (t.startsWith('<?xml') && t.substring(0, 200).includes('<svg'));
  }

  static likelyImageUrl(s) {
    if (s == null)
      return false;
    s = s.toLowerCase();
    return (s.startsWith('http') || s.startsWith('system:')) && ['png', 'jpg', 'jpeg'].some((x) => s.endsWith(x));
  }

  /* Checks the file signature of the formats browsers decode: PNG, JPEG, GIF, BMP, WebP, ICO and AVIF.
   * @param b {Uint8Array}
   * @returns {boolean} */
  static likelyImageBytes(b) {
    const starts = (...signature) => signature.every((x, i) => b[i] === x);
    const ascii = (from, s) => [...s].every((c, i) => b[from + i] === c.charCodeAt(0));
    return b.length > 12 && (starts(0x89, 0x50, 0x4E, 0x47) || starts(0xFF, 0xD8, 0xFF) || ascii(0, 'GIF8') ||
      ascii(0, 'BM') || (ascii(0, 'RIFF') && ascii(8, 'WEBP')) || starts(0, 0, 1, 0) || starts(0, 0, 2, 0) ||
      (ascii(4, 'ftyp') && /avi[fs]/.test(String.fromCharCode(...b.subarray(8, 32)))));
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectSvg(col) {
    if (col.type === DG.COLUMN_TYPE.STRING &&
        DG.Detector.sampleCategories(col, (s) => PowerGridPackageDetectors.likelySvg(s), 1))
      return 'SVG';
    return null;
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
    if (col.type !== DG.COLUMN_TYPE.BYTE_ARRAY)
      return null;
    let checked = 0;
    for (let i = 0; i < Math.min(col.length, 10000) && checked < 100; i++) {
      const bytes = col.get(i);
      if (bytes == null || bytes.length === 0)
        continue;
      if (PowerGridPackageDetectors.likelyImageBytes(bytes))
        return 'BinaryImage';
      checked++;
    }
    // an empty column stays an image column, to be filled by double-clicking its cells
    return checked === 0 ? 'BinaryImage' : null;
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
