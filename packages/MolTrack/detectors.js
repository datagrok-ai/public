/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class MoltrackPackageDetectors extends DG.Package {
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectGrokId(col) {
    if (!DG.Detector.sampleCategories(col, (s) => /^DG-\d{6}$/.test(s), 1))
      return null;
    col.semType = 'Grok ID';
    return col.semType;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectGrokBatchId(col) {
    if (!DG.Detector.sampleCategories(col, (s) => /^DGB-\d{6}$/.test(s), 1))
      return null;
    col.semType = 'Grok Batch ID';
    return col.semType;
  }
}
