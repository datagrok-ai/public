/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class BiostructureviewerPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPDB(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.includes('COMPND') && s.includes('ATOM') && s.includes('END'), 1))
      return 'xray';

    return null;
  }
}
