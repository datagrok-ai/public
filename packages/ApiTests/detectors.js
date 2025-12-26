/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class ApiTestsPackageDetectors extends DG.Package {

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectCountry(col) {
    if (col.name === 'countries' && col.type === DG.TYPE.STRING) {
      col.semType = 'country';
      return col.semType;
    }
    return null;
  }
}
