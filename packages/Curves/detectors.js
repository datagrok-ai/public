/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
const FIT_SEM_TYPE = 'fit';
class CurvesPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectXMLCurveChart(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.startsWith('<chart>') && s.endsWith('</chart>'), 1)) {
      col.setTag('.fitChartFormat', '3dx');
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }
}
