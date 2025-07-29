/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
const FIT_SEM_TYPE = 'fit';
const TAG_FIT_CHART_FORMAT = '.fitChartFormat';
const TAG_FIT_CHART_FORMAT_3DX = '3dx';
const RAW_PNG_SEM_TYPE = 'rawPng';
const likelyPNGColNames = ['image', 'png', 'thumbnail', 'img'];
/// <reference path="../../globals.d.ts" />
class CurvesPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectXMLCurveChart(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.startsWith('<chart>') && s.endsWith('</chart>'), 1)) {
      col.setTag(TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX);
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectFit(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      return s.includes('series') && s.includes('points');
    }, 1)) {
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectRawPng(col) {
    const colName = col.name.toLowerCase();
    if (col.type === DG.COLUMN_TYPE.STRING &&
      likelyPNGColNames.some((name) => colName.includes(name)) && DG.Detector.sampleCategories(col, (s) => {
      return s && s.length > 1000 && s.startsWith('iVBORw0KGgo');
    }, 1))
      return RAW_PNG_SEM_TYPE;
    return null;
  }
}
