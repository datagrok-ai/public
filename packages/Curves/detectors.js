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
const TAG_CURVE_FORMAT = '.%curve-format';
const likelyPNGColNames = ['image', 'png', 'thumbnail', 'img', 'glyph', 'icon', 'picture', 'photo'];
/// <reference path="../../globals.d.ts" />
class CurvesPackageDetectors extends DG.Package {
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectXMLCurveChart(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.startsWith('<chart>') && s.endsWith('</chart>'), 1)) {
      col.setTag(TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX);
      col.setTag(TAG_CURVE_FORMAT, '3dx');
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
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

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectPzfxCurveChart(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      return s.includes('<Table') && s.includes('TableType="XY"') && s.includes('<XColumn');
    }, 1)) {
      col.setTag(TAG_CURVE_FORMAT, 'pzfx');
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectCompactDoseResponse(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      try {
        if (!s.includes('"p"') || !s.includes('"poi"'))
          return false;
        const obj = JSON.parse(s);
        return obj.p !== undefined && Array.isArray(obj.p) && obj.poi !== undefined;
      } catch (_e) {
        return false;
      }
    }, 1)) {
      col.setTag(TAG_CURVE_FORMAT, 'compact-dr');
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }
}
