/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */

const GROK_ID_SEM_TYPE = 'Grok ID';
const GROK_BATCH_ID_SEM_TYPE = 'Grok Batch ID';

class MoltrackPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGrokId(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      return s.includes('DG-');
    }, 1)) {
      col.semType = GROK_ID_SEM_TYPE;
      return col.semType;
    }
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGrokBatchId(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      return s.includes('DGB-');
    }, 1)) {
      col.semType = GROK_BATCH_ID_SEM_TYPE;
      return col.semType;
    }
    return null;
  }
}
