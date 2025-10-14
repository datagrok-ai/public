/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */

TEST_SEM_TYPE = 'test';
class UsageAnalysisPackageDetectors extends DG.Package {
    //tags: semTypeDetector
    //input: column col
    //output: string semType 
    //meta.skipTest: test
    detectTest(col) {
      if ((col.name.toLowerCase().includes('test') | col.name.toLowerCase().includes('name')) && DG.Detector.sampleCategories(col, (s) => {
        return new RegExp(".*:.*:.*").test(s);
      }, 5)) {
        col.semType = TEST_SEM_TYPE;
        return col.semType;
      }
      return null;
    }
  }