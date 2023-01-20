/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class BiostructureViewerPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPdb(col) {
    let res = null;
    if (DG.Detector.sampleCategories(col,
      (s) => s.includes('COMPND') && s.includes('ATOM') && s.includes('END'), 1)
    ) {
      res = 'Molecule3D';
    }
    return res;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPdbId(col) {
    let res = null;
    if (col.type === DG.TYPE.STRING &&
      col.name.toLowerCase().includes('pdb') &&
      DG.Detector.sampleCategories(col, (s) => s.length === 4)
    ) {
      res = 'PDB_ID';
    }
    return res;
  }
}
