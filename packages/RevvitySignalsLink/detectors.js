/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
const REVVITY_SIGNALS_ID_SEM_TYPE = 'RevvitySignalsId';
class RevvitysignalslinkPackageDetectors extends DG.Package {

    static entityTypes = [
        'experiment',
        'request',
        'task',
        'journal',
        'chemicalDrawing',
        'grid',
        'materialsTable',
        'acronym',
        'text',
        'pdf',
        'viewonly',
        'presentation',
        'excel',
        'imageResource',
        'uploadedResource',
        'spotfiredxp',
        'sample',
        'assetType',
        'asset',
        'batch',
        'paraexp',
        'parasubexp',
        'sampleContainer',
        'worksheet',
    ];

    static re = new RegExp(`^(${RevvitysignalslinkPackageDetectors.entityTypes.join('|')}):[0-9a-fA-F]+(-[0-9a-fA-F]+)*$`);

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectRevvitySignalsId(col) {
      const REVVITY_SIGNALS_ID_SEM_TYPE = 'RevvitySignalsId';
      if (col.type === DG.TYPE.STRING && DG.Detector.sampleCategories(col, (s) => {
        return RevvitysignalslinkPackageDetectors.re.test(s);
      }, 1)) {
        col.semType = REVVITY_SIGNALS_ID_SEM_TYPE;
        return col.semType;
      }
      return null;
    }
}
