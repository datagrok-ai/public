/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */

const COL_NAMES = {
  CHEMISTRY: 'Chemistry',
  NUMBER: 'Number',
  TYPE: 'Type',
  CHEMISTRY_NAME: 'Chemistry Name',
  INTERNAL_COMPOUND_ID: 'Internal compound ID',
  IDP: 'IDP',
  SEQUENCE: 'Sequence',
  COMPOUND_NAME: 'Compound Name',
  COMPOUND_COMMENTS: 'Compound Comments',
  SALT: 'Salt',
  EQUIVALENTS: 'Equivalents',
  PURITY: 'Purity',
  COMPOUND_MOL_WEIGHT: 'Cpd MW',
  SALT_MOL_WEIGHT: 'Salt MW',
  SALT_MASS: 'Salt mass',
  BATCH_MOL_WEIGHT: 'Batch MW',
  SOURCE: 'Source',
  ICD: 'ICD',
  OWNER: 'Owner',
};

class SequenceTranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotide(col) {
    return null;
  }

  isMatchView(view) {
    if (view.type !== DG.VIEW_TYPE.TABLE_VIEW) return false;
    const colNameList = view.dataFrame.columns.names();
    return [COL_NAMES.TYPE].every((requiredColName) => colNameList.includes(requiredColName));
  }

  //name: autostart
  //tags: autostart
  //description: SequenceTranslator bootstrap
  autostart() {
    console.debug('ST: detectors.js: autostart()');

    const eventsOnViewAdded = (view) => {
      console.debug('ST: detectors.js: eventsOnViewAdded() handler start');
      if (this.isMatchView(view)) {
        // Do not await view engaging
        grok.functions.call('SequenceTranslator:engageViewForOligoSdFile', {view: view})
          .catch((err) => {
            const errMsg = err.hasOwnProperty('message') ? err.message : err.toString();
            grok.shell.error('Engaging view for SequenceTranslator OligoSdFile error: ' + errMsg);
          });
      }
    };

    const viewList = [...grok.shell.views];
    for (const view of viewList) {
      console.debug('ST: detectors.js: autostart() view found \'' + view.name + '\'');
      eventsOnViewAdded(view);
    }
    grok.events.onViewAdded.subscribe(eventsOnViewAdded);
  }
}
