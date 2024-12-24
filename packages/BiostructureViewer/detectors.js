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
  //meta.skipTest: #2596, Fix for test data in the utils library
  detectPdb(col) {
    if (DG.Detector.sampleCategories(col,
      // (s) => s.includes('COMPND') && s.includes('ATOM') && s.includes('END'), 1)
      (s) => s.match(/^COMPND/m) && s.match(/^END/m) &&
        (s.match(/^ATOM/m) || s.match(/^HETATM/m)), 1
    )) {
      col.meta.units = 'pdb';
      return 'Molecule3D';
    } else if (DG.Detector.sampleCategories(col,
      (s) => s.match(/^MODEL/m) && s.match(/^ENDMDL/m) &&
        (s.match(/^ATOM/m) || s.match(/^HETATM/m)),
      1)
    ) {
      col.meta.units = 'pdbqt';
      return 'Molecule3D';
    }

    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPdbId(col) {
    let res = null;
    if (col.type === DG.TYPE.STRING &&
      col.name.toLowerCase().includes('pdb') &&
      DG.Detector.sampleCategories(col, (s) => s.length === 4)
    )
      res = 'PDB_ID';
    return res;
  }

  //name: autostart
  //tags: autostart
  //description: BiostructureViewer bootstrap
  autostart() {
    this.logger.debug('BsV: detectors.js: autostart()');

    this.autostartContextMenu();
  }

  autostartContextMenu() {
    const logPrefix = 'BsV: detectors.js: autostartContextMenu()';
    this.logger.debug(`${logPrefix}, start`);

    const catchError = (err) => {
      this.logger.error(`${logPrefix} error`);
      this.logger.error(err?.message ?? String(err));
    };

    grok.events.onContextMenu.subscribe((event) => {
      this.logger.debug(`${logPrefix}, onContextMenu, start`);
      try {
        if (!event || !event.args || !event.args.item || !event.args.menu)
          return;
        const item = event.args.item;
        const menu = event.args.menu;
        const packageName = this.name;
        if (item.tableColumn && item.cell) {
          // if its grid cell
          if ((item instanceof DG.GridCell || item.constructor?.name === 'GridCell') &&
          item.tableColumn.semType === DG.SEMTYPE.MOLECULE3D) {
            menu.item('Copy', () => {
              grok.functions.call(`${packageName}:copyRawBiostructureValue`, {gridCell: item}).catch(catchError);
            });
            menu.item('Download', () => {
              grok.functions.call(`${packageName}:downloadRawBiostructureValue`, {gridCell: item}).catch(catchError);
            });
            const group = menu.group('Show');
            group.item('Biostructure', () => {
              grok.functions.call(`${packageName}:showBiostructureViewerMenuItem`, {gridCell: item}).catch(catchError);
            }, null,
            {description: 'Show with Biostructure (mol*) viewer'});
            group.item('NGL', () => {
              grok.functions.call(`${packageName}:showNglViewerMenuItem`, {gridCell: item}).catch(catchError);
            }, null,
            {description: 'Show with NGL viewer'});
          }
        } else if ((item instanceof DG.FileInfo || item.constructor.name == 'FileInfo') &&
          item.extension.toLowerCase() === 'pdb') {
          menu.item('Open table residues', () => {
            grok.functions.call(`${packageName}:openTableResiduesMenuItem`, {fi: item}).catch(catchError);
          });
        }
      } catch (error) {
        this.logger.error(`${logPrefix}, error`);
        if (!window.$biostructureViewer) window.$biostructureViewer = {};
        window.$biostructureViewer.contextMenuError = error;
      }
      this.logger.debug(`${logPrefix}, onContextMenu, end`);
    });
  }
}
