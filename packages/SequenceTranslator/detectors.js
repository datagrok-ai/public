/* eslint-disable max-len */
/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
/// <reference path="../../globals.d.ts" />
class SequenceTranslatorPackageDetectors extends DG.Package {
  //name: autostart
  //meta.role: autostart
  //tags: autostart
  //description: SequenceTranslator bootstrap
  autostart() {
    this.logger.debug('ST: detectors.js: autostart()');

    this.autostartContextMenu();
  }

  autostartContextMenu() {
    const logPrefix = `ST: detectors.js: autostartContextMenu()`;
    this.logger.debug(`${logPrefix}, start`);
    grok.events.onContextMenu.subscribe((event) => {
      try {
        this.logger.debug(`${logPrefix}, onContextMenu, start`);
        if (!event || !event.args || !event.args.item || !event.args.menu)
          return;
        const item = event.args.item;
        const menu = event.args.menu;

        const catchError = (er) => {
          this.logger.error('enumerator error');
          this.logger.error(er?.toString());
        };

        if ((item instanceof DG.GridCell || item.constructor?.name === 'GridCell') && item.tableColumn && item.cell) {
          const packageName = this.name;
          const col = item.tableColumn;
          switch (col.semType) {
          case DG.SEMTYPE.MACROMOLECULE: {
            menu.item('PolyTool-Enumerate', () => { grok.functions.call(`${packageName}:getPtHelmEnumeratorDialog`, {cell: item.cell}).catch(catchError); });

            // Oligo submenu — all oligo-related actions live under one group
            // so they stay discoverable but don't clutter the top context menu.
            const oligoMenu = menu.group('Oligo');

            // HELM → OligoNucleotide: only offered for Macromolecule columns whose units are HELM.
            const units = col.meta?.units || col.tags?.['units'];
            if (units === 'helm') {
              oligoMenu.item('Convert HELM to Oligo', () => {
                grok.functions.call(`${packageName}:convertHelmToOligoNucleotide`, {
                  table: col.dataFrame,
                  helmCol: col,
                }).catch(catchError);
              });
            }

            // Combine sense + antisense → opens the function editor so the user
            // can pick the antisense column (and confirm sense).
            oligoMenu.item('Combine sense+antisense to Oligo...', () => {
              try {
                const fns = DG.Func.find({package: packageName, name: 'combineSenseAntisenseToOligoNucleotide'});
                if (fns.length > 0)
                  fns[0].prepare({table: col.dataFrame, senseCol: col}).edit();
              } catch (er) { catchError(er); }
            });

            oligoMenu.endGroup();
            break;
          }
          case DG.SEMTYPE.MOLECULE: {
            menu.item('PolyTool-Enumerate', () => { grok.functions.call(`${packageName}:getPtChemEnumeratorDialog`, {cell: item.cell}).catch(catchError); });
            break;
          }
          }
        }
        this.logger.debug(`${logPrefix}, onContextMenu, end`);
      } catch (error) {
        this.logger.error(`${logPrefix}, error`);
        this.logger.error(error?.toString());
        if (!window.$sequenceTranslator) window.$sequenceTranslator = {};
        window.$sequenceTranslator.contextMenuError = error;
      }
    });
  }

  // //meta.role: semTypeDetector
  // //input: column col
  // //output: string semType
  // detectHarmonizedSequence(col) {
  //   // To prevent "Error initializing SequenceTranslator 1.4.4.X-ba98e6f7: ReferenceError: SequenceTranslatorPackageDetectors is not defined"
  //   return null;
  // }


  //name: refineNotationProviderForHarmonizedSequence
  //meta.role: notationRefiner
  //tags: notationRefiner
  //input: column col
  //input: object stats
  //input: string separator = null { optional: true }
  //output: bool result
  async refineNotationProviderForHarmonizedSequence(col, stats, separator) {
    const isCyclized = Object.keys(stats.freq).some((om) => om.match(/^.+\(\d+\)$/));
    const isDimerized = Object.keys(stats.freq).some((om) => om.match(/^\(#\d\).+$/));

    if (!isCyclized && !isDimerized) return false;

    await grok.functions.call('SequenceTranslator:applyNotationProviderForCyclized', {col, separator});
    return true;
  }
}
