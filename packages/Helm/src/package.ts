/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import '@datagrok-libraries/bio/src/types/dojo';
// import * as dojo from 'DOJO';
import '@datagrok-libraries/bio/src/types/helm';
import * as scil from 'scil';
import * as org from 'org';
import '@datagrok-libraries/bio/src/types/jsdraw2';
import * as JSDraw2 from 'JSDraw2';

import $ from 'cash-dom';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GapOriginals, SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {getRS} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';

import {findMonomers, parseHelm} from './utils';
import {HelmCellRenderer} from './cell-renderer';
import {HelmHelper} from './helm-helper';
import {getPropertiesWidget} from './widgets/properties-widget';
import {HelmGridCellRenderer, HelmGridCellRendererBack} from './utils/helm-grid-cell-renderer';
import {_getHelmService, initHelmPatchDojo, initHelmPatchPistoia} from './package-utils';

let monomerLib: IMonomerLib | null = null;

export const _package = new DG.Package();

//tags: init
export async function initHelm(): Promise<void> {
  const logPrefix: string = 'Helm: initHelm()';
  _package.logger.debug(`${logPrefix}, start`);
  org.helm.webeditor.kCaseSensitive = true; // GROK-13880

  return Promise.all([
    new Promise((resolve, reject) => {
      // @ts-ignore
      dojo.ready(function() { resolve(null); });
    }).then(() => {
      initHelmPatchDojo();
    }),
    grok.functions.call('Bio:getBioLib'),
  ])
    .then(([_, lib]: [unknown, IMonomerLib]) => {
      _package.logger.debug(`${logPrefix}, then(), lib loaded`);
      monomerLib = lib;
      rewriteLibraries(); // initHelm()
      initHelmPatchPistoia(monomerLib, _package.logger);
      monomerLib.onChanged.subscribe((_) => {
        try {
          const logPrefixInt = `${logPrefix} monomerLib.onChanged()`;
          _package.logger.debug(`${logPrefixInt}, start, org,helm.webeditor.Monomers updating`);
          rewriteLibraries(); // initHelm()

          const libSummary: string = monomerLib!.getSummary()
            .replaceAll('\n', '<br />\n')
            .replaceAll(',', ' ');
          const libMsg: string = `Monomer lib updated:<br /> ${libSummary}'`;
          grok.shell.info(libMsg);

          _package.logger.debug(`${logPrefixInt}, end, org.helm.webeditor.Monomers completed`);
        } catch (err: any) {
          const errMsg = errorToConsole(err);
          console.error('Helm: initHelm monomerLib.onChanged() error:\n' + errMsg);
          // throw err; // Prevent disabling event handler
        }
      });
    })
    .catch((err: any) => {
      const errMsg: string = err instanceof Error ? err.message : !!err ? err.toString() : 'Exception \'undefined\'';
      grok.shell.error(`Package \'Helm\' init initHelm() error: ${errMsg}`);
      const errRes = new Error(errMsg);
      errRes.stack = err.stack;
      throw errRes;
    });
}

export function getMonomerLib(): IMonomerLib {
  return monomerLib!;
}

/** Fills org.helm.webeditor.Monomers dictionary for WebEditor */
function rewriteLibraries() {
  org.helm.webeditor.Monomers.clear();
  monomerLib!.getPolymerTypes().forEach((polymerType) => {
    const monomerSymbols = monomerLib!.getMonomerSymbolsByType(polymerType);
    monomerSymbols.forEach((monomerSymbol) => {
      let isBroken = false;
      const monomer: Monomer = monomerLib!.getMonomer(polymerType, monomerSymbol)!;
      const webEditorMonomer: org.helm.WebEditorMonomer = {
        id: monomerSymbol,
        m: monomer.molfile,
        n: monomer.name,
        na: monomer.naturalAnalog,
        rs: monomer.rgroups.length,
        type: monomer.polymerType,
        mt: monomer.monomerType,
        at: {}
      };

      if (monomer.rgroups.length > 0) {
        webEditorMonomer.rs = monomer.rgroups.length;
        const at: { [prop: string]: any } = {};
        monomer.rgroups.forEach((it) => {
          at[it[RGROUP_LABEL]] = it[RGROUP_CAP_GROUP_NAME];
        });
        webEditorMonomer.at = at;
      } else if (monomer[SMILES] != null) {
        webEditorMonomer.rs = Object.keys(getRS(monomer[SMILES].toString())).length;
        webEditorMonomer.at = getRS(monomer[SMILES].toString());
      } else
        isBroken = true;

      if (!isBroken) {
        org.helm.webeditor.Monomers.addOneMonomer(webEditorMonomer);
      }
    });
  });

  // Obsolete
  const grid: DG.Grid = grok.shell.tv?.grid;
  if (grid) grid.invalidate();
}

//name: getHelmService
//output: object result
export function getHelmService(): HelmServiceBase {
  return _getHelmService();
}

//name: helmCellRenderer
//tags: cellRenderer
//meta.cellType: helm
//meta.columnTags: quality=Macromolecule, units=helm
//output: grid_cell_renderer result
export function helmCellRenderer(): HelmCellRenderer {
  // return new HelmCellRenderer(); // old
  return new HelmGridCellRenderer(); // new
}

function checkMonomersAndOpenWebEditor(cell: DG.Cell, value?: string, units?: string) {
  const cellValue: string = !!cell && units === undefined ? cell.value : value;
  const monomerList: string[] = parseHelm(cellValue);
  const missedMonomerSet = findMonomers(monomerList);
  if (missedMonomerSet.size === 0)
    webEditor(cell, value, units);
  else if (missedMonomerSet.size === 1 && missedMonomerSet.has(GapOriginals[NOTATION.HELM]))
    grok.shell.warning(`WebEditor doesn't support Helm with gaps '${GapOriginals[NOTATION.HELM]}'.`);
  else {
    grok.shell.warning(
      `Monomers ${Array.from(missedMonomerSet).map((m) => `'${m}'`).join(', ')} are not found. <br/>` +
      `Please specify the monomer library to use. <br/>` +
      `<a href="https://datagrok.ai/help/datagrok/solutions/domains/bio/#manage-monomer-libraries" target="_blank">Learn more</a>`);
  }
}

//tags: cellEditor
//description: Macromolecule
//input: grid_cell cell
//meta.columnTags: quality=Macromolecule, units=helm
export function editMoleculeCell(cell: DG.GridCell): void {
  checkMonomersAndOpenWebEditor(cell.cell, undefined, undefined);
}

//name: Open Helm Web Editor
//description: Adds editor
//meta.action: Open Helm Web Editor
//input: string mol { semType: Macromolecule }
export function openEditor(mol: string): void {
  const df = grok.shell.tv.grid.dataFrame;
  const col = df.columns.bySemType('Macromolecule')! as DG.Column<string>;
  const colSh = SeqHandler.forColumn(col);
  const colUnits = col.getTag(DG.TAGS.UNITS);
  if (colUnits === NOTATION.HELM)
    checkMonomersAndOpenWebEditor(df.currentCell, undefined, undefined);
  const convert = colSh.getConverter(NOTATION.HELM);
  const helmMol = convert(mol);
  checkMonomersAndOpenWebEditor(df.currentCell, helmMol, col.getTag(DG.TAGS.UNITS));
}

//name: Properties
//tags: panel, bio, helm, widgets
//input: semantic_value sequence {semType: Macromolecule}
//output: widget result
export function propertiesWidget(sequence: DG.SemanticValue): DG.Widget {
  return getPropertiesWidget(sequence);
}

function webEditor(cell: DG.Cell, value?: string, units?: string) {
  const view = ui.div();
  // const df = grok.shell.tv.grid.dataFrame;
  // const col = df.columns.bySemType('Macromolecule')!;
  const col = cell.column as DG.Column<string>;
  const sh = SeqHandler.forColumn(col);
  const rowIdx = cell.rowIndex;
  org.helm.webeditor.MolViewer.molscale = 0.8;
  const app = new org.helm.webeditor.App(view, {
    showabout: false,
    mexfontsize: '90%',
    mexrnapinontab: true,
    topmargin: 20,
    mexmonomerstab: true,
    sequenceviewonly: false,
    mexfavoritefirst: true,
    mexfilter: true
  });
  const sizes = app.calculateSizes();
  app.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
  let s = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
  scil.apply(app.sequence.style, s);
  scil.apply(app.notation.style, s);
  s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + app.toolbarheight) + 'px'};
  scil.apply(app.properties.parent.style, s);
  app.structureview.resize(sizes.rightwidth, sizes.bottomheight + app.toolbarheight);
  app.mex.resize(sizes.topheight - 80);
  setTimeout(() => {
    if (!!cell && units === undefined)
      app.canvas.helm.setSequence(cell.value, 'HELM');
    else
      app.canvas.helm.setSequence(value!, 'HELM');
  }, 200);
  ui.dialog({showHeader: false, showFooter: true})
    .add(view)
    .onOK(() => {
      const helmValue: string = app.canvas.getHelm(true).replace(/<\/span>/g, '')
        .replace(/<span style='background:#bbf;'>/g, '');
      if (!!cell) {
        if (units === undefined)
          cell.value = helmValue;
        else {
          const convertedRes = sh.convertHelmToFastaSeparator(helmValue, units!, sh.separator);
          cell.value = convertedRes;
        }
      }
    }).show({modal: true, fullScreen: true});
}

//name: getMolfiles
//input: column col {semType: Macromolecule}
//output: column res
export function getMolfiles(col: DG.Column): DG.Column {
  const res = DG.Column.string('mols', col.length);

  const host = ui.div([], {style: {width: '0', height: '0'}});
  document.documentElement.appendChild(host);
  try {
    const editor = new JSDraw2.Editor(host, {viewonly: true});
    res.init((i) => {
      editor.setHelm(col.get(i));
      const mol = editor.getMolfile();
      return mol;
    });
    return res;
  } finally {
    $(host).empty();
    host.remove();
  }
}

// -- Utils --

//name: getHelmHelper
//output: object result
export async function getHelmHelper(): Promise<IHelmHelper> {
  return HelmHelper.getInstance();
}

import {testEvent} from '@datagrok-libraries/utils/src/test';
import {CellRendererBackAsyncBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from './constants';

//name: measureCellRenderer
export async function measureCellRenderer(): Promise<void> {
  const grid = grok.shell.tv.grid;
  const gridCol = grid.columns.byName('sequence')!;
  const back = gridCol.temp['rendererBack'] as HelmGridCellRendererBack;

  let etSum: number = 0;
  let etCount: number = 0;
  for (let i = 0; i < 20; ++i) {
    const t1 = window.performance.now();
    let t2: number;
    if (!back.cacheEnabled) {
      await testEvent(back.onRendered, () => {
        t2 = window.performance.now();
        _package.logger.info(`measureCellRenderer() cache disabled , ET: ${t2 - t1} ms`);
      }, () => {
        back.invalidate(); // grid.invalidate();
      }, 5000);
    } else {
      await testEvent(grid.onAfterDrawContent, () => {
        t2 = window.performance.now();
        _package.logger.info(`measureCellRenderer() cache enabled, ET: ${t2! - t1} ms`);
      }, () => {
        grid.invalidate();
      }, 5000);
    }

    etSum += (t2! - t1);
    etCount++;
  }
  _package.logger.info(`measureCellRenderer(), avg ET: ${etSum / etCount} ms`);
}
