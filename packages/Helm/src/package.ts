/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GapOriginals, SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';

import {findMonomers, parseHelm} from './utils';
import {HelmCellRenderer} from './cell-renderer';
import {HelmHelper} from './helm-helper';
import {getPropertiesWidget} from './widgets/properties-widget';
import {_getHelmService} from './package-utils';

import {WebEditorMonomer, RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from './constants';

let monomerLib: IMonomerLib | null = null;

export const _package = new DG.Package();

function initHelmPatchDojo(): void {
  // patch window.dojox.gfx.svg.Text.prototype.getTextWidth hangs
  /** get the text width in pixels */
  // @ts-ignore
  window.dojox.gfx.svg.Text.prototype.getTextWidth = function() {
    const rawNode = this.rawNode;
    const oldParent = rawNode.parentNode;
    const _measurementNode = rawNode.cloneNode(true);
    _measurementNode.style.visibility = 'hidden';

    // solution to the "orphan issue" in FF
    let _width = 0;
    const _text = _measurementNode.firstChild.nodeValue;
    oldParent.appendChild(_measurementNode);

    // solution to the "orphan issue" in Opera
    // (nodeValue == "" hangs firefox)
    if (_text != '') {
      let watchdogCounter = 100;
      while (!_width && --watchdogCounter > 0) { // <-- hangs
        //Yang: work around svgweb bug 417 -- http://code.google.com/p/svgweb/issues/detail?id=417
        if (_measurementNode.getBBox)
          _width = parseInt(_measurementNode.getBBox().width);
        else
          _width = 68;
      }
    }
    oldParent.removeChild(_measurementNode);
    return _width;
  };
}

//tags: init
export async function initHelm(): Promise<void> {
  //@ts-ignore
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
      monomerLib = lib;
      rewriteLibraries(); // initHelm()
      monomerLib.onChanged.subscribe((_) => {
        try {
          rewriteLibraries(); // initHelm()

          const monTypeList: string[] = monomerLib!.getPolymerTypes();
          const msgStr: string = 'Monomer lib updated:<br />' + (
            monTypeList.length == 0 ? 'empty' : monTypeList.map((monType) => {
              return `${monType} ${monomerLib!.getMonomerSymbolsByType(monType).length}`;
            }).join('<br />'));

          grok.shell.info(msgStr);
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
  // @ts-ignore
  org.helm.webeditor.Monomers.clear();
  monomerLib!.getPolymerTypes().forEach((polymerType) => {
    const monomerSymbols = monomerLib!.getMonomerSymbolsByType(polymerType);
    monomerSymbols.forEach((monomerSymbol) => {
      let isBroken = false;
      const monomer: Monomer = monomerLib!.getMonomer(polymerType, monomerSymbol)!;
      const webEditorMonomer: WebEditorMonomer = {
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
        // @ts-ignore
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
  return new HelmCellRenderer(); // old
  // return new HelmGridCellRenderer(); // new
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
      `Monomers ${Array.from(missedMonomerSet).map((m) => `'${m}'`).join(', ')} are absent! <br/>` +
      `Please, upload the monomer library! <br/>` +
      `<a href="https://datagrok.ai/help/domains/bio/macromolecules" target="_blank">Learn more</a>`);
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
  // @ts-ignore
  org.helm.webeditor.MolViewer.molscale = 0.8;
  // @ts-ignore
  const app = new scil.helm.App(view, {
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
  //@ts-ignore
  scil.apply(app.sequence.style, s);
  //@ts-ignore
  scil.apply(app.notation.style, s);
  s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + app.toolbarheight) + 'px'};
  //@ts-ignore
  scil.apply(app.properties.parent.style, s);
  app.structureview.resize(sizes.rightwidth, sizes.bottomheight + app.toolbarheight);
  app.mex.resize(sizes.topheight - 80);
  setTimeout(() => {
    if (!!cell && units === undefined)
      app.canvas.helm.setSequence(cell.value, 'HELM');
    else
      app.canvas.helm.setSequence(value, 'HELM');
  }, 200);
  //@ts-ignore
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

function getRS(smiles: string) {
  const newS = smiles.match(/(?<=\[)[^\][]*(?=])/gm);
  const res: { [name: string]: string } = {};
  let el = '';
  let digit;
  if (!!newS) {
    for (let i = 0; i < newS.length; i++) {
      if (newS[i] != null) {
        if (/\d/.test(newS[i])) {
          digit = newS[i][newS[i].length - 1];
          newS[i] = newS[i].replace(/[0-9]/g, '');
          for (let j = 0; j < newS[i].length; j++) {
            if (newS[i][j] != ':')
              el += newS[i][j];
          }
          res['R' + digit] = el;
          el = '';
        }
      }
    }
  }
  return res;
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
