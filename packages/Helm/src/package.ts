/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {WebEditorMonomer, RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from './constants';
import {HelmWebEditor} from './helm-web-editor';
import {HelmCellRenderer} from './cell-renderer';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {NotationConverter} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {findMonomers} from './utils';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

export const _package = new DG.Package();
let monomerLib: IMonomerLib | null = null;

//tags: init
export async function initHelm(): Promise<void> {
  return Promise.all([new Promise((resolve, reject) => {
    // @ts-ignore
    dojo.ready(function() { resolve(null); });
  }), await grok.functions.call('Bio:getBioLib')]).then(([_, lib]: [void, IMonomerLib]) => {
    monomerLib = lib;
    rewriteLibraries(); // initHelm()
    monomerLib.onChanged.subscribe((_) => {
      try {
        rewriteLibraries(); // initHelm()

        const monTypeList: string[] = monomerLib.getTypes();
        const msgStr: string = 'Monomer lib updated:<br />' + (
          monTypeList.length == 0 ? 'empty' : monTypeList.map((monType) => {
            return `${monType} ${monomerLib.getMonomerNamesByType(monType).length}`;
          }).join('<br />'));

        grok.shell.info(msgStr);
      } catch (err: any) {
        const errMsg = errorToConsole(err);
        console.error('Helm: initHelm monomerLib.onChanged() error:\n' + errMsg);
        // throw err; // Prevent disabling event handler
      }
    });
  });
}

function rewriteLibraries() {
  org.helm.webeditor.Monomers.clear();
  monomerLib.getTypes().forEach((type) => {
    const mms = monomerLib.getMonomerNamesByType(type);
    mms.forEach((symbol) => {
      let isBroken = false;
      const monomer: Monomer = monomerLib.getMonomer(type, symbol);
      const webEditorMonomer: WebEditorMonomer = {
        id: symbol,
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
        const at = {};
        monomer.rgroups.forEach(it => {
          at[it[RGROUP_LABEL]] = it[RGROUP_CAP_GROUP_NAME];
        });
        webEditorMonomer.at = at;
      } else if (monomer.data[SMILES] != null) {
        webEditorMonomer.rs = Object.keys(getRS(monomer.data[SMILES].toString())).length;
        webEditorMonomer.at = getRS(monomer.data[SMILES].toString());
      } else {
        isBroken = true;
      }

      if (!isBroken)
        org.helm.webeditor.Monomers.addOneMonomer(webEditorMonomer);
    });
  });

  // Obsolete
  let grid: DG.Grid = grok.shell.tv.grid;
  if (grid) grid.invalidate();
}

//name: helmCellRenderer
//tags: cellRenderer
//meta.cellType: helm
//meta.columnTags: units=helm
//output: grid_cell_renderer result
export function helmCellRenderer(): HelmCellRenderer {
  return new HelmCellRenderer();
}

function checkMonomersAndOpenWebEditor(cell?: DG.Cell, value?: string, units?: string) {
  const cellValue = typeof units === 'undefined' ? cell.value : value;
  const monomers = findMonomers(cellValue);
  if (monomers.size == 0) 
    webEditor(cell, value, units);
  else
    grok.shell.warning(`Monomers ${Array.from(monomers).join(', ')} are absent! <br/>` +
    `Please, upload the monomer library! <br/>` +
    `<a href="https://datagrok.ai/help/domains/bio/macromolecules" target="_blank">Learn more</a>`);
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
  let df = grok.shell.tv.grid.dataFrame;
  let col = df.columns.bySemType('Macromolecule');
  let converter = new NotationConverter(col);
  const resStr = converter.convertStringToHelm(mol, '/');
  checkMonomersAndOpenWebEditor(df.currentCell, resStr, col.getTag(DG.TAGS.UNITS));
}

//name: Properties
//tags: panel, widgets
//input: string helmString {semType: Macromolecule}
//output: widget result
export async function propertiesPanel(helmString: string) {
  let grid = grok.shell.tv.grid;
  let parent = grid.root.parentElement;
  const host = ui.div([]);
  parent.appendChild(host);
  let editor = new JSDraw2.Editor(host, {viewonly: true});
  host.style.width = '0px';
  host.style.height = '0px';
  editor.setHelm(helmString);
  let formula = editor.getFormula(true);
  let molWeight = Math.round(editor.getMolWeight() * 100) / 100;
  let coef = Math.round(editor.getExtinctionCoefficient(true) * 100) / 100;
  parent.lastChild.remove();
  return new DG.Widget(
    ui.tableFromMap({
      'formula': formula.replace(/<sub>/g, '').replace(/<\/sub>/g, ''),
      'molecular weight': molWeight,
      'extinction coefficient': coef,
    })
  );
}

function webEditor(cell?: DG.Cell, value?: string, units?: string) {
  let view = ui.div();
  let df = grok.shell.tv.grid.dataFrame;
  let converter = new NotationConverter(df.columns.bySemType('Macromolecule'));
  org.helm.webeditor.MolViewer.molscale = 0.8;
  let app = new scil.helm.App(view, {
    showabout: false,
    mexfontsize: '90%',
    mexrnapinontab: true,
    topmargin: 20,
    mexmonomerstab: true,
    sequenceviewonly: false,
    mexfavoritefirst: true,
    mexfilter: true
  });
  var sizes = app.calculateSizes();
  app.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
  var s = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
  //@ts-ignore
  scil.apply(app.sequence.style, s);
  //@ts-ignore
  scil.apply(app.notation.style, s);
  s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + app.toolbarheight) + 'px'};
  //@ts-ignore
  scil.apply(app.properties.parent.style, s);
  app.structureview.resize(sizes.rightwidth, sizes.bottomheight + app.toolbarheight);
  app.mex.resize(sizes.topheight - 80);
  setTimeout(function() {
    if (typeof units === 'undefined') {
      app.canvas.helm.setSequence(cell.value, 'HELM');
    } else {
      app.canvas.helm.setSequence(value, 'HELM');
    }
  }, 200);
  //@ts-ignore
  ui.dialog({showHeader: false, showFooter: true})
    .add(view)
    .onOK(() => {
      let helmValue = app.canvas.getHelm(true).replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
      if (typeof units === 'undefined') {
        cell.value = helmValue;
      } else {
        let convertedRes = converter.convertHelmToFastaSeparator(helmValue, units);
        cell.value = convertedRes;
      }
    }).show({modal: true, fullScreen: true});
}

function getRS(smiles: string) {
  const newS = smiles.match(/(?<=\[)[^\][]*(?=])/gm);
  const res = {};
  let el = '';
  let digit;
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
  return res;
}

//name: getMolfiles
//input: column col {semType: Macromolecule}
//output: column res 
export function getMolfiles(col: DG.Column): DG.Column {
  let grid = grok.shell.tv.grid;
  let parent = grid.root.parentElement;
  let res = DG.Column.string('mols', col.length);
  const host = ui.div([]);
  parent.appendChild(host);
  let editor = new JSDraw2.Editor(host, {viewonly: true});
  host.style.width = '0px';
  host.style.height = '0px';
  res.init((i) => {
    editor.setHelm(col.get(i));
    let mol = editor.getMolfile();
    return mol;
  });
  parent.lastChild.remove();
  return res;
}

//name: helmWebEditor
//output: object
export function helmWebEditor(): HelmWebEditor {
  return new HelmWebEditor();
}
