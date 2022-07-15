/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//import {lru} from '../../Bio/src/utils/cell-renderer';
import {NotationConverter, NOTATION} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import { createJsonMonomerLibFromSdf } from './utils';
import { MONOMER_MANAGER_MAP, RGROUPS, RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES } from './constants';
//import {ConverterFunc, DfReaderFunc} from '../../Bio/src/tests/types';

export const _package = new DG.Package();

const lru = new DG.LruCache<any, any>();
const STORAGE_NAME = 'Libraries';
let i = 0;


//tags: init
export async function initChem(): Promise<void> {
  // apparently HELMWebEditor requires dojo to be initialized first
  return new Promise((resolve, reject) => {
    // @ts-ignore
    dojo.ready(function() { resolve(null); });
  });
}

//tags: cellEditor
//description: Macromolecule
//input: grid_cell cell
export function editMoleculeCell(cell: DG.GridCell): void {
  const view = ui.div();
  org.helm.webeditor.MolViewer.molscale = 0.8;
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
  app.canvas.resize(sizes.rightwidth, sizes.topheight - 210);
  let s = {width: sizes.rightwidth + 'px', height: sizes.bottomheight + 'px'};
  //@ts-ignore
  scil.apply(app.sequence.style, s);
  //@ts-ignore
  scil.apply(app.notation.style, s);
  s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + app.toolbarheight) + 'px'};
  //@ts-ignore
  scil.apply(app.properties.parent.style, s);
  app.structureview.resize(sizes.rightwidth, sizes.bottomheight + app.toolbarheight);
  app.mex.resize(sizes.topheight - 80);

  let dialogValue: string;
  if (cell.gridColumn.column.tags[DG.TAGS.UNITS] === 'HELM') {
    dialogValue = cell.cell.value;
  } else {
    const converter: NotationConverter = new NotationConverter(cell.cell.column);
    dialogValue = converter.convertStringToHelm(cell.cell.value);
  }

  setTimeout(function() {
    app.canvas.helm.setSequence(dialogValue, 'HELM');
  }, 200);
  //@ts-ignore
  ui.dialog({showHeader: false, showFooter: true})
    .add(view)
    .onOK(() => {
      cell.cell.value;
    }).show({modal: true, fullScreen: true});
}

//name: Details
//tags: panel, widgets
//input: string helmString {semType: Macromolecule}
//output: widget result
export function detailsPanel(helmString: string) {
  //return new DG.Widget(ui.divText(lru.get(helmString)));
  const result = lru.get(helmString).split(',');
  return new DG.Widget(
    ui.tableFromMap({
      'formula': result[0].replace(/<sub>/g, '').replace(/<\/sub>/g, ''),
      'molecular weight': result[1],
      'extinction coefficient': result[2],
      'molfile': result[3],
      'fasta': ui.wait(async () => ui.divText(await helmToFasta(helmString))),
      'rna analogue sequence': ui.wait(async () => ui.divText(await helmToRNA(helmString))),
      'smiles': ui.wait(async () => ui.divText(await helmToSmiles(helmString))),
      //'peptide analogue sequence': ui.wait(async () => ui.divText(await helmToPeptide(helmString))),
    })
  );
}


async function loadDialog () {
  //@ts-ignore
  let res = await grok.dapi.files.list('System:AppData/Helm/libraries', false, '');
  //@ts-ignore
  res = res.map((e) => e.path);
  //@ts-ignore
  let FilesList = await ui.choiceInput('Monomer Libraries', ' ', res);
  let grid = grok.shell.tv.grid;
  ui.dialog('Load library from file')
  .add(FilesList)
  .onOK(async () => {
    await monomerManager(FilesList.value);
    grok.dapi.userDataStorage.postValue(STORAGE_NAME, i.toString(), FilesList.value, true);
    i += 1;
    grid.invalidate()
  }).show();
};

//name: Manage Libraries
//tags: panel, widgets
//input: column helmColumn {semType: Macromolecule}
//output: widget result
export async function libraryPanel(helmColumn: DG.Column) {
  //@ts-ignore
  let loadButton = ui.button('Load Library');
  loadButton.addEventListener('click', loadDialog);
  let list = (await grok.dapi.files.list('System:AppData/Helm/libraries', false, '')).toString().split(',');
  let usedLibraries = [];
  for (let j = 0; j <= i; j++) {
    usedLibraries.push(await grok.dapi.userDataStorage.getValue(STORAGE_NAME, j.toString(), true));
  }
  let unusedLibraries = list.filter(x => !usedLibraries.includes(x));
  return new DG.Widget(
    ui.tableFromMap({
      'Uploaded libraries': ui.divText(usedLibraries.toString()),
      'Not used libraries': ui.divText(unusedLibraries.toString()),
      'Upload new library': ui.divH([loadButton]),
    })
  );
}

async function accessServer(url: string, key: string) {
  const params: RequestInit = {
    method: 'GET',
    headers: {
      'Accept': 'application/json',
    }
  };
  const response = await fetch(url, params);
  const json = await response.json();
  return json[key];
}


//name: helmToFasta
//input: string helmString {semType: Macromolecule}
//output: string res
export async function helmToFasta(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Produce/${helmString}`;
  return await accessServer(url, 'FastaFile');
}

//name: helmToRNA
//description: converts to rna analogue sequence
//input: string helmString {semType: Macromolecule}
//output: string res
export async function helmToRNA(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Convert/RNA/${helmString}`;
  return await accessServer(url, 'Sequence');
}


//name: helmToPeptide
//description: converts to peptide analogue sequence
//input: string helmString {semType: Macromolecule}
//output: string res
export async function helmToPeptide(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Convert/PEPTIDE/${helmString}`;
  return await accessServer(url, 'Sequence');
}

//name: helmToSmiles
//description: converts to smiles
//input: string helmString {semType: Macromolecule}
//output: string smiles {semType: Molecule}
export async function helmToSmiles(helmString: string) {
  const url = `http://localhost:8081/WebService/service/SMILES/${helmString}`;
  return await accessServer(url, 'SMILES');
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

async function monomerManager(value: string) {
  let df: any[];
  if(value.endsWith('.sdf')) {
    const file = await _package.files.readAsBytes(value.split('/')[1]);
    const dfSdf = await grok.functions.call('Chem:importSdf', {bytes: file});
    df = createJsonMonomerLibFromSdf(dfSdf[0]);

  } else {
    const file = await _package.files.readAsText(value.split('/')[1]);
    df = JSON.parse(file);
  }
  
  df.forEach(monomer => {
    const m = {};
    Object.keys(MONOMER_MANAGER_MAP).forEach(field => {
      m[field] = monomer[MONOMER_MANAGER_MAP[field]] ?? '';
    });
    if (monomer[RGROUPS]) {
      m['rs'] = monomer[RGROUPS].length;
      const at = {};
      monomer[RGROUPS].forEach(it => {
        at[it[RGROUP_LABEL]] = it[RGROUP_CAP_GROUP_NAME];
      });
      m['at'] = at;
    } else {
      m['rs'] = Object.keys(getRS(df[SMILES].toString())).length;
      m['at'] = getRS(df[SMILES].toString());
    }
    org.helm.webeditor.Monomers.addOneMonomer(m);
  })
}

//name: helmColumnToSmiles
//input: column helmColumn {semType: Macromolecule}
export function helmColumnToSmiles(helmColumn: DG.Column) {
  //todo: add column with smiles to col.dataFrame.
}

// Synonym to overcome eslint error
const GetMonomerSet = scil.helm.Monomers.getMonomerSet;

//name: findMonomers
//input: string helmString { semType: Macromolecule }
export async function findMonomers(helmString: string) {
  const types = Object.keys(org.helm.webeditor.monomerTypeList());
  const monomers = [];
  const monomerNames = [];
  for (let i = 0; i < types.length; i++) {
    monomers.push(new GetMonomerSet(types[i]));
    Object.keys(monomers[i]).forEach((k) => {
      monomerNames.push(monomers[i][k].id);
    });
  }
  const splitString = WebLogo.splitterAsHelm(helmString);
  return new Set(splitString.filter((val) => !monomerNames.includes(val)));
}

class HelmCellRenderer extends DG.GridCellRenderer {
  get name() { return 'macromolecule'; }

  get cellType() { return 'macromolecule'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const host = ui.div([], {style: {width: `${w}px`, height: `${h}px`}});
    host.setAttribute('dataformat', 'helm');
    host.setAttribute('data', gridCell.cell.value);

    gridCell.element = host;
    const canvas = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true});
    const formula = canvas.getFormula(true);
    const molWeight = Math.round(canvas.getMolWeight() * 100) / 100;
    const coef = Math.round(canvas.getExtinctionCoefficient(true) * 100) / 100;
    const molfile = canvas.getMolfile();
    const result = formula + ', ' + molWeight + ', ' + coef + ', ' + molfile;
    lru.set(gridCell.cell.value, result);
  }
}
