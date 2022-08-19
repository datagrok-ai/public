/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NotationConverter} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {createJsonMonomerLibFromSdf} from './utils';
import {MONOMER_MANAGER_MAP, RGROUPS, RGROUP_CAP_GROUP_NAME, RGROUP_LABEL, SMILES} from './constants';
import {printLeftOrCentered} from '@datagrok/bio/src/utils/cell-renderer';

export const _package = new DG.Package();

const lru = new DG.LruCache<any, any>();
const molfiles = new DG.LruCache<any, any>();
const STORAGE_NAME = 'Libraries';
let i = 0;
const LIB_PATH = 'libraries/';


//tags: init
export async function initHelm(): Promise<void> {
  let grid = grok.shell.tv.grid;
  for (let j = 0; j <= i; j++) {
    await monomerManager(await grok.dapi.userDataStorage.getValue(STORAGE_NAME, j.toString(), true));
  }
  grid.invalidate();
  
  return new Promise((resolve, reject) => {
    // @ts-ignore
    dojo.ready(function() { resolve(null); });
  });
}

//name: helmCellRenderer
//tags: cellRenderer
//meta.cellType: helm
//meta.columnTags: units=helm
//output: grid_cell_renderer result
export function helmCellRenderer(): HelmCellRenderer {
  return new HelmCellRenderer();
}


export function webEditor(value: string) {
  let view = ui.div();
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
    app.canvas.helm.setSequence(value, 'HELM');
  }, 200);
  //@ts-ignore
  ui.dialog({showHeader: false, showFooter: true})
  .add(view)
  .onOK(() => {
    const val = app.canvas.getHelm(true).replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
    value = val;
  }).show({modal: true, fullScreen: true});
}


//tags: cellEditor
//description: Macromolecule
//input: grid_cell cell
export function editMoleculeCell(cell: DG.GridCell): void {
  let grid = grok.shell.tv.grid;
  webEditor(cell.cell.value);
  grid.invalidate();
}

//name: Open Helm Web Editor
//description: Adds editor
//meta.action: Open Helm Web Editor
//input: string mol { semType: Macromolecule }
export function openEditor(mol: string): void {
  let df = grok.shell.tv.grid.dataFrame;
  let converter = new NotationConverter(df.col('sequence'));
  const resStr = converter.convertStringToHelm(mol, '/');
  webEditor(resStr);
}

//name: Properties
//tags: panel, widgets
//input: string helmString {semType: Macromolecule}
//output: widget result
export async function propertiesPanel(helmString: string) {
  //const lru = await getLru();
  const result = lru.get(helmString).split(',');
  return new DG.Widget(
    ui.tableFromMap({
      'formula': result[0].replace(/<sub>/g, '').replace(/<\/sub>/g, ''),
      'molecular weight': result[1],
      'extinction coefficient': result[2],
      /*'fasta': ui.wait(async () => ui.divText(await helmToFasta(helmString))),
      'rna analogue sequence': ui.wait(async () => ui.divText(await helmToRNA(helmString))),
      'smiles': ui.wait(async () => ui.divText(await helmToSmiles(helmString))),
      //'peptide analogue sequence': ui.wait(async () => ui.divText(await helmToPeptide(helmString))),*/
    })
  );
}

//name: Molfile
//tags: panel, widgets
//input: string helmString {semType: Macromolecule}
//output: widget result
export async function molfilePanel(helmString: string) {
  return ui.textInput('', await getMolfile(helmString));
}


//name: loadDialog
export async function loadDialog() {
  let res = (await _package.files.list(`${LIB_PATH}`, false, '')).map(it => it.fileName);
  let FilesList = await ui.choiceInput('Monomer Libraries', ' ', res);
  let grid = grok.shell.tv.grid;
  ui.dialog('Load library from file')
    .add(FilesList)
    .onOK(async () => {
      await monomerManager(FilesList.value);
      grok.dapi.userDataStorage.postValue(STORAGE_NAME, i.toString(), FilesList.value, true);
      i += 1;
      grid.invalidate();
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
  //@ts-ignore
  let filesButton = ui.button('Manage');
  filesButton.addEventListener('click', manageFiles);
  let list = (await _package.files.list(`${LIB_PATH}`, false, '')).map(it => it.fileName);
  let usedLibraries = [];
  for (let j = 0; j <= i; j++) {
    usedLibraries.push(await grok.dapi.userDataStorage.getValue(STORAGE_NAME, j.toString(), true));
  }
  let unusedLibraries = list.filter(x => !usedLibraries.includes(x));
  return new DG.Widget(ui.splitV([
    ui.splitH([
      ui.tableFromMap({
        'Uploaded libraries': ui.divText(usedLibraries.toString()),
        'Not used libraries': ui.divText(unusedLibraries.toString()),
      })
    ]),
    ui.splitH([
      ui.divH([loadButton]),
      ui.divH([filesButton]),
    ])
  ]));
}

//name: manageFiles
export async function manageFiles() {
  const a = ui.dialog({title: 'Manage files'})
    //@ts-ignore
    .add(ui.fileBrowser({path: 'System:AppData/Helm/libraries'}).root)
    .addButton('OK', () => a.close())
    .show();
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

//name: monomerManager
//input: string value
export async function monomerManager(value: string) {
  let df: any[];
  let file;
  let dfSdf;
  if (value.endsWith('.sdf')) {
    const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'importSdf'});
    console.debug(`Helm: initHelm() funcList.length = ${funcList.length}`);
    if (funcList.length === 1) {
      file = await _package.files.readAsBytes(`${LIB_PATH}${value}`);
      dfSdf = await grok.functions.call('Chem:importSdf', {bytes: file});
      df = createJsonMonomerLibFromSdf(dfSdf[0]);
    } else {
      grok.shell.warning("Chem package is not installed");
    }
  } else {
    const file = await _package.files.readAsText(`${LIB_PATH}${value}`);
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
  });
}

//name: helmColumnToSmiles
//input: column helmColumn {semType: Macromolecule}
export function helmColumnToSmiles(helmColumn: DG.Column) {
  //todo: add column with smiles to col.dataFrame.
}

//name: getMolfile
//input: string helmString {semType: Macromolecule}
export async function getMolfile(helmString) {
  return molfiles.get(helmString);
}

function split(s: string, sep: string) {
  var ret = [];
  var frag = "";
  var parentheses = 0;
  var bracket = 0;
  var braces = 0;
  var quote = 0;
  for (var i = 0; i < s.length; ++i) {
      var c = s.substring(i, i + 1);
      if (c == sep && bracket == 0 && parentheses == 0 && braces == 0 && quote == 0) {
          ret.push(frag);
          frag = "";
      }
      else {
          frag += c;
          if (quote > 0) {
              if (c == '\\' && i + 1 < s.length) {
                  ++i;
                  var c2 = s.substring(i, i + 1);
                  frag += c2;
                  c += c2;
              }
          }
          if (c == '\"') {
              if (!(i > 0 && s.substring(i - 1, i) == '\\'))
                  quote = quote == 0 ? 1 : 0;
          }
          else if (c == '[')
              ++bracket;
          else if (c == ']')
              --bracket;
          else if (c == '(')
              ++parentheses;
          else if (c == ')')
              --parentheses;
          else if (c == '{')
              ++braces;
          else if (c == '}')
              --braces;
      }
  }
  ret.push(frag);
  return ret;
}

function detachAnnotation(s: string) {
  var ret = _detachAppendix(s, '\"');
  if (ret.tag != null)
      return ret;

  var r = _detachAppendix(s, '\'');
  return { tag: ret.tag, repeat: r.tag, str: r.str };
}

function _detachAppendix(s: string, c: string) {
  var tag = null;
  //@ts-ignore
  if (scil.Utils.endswith(s, c)) {
      var p = s.length - 1;
      while (p > 0) {
          p = s.lastIndexOf(c, p - 1);
          if (p <= 0 || s.substring(p - 1, p) != '\\')
              break;
      }

      if (p > 0 && p < s.length - 1) {
          tag = s.substring(p + 1, s.length - 1);
          s = s.substring(0, p);
      }
  }
  if (tag != null)
      tag = tag.replace(new RegExp("\\" + c, "g"), c);
  return { tag: unescape(tag), str: s };
}

function unescape(s: string) {
  //@ts-ignore
  if (scil.Utils.isNullOrEmpty(s))
      return s;

  return s.replace(/[\\]./g, function (m) {
      switch (m) {
          case "\\r":
              return "\r";
          case "\\n":
              return "\n";
          case "\\t":
              return "\t";
          default:
              return m.substring(1);
      }
  });
}

function parseHelm(s: string) {
  var sections = split(s, '$');
  s = sections[0];
  var monomers = [];
  //@ts-ignore
  if (!scil.Utils.isNullOrEmpty(s)) {
      var seqs = split(s, '|');
      for (var i = 0; i < seqs.length; ++i) {
          var e = detachAnnotation(seqs[i]);
          s = e.str;

          var p = s.indexOf("{");
          
          s = s.substring(p + 1);
          p = s.indexOf('}');
          s = s.substring(0, p);

          var ss = split(s, '.');
          for (var monomer of ss) {
            if (monomer.includes('(') && monomer.includes(')')) {
                var elements = monomer.replace(/[()]/g, "").split("");
                for (var el of elements) {
                    monomers.push(el);
                }
            } else {
                monomers.push(monomer);
            }
          }
      }
  }
  return monomers;
}

//name: findMonomers
//input: string helmString
export function findMonomers(helmString: string) {
  //@ts-ignore
  const types = Object.keys(org.helm.webeditor.monomerTypeList());
  const monomers: any = [];
  const monomer_names: any = [];
  for (var i = 0; i < types.length; i++) {
    //@ts-ignore
    monomers.push(new scil.helm.Monomers.getMonomerSet(types[i]));
    Object.keys(monomers[i]).forEach(k => {
      monomer_names.push(monomers[i][k].id);
    });
  }
  const split_string = parseHelm(helmString);
  return new Set(split_string.filter(val => !monomer_names.includes(val)));
}

class HelmCellRenderer extends DG.GridCellRenderer {
  get name() { return 'helm'; }

  get cellType() { return 'helm'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
    const undefinedColor = 'rgb(100,100,100)';
    const grayColor = '#808080';
    const monomers = findMonomers(gridCell.cell.value);
    if (monomers.size == 0) {
      const host = ui.div([], {style: {width: `${w}px`, height: `${h}px`}});
      host.setAttribute('dataformat', 'helm');
      host.setAttribute('data', gridCell.cell.value);
      gridCell.element = host;
        //@ts-ignore
        const canvas = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true});
        const formula = canvas.getFormula(true);
        if (!formula) {
          gridCell.element = ui.divText(gridCell.cell.value, {style: {color: 'red'}});
        }
        const molWeight = Math.round(canvas.getMolWeight() * 100) / 100;
        const coef = Math.round(canvas.getExtinctionCoefficient(true) * 100) / 100;
        const molfile = canvas.getMolfile();
        const result = formula + ', ' + molWeight + ', ' + coef;
        lru.set(gridCell.cell.value, result);
        molfiles.set(gridCell.cell.value, molfile);
        return;
      }
      if (monomers.size > 0) {
        w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
        g.save();
        g.beginPath();
        g.rect(x, y, w, h);
        g.clip();
        g.font = '12px monospace';
        g.textBaseline = 'top';
        let x1 = x;
        const s: string = gridCell.cell.value ?? '';
        let subParts: string[] = parseHelm(s);
        subParts.forEach((amino, index) => {
          let color = monomers.has(amino) ? 'red' : grayColor;
          g.fillStyle = undefinedColor;
          let last = index === subParts.length - 1;
          x1 = printLeftOrCentered(x1, y, w, h, g, amino, color, 0, true, 1.0, '/', last);
        });
        g.restore();
        return;
      }
  }
}

//name: getMonomerLib
//input: string type
//output: string monomerSetJSON
export function getMonomerLib(type: string): string {
  return JSON.stringify(scil.helm.Monomers.getMonomerSet(type));
}
