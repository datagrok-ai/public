/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


export const _package = new DG.Package();

//tags: fileViewer, fileViewer-nwk
//input: file file
//output: view v
export async function nwkTreeViewer(file) {
  const newick = await file.readAsString();

  const view = DG.View.create();
  const host = ui.div([
    ui.button('Load dataframe', () => exportToDf(file), 'Export to dataframe'),
    (await grok.functions.call('TreeViewer:treePreview', { newick })) // DG.Viewer
  ], 'd4-ngl-viewer');

  view.append(host);
  return view;
}

//input: string s {semType: newick}
//output: view v
export function newickTableView(s) {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'node', ['node']),
    DG.Column.fromList('string', 'parent', ['parent']),
    DG.Column.fromList('double', 'distance', [0.5]),
  ]);
  return DG.TableView.create(df);
}

// https://github.com/jasondavies/newick.js
function parseNewick(a) {for(var e=[],r={},s=a.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){var n=s[t];switch(n){case"(":var c={};r.branchset=[c],e.push(r),r=c;break;case",":var child={};e[e.length-1].branchset.push(child),r=child;break;case")":r=e.pop();break;case":":break;default:var h=s[t-1];")"==h||"("==h||","==h?r.name=n:":"==h&&(r.length=parseFloat(n))}}return r};

function newickToDf(newick, filename) {
  let parent = 'root';
  let i = 0;
  const obj = parseNewick(newick);
  const nodes = [], parents = [], distances = [];

  function traverse(obj) {
    if (obj === null || typeof obj != "object" ) return;
    if (!Array.isArray(obj)) {
      let name = obj.name;
      if (!name) name = `node-${i}`, i += 1;
      nodes.push(name);
      distances.push(obj.length);
      parents.push(parent);
      parent = name;
    }
    Object.values(obj).forEach(value => traverse(value));
  }
  traverse(obj);

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'node', nodes),
    DG.Column.fromList('string', 'parent', parents),
    DG.Column.fromList('double', 'distance', distances),
  ]);
  df.name = `df-${filename.slice(0, -4)}`;
  return df;
};


//input: file file
//output: view view
export async function exportToDf(file) {
  const currentFolderPath = file.fullPath.slice(0, -file.name.length);
  const fileInfos = (await grok.dapi.files.list(currentFolderPath, false, '')).filter(f => f.extension === 'nwk');
  const newickStrings = await Promise.all(fileInfos.map(f => f.readAsString()));

  const rowCount = newickStrings.length;
  const fileNameCol = DG.Column.fromList('string', 'file', fileInfos.map(f => f.fileName));
  const newickCol = DG.Column.fromList('string', 'newick', newickStrings);
  const dfCol = DG.Column.dataFrame('dataframe', rowCount);
  const linkCol = DG.Column.string('links', rowCount);
  const df = DG.DataFrame.fromColumns([fileNameCol, newickCol, dfCol, linkCol]);

  newickCol.semType = 'newick';
  for (let i = 0; i < rowCount; i++) dfCol.set(i, newickToDf(newickCol.get(i), fileNameCol.get(i)));

  const view = grok.shell.addTableView(df);
  view.grid.columns.byName('links').cellType = 'html';
  view.grid.onCellPrepare(gc => {
    if (gc.isTableCell && gc.gridColumn.name === 'links') {
      gc.style.element = ui.button('open', () => grok.shell.addTableView(DG.toJs(gc.tableRow.dataframe)), 'Open dataframe');
    }
  });
}
