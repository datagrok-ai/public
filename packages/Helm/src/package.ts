/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

const lru = new DG.LruCache<any, any>(); 

//tags: init
export async function initChem(): Promise<void> {
  // apparently HELMWebEditor requires dojo to be initialized first
  return new Promise((resolve, reject) => {
    // @ts-ignore
    dojo.ready(function () { resolve(null); });
  });
}


//name: helmCellRenderer
//tags: cellRenderer,cellRenderer-HELM
//meta.cellType: HELM
//meta-cell-renderer-sem-type: HELM
//output: grid_cell_renderer result
export function helmCellRenderer(): DG.GridCellRenderer {
  return new HelmCellRenderer();
}

//tags: cellEditor
//description: HELM
//input: grid_cell cell
export function editMoleculeCell(cell: DG.GridCell): void {
  let view = ui.div();
  //@ts-ignore
  org.helm.webeditor.MolViewer.molscale = 0.8;
  //@ts-ignore
  let app = new scil.helm.App(view, { showabout: false, mexfontsize: "90%", mexrnapinontab: true, topmargin: 20, mexmonomerstab: true, sequenceviewonly: false, mexfavoritefirst: true, mexfilter: true });
  setTimeout(function() {
    //@ts-ignore
    app.canvas.helm.setSequence(cell.cell.value , 'HELM');
  }, 200);
  //@ts-ignore
  ui.dialog({ showHeader: false, showFooter: true })
  .add(view)
  .onOK(() => {
    cell.cell.value = app.canvas.getHelm(true);
  }).show({ modal: true, fullScreen: true});
}

//name: Details
//tags: panel, widgets
//input: string helmString {semType: HELM}
//output: widget result
export function detailsPanel(helmString: string){
  //return new DG.Widget(ui.divText(lru.get(helmString)));
  var result = lru.get(helmString).split(',');
  return new DG.Widget(
    ui.tableFromMap({
      'formula': result[0].replace(/<sub>/g, '').replace(/<\/sub>/g, ''),
      'molecular weight': result[1],
      'extinction coefficient': result[2],
      'molfile': result[3],
      'fasta': ui.wait(async () => ui.divText(await helmToFasta(helmString))),
      'rna analogue sequence': ui.wait(async () => ui.divText(await helmToRNA(helmString))),
      'smiles': ui.wait(async () => ui.divText(await helmToSmiles(helmString))),
      'peptide analogue sequence': ui.wait(async () => ui.divText(await helmToPeptide(helmString))),
    })
  )
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
//input: string helmString {semType: HELM}
//output: string res
export async function helmToFasta(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Produce/${helmString}`
  return await accessServer(url, 'FastaFile');
}

//name: helmToRNA
//description: converts to rna analogue sequence
//input: string helmString {semType: HELM}
//output: string res
export async function helmToRNA(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Convert/RNA/${helmString}`
  return await accessServer(url, 'Sequence');
}


//name: helmToPeptide
//description: converts to peptide analogue sequence
//input: string helmString {semType: HELM}
//output: string res
export async function helmToPeptide(helmString: string) {
  const url = `http://localhost:8081/WebService/service/Fasta/Convert/PEPTIDE/${helmString}`
  return await accessServer(url, 'Sequence');
}

//name: helmToSmiles
//description: converts to smiles
//input: string helmString {semType: HELM}
//output: string smiles {semType: Molecule}
export async function helmToSmiles(helmString: string) {
  const url = `http://localhost:8081/WebService/service/SMILES/${helmString}`
  return await accessServer(url, 'SMILES');
}


//name: helmColumnToSmiles
//input: column helmColumn {semType: HELM}
export function helmColumnToSmiles(helmColumn: DG.Column) {
  //todo: add column with smiles to col.dataFrame.
}


class HelmCellRenderer extends DG.GridCellRenderer {

  get name() { return 'helm'; }
  get cellType() { return 'helm'; }
  get defaultWidth(): number | null { return 400; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    let host = ui.div([], { style: { width: `${w}px`, height: `${h}px`}});
    host.setAttribute('dataformat', 'helm');
    host.setAttribute('data', gridCell.cell.value);

    gridCell.element = host;
    // @ts-ignore
    var canvas = new JSDraw2.Editor(host, { width: w, height: h, skin: "w8", viewonly: true });
    var formula = canvas.getFormula(true);
    var molWeight = Math.round(canvas.getMolWeight() * 100) / 100;
    var coef = Math.round(canvas.getExtinctionCoefficient(true) * 100) / 100;
    var molfile = canvas.getMolfile();
    var result = formula + ', ' + molWeight + ', ' + coef + ', ' + molfile;
    lru.set(gridCell.cell.value, result);
  }
}
