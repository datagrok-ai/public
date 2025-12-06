/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {merge} from 'rxjs';

export const DSM_MOLECULE_ROW_SEMTYPE = 'DSM_MOLECULE_ROW';
export const getMolFuncTemp = 'getMolTempFunc';
export const DSM_MOLECULE_COLUMN_VID_TAG = 'DSM_MOLECULE_COLUMN_VID';

function getEmbeddings(_analysisName: string, rowCount: number) {
  // dummy function for getting embeddings
  console.log('Getting embeddings for analysis:', _analysisName);
  return {
    x: new Float32Array(rowCount).fill(0).map(() => Math.random()),
    y: new Float32Array(rowCount).fill(0).map(() => Math.random()),
  };
}

export async function chemSpaceApp(_analysisName: string) {
  const smilesCache = new DG.LruCache<number, string>(100);
  // _analysisName is not used for now
  const demoData = grok.data.demo.molecules(1000);
  const molecules = demoData.col('smiles')!.toList();
  // demo function
  const getMol = async (i: number): Promise<string> => {
    return smilesCache.getOrCreate(i, (k) => {
      console.log('Fetching molecule at index:', i);
      return molecules[k];
    });
  };
  const rowCount = molecules.length;
  const embeddings = getEmbeddings(_analysisName, rowCount);
  const xEmb = embeddings.x;
  const yEmb = embeddings.y;
  const props = ['logP', 'TPSA', 'MW', 'HBA', 'HBD'];
  const propVals = props.map((p) => new Float32Array(rowCount).fill(0).map(() => Math.random() * 10));

  const analysisDf = DG.DataFrame.create(rowCount);

  // add column that will have only row numbers in it, so that we can use it to fetch molecules
  const virtualMolCol = analysisDf.columns.addNewInt('Molecule');
  virtualMolCol.init((i) => i);
  virtualMolCol.semType = DSM_MOLECULE_ROW_SEMTYPE;
  virtualMolCol.setTag(DG.TAGS.CELL_RENDERER, DSM_MOLECULE_ROW_SEMTYPE);
  virtualMolCol.temp[getMolFuncTemp] = getMol;
  virtualMolCol.setTag(DSM_MOLECULE_COLUMN_VID_TAG, `chemspace-mol-vid-${Date.now()}-${Math.random()}`);

  analysisDf.columns.addNewFloat('x').init((u) => xEmb[u]);
  analysisDf.columns.addNewFloat('y').init((u) => yEmb[u]);
  props.forEach((p, i) => {
    analysisDf.columns.addNewFloat(p).init((u) => propVals[i][u]);
  });

  // add categorical columns
  const catColumnNames = ['Source', 'Class', 'Smell'];
  const catColumnValues = [
    ['PubChem', 'ChEMBL', 'ZINC', 'DrugBank', 'Other'],
    ['Alcohol', 'Ketone', 'Ester', 'Amine', 'Hydrocarbon', 'Acid', 'Other'],
    ['Pleasant', 'Unpleasant', 'Neutral', 'Strong', 'Smells like chicken', 'No smell']
  ];

  catColumnNames.forEach((colName, i) => {
    analysisDf.columns.addNewString(colName).init((u) => {
      const values = catColumnValues[i];
      return values[Math.floor(Math.random() * values.length)];
    });
  });

  // add hidden size column for scatterplot
  const sizeCol = analysisDf.columns.addNewFloat('~point_size').init((_i) => 6);

  const scatterplot = analysisDf.plot.scatter({
    x: 'x', y: 'y', color: 'Source', showXAxis: false, showYAxis: false,
    showXSelector: false, showYSelector: false, showSizeSelector: false,
    markerMinSize: 8, markerMaxSize: 8, sizeColumnName: '~point_size', markerDefaultSize: 8,
    rowSource: 'All', showTooltip: 'show custom tooltip', rowTooltip: 'Molecule',
    dataValues: 'Do not add'} as Partial<DG.IScatterPlotSettings>);


  // create grid to init renderer
  const tableView = DG.TableView.create(analysisDf, false);
  const view = ChemSpaceView.createWithGrid(null);


  view.name = 'ChemSpace Analysis';
  scatterplot.root.style.height = '100%';

  const filters = tableView.filters();


  view.root.appendChild(tableView.root);

  const p = view.getRibbonPanels();
  const downloadSelectionButton = ui.button('Download Selection', async () => {
    const selection = analysisDf.selection.clone().and(analysisDf.filter);
    if (selection.trueCount === 0) {
      grok.shell.warning('No filtered compounds selected for download.');
      return;
    }
    if (selection.trueCount > 20)
      grok.shell.warning('Maximum allowed compounds for download is 20. Reducing selection to first 20 compounds.');
    const newTable = DG.DataFrame.create(Math.min(selection.trueCount, 20));
    const indexes = selection.getSelectedIndexes();
    for (const col of analysisDf.columns) {
      if (col.name === 'x' || col.name === 'y' || col.name.startsWith('~'))
        continue;
      if (col.name === virtualMolCol.name) {
        const mols = await Promise.all(
          Array.from(indexes.slice(0, newTable.rowCount).map((i) => col.get(i) as number)).map(async (rowIdx) => await getMol(rowIdx)));
        newTable.columns.addNewString(col.name).init((i) => mols[i]);
      } else
        newTable.columns.addNew(col.name, col.type).init((i) => col.get(indexes[i]));
    }
    DG.Utils.download('selected-compounds.csv', newTable.toCsv());
  });

  let downloadTooltip = '';
  const updateDownloadButton = () => {
    const selectionCount = analysisDf.selection.clone().and(analysisDf.filter).trueCount;
    downloadTooltip = selectionCount > 0 ? `Download selected compounds as CSV (${selectionCount})` : 'No filtered compounds selected for download';
    downloadSelectionButton.disabled = selectionCount === 0;
    downloadSelectionButton.style.pointerEvents = 'auto';
  };
  updateDownloadButton();

  DG.debounce(analysisDf.onSelectionChanged, 10).subscribe(() => updateDownloadButton());
  ui.tooltip.bind(downloadSelectionButton, () => downloadTooltip);
  const minMarkerSizeInput = ui.input.slider('Filtered Out Size', {min: 1, max: 8, step: 1, value: 2,
    onValueChanged: (v) => {
      if (analysisDf.filter.anyFalse)
        scatterplot.props.markerMinSize = v;
    }});
  const maxMarkerSizeInput = ui.input.slider('Filtered Size', {min: 8, max: 20, step: 1, value: 8,
    onValueChanged: (v) => {
      if (analysisDf.filter.anyFalse)
        scatterplot.props.markerMaxSize = v;
      else {
        scatterplot.props.markerMinSize = v;
        scatterplot.props.markerMaxSize = v;
      }
    }});


  // when filter changes, update the size column
  DG.debounce(analysisDf.onFilterChanged, 50).subscribe(() => {
    const filter = analysisDf.filter;
    sizeCol.init((i) => filter.get(i) ? 6 : 1);
    if (!filter.anyFalse)
      scatterplot.props.markerMinSize = maxMarkerSizeInput.value!;
    else
      scatterplot.props.markerMinSize = minMarkerSizeInput.value!;
    scatterplot.props.markerMaxSize = maxMarkerSizeInput.value!;
  });

  DG.debounce(analysisDf.onCurrentRowChanged, 300).subscribe(() => {
    const currRow = analysisDf.currentRowIdx;
    if (currRow == null || currRow < 0 || currRow >= analysisDf.rowCount)
      return;
    getMol(virtualMolCol.get(currRow) as number).then((mol) => {
      if (mol)
        grok.shell.o = DG.SemanticValue.fromValueType(mol, DG.SEMTYPE.MOLECULE);
    }).catch((err) => console.error(err));
  });

  p.push([minMarkerSizeInput.root, maxMarkerSizeInput.root, downloadSelectionButton]);


  view.setRibbonPanels(p);
  grok.shell.addView(view);
  tableView._onAdded();
  view.grid = tableView.grid;
  tableView.dockManager.close(tableView.grid.root);
  const scNode = tableView.dockManager.dock(scatterplot, DG.DOCK_TYPE.FILL, null, 'Space');
  tableView.dockManager.dock(filters, DG.DOCK_TYPE.LEFT, tableView.dockManager.findNode(scatterplot.root), 'Filters', 0.2);
  const formsViewer = tableView.addViewer('Forms', {fieldsColumnNames: analysisDf.columns.names().filter((n) => n !== 'x' && n !== 'y' && !n.startsWith('~')), moleculeSize: 'normal', showCurrentRow: false, showMouseOverRow: false});
  console.log(formsViewer);
  tableView.dockManager.dock(formsViewer, DG.DOCK_TYPE.DOWN, scNode, 'Details', 0.3);

  // clicking without ctrl should also select the point
  analysisDf.onCurrentRowChanged.subscribe(() => {
    const currRow = analysisDf.currentRowIdx;
    if (currRow == null || currRow < 0 || currRow >= analysisDf.rowCount)
      return;
    const isSelected = analysisDf.selection.get(currRow);
    analysisDf.selection.set(currRow, !isSelected);
  });

  // forms viewer needs to have its own table. there are some hacks in, but don't mind them for now
  const resetFormsViewerTable = () => {
    const selectionWithFilter = analysisDf.selection.clone().and(analysisDf.filter);
    // limit to 20 rows
    const fvTable = analysisDf.clone(selectionWithFilter);
    if (selectionWithFilter.trueCount > 20)
      fvTable.rows.removeWhereIdx((i) => i >= 20);

    // set selection to all
    fvTable.selection.setAll(true);
    if (fvTable.rowCount === 0)
      fvTable.rows.addNew(); // forms viewer cannot work with empty tables ....

    fvTable.col(virtualMolCol.name)!.semType = DSM_MOLECULE_ROW_SEMTYPE;
    fvTable.col(virtualMolCol.name)!.setTag(DG.TAGS.CELL_RENDERER, DSM_MOLECULE_ROW_SEMTYPE);
    fvTable.col(virtualMolCol.name)!.temp[getMolFuncTemp] = getMol;
    fvTable.col(virtualMolCol.name)!.setTag(DSM_MOLECULE_COLUMN_VID_TAG, virtualMolCol.getTag(DSM_MOLECULE_COLUMN_VID_TAG));
    //formsViewer.dataFrame = fvTable;

    const g = fvTable.plot.grid();
    if ((formsViewer as any)?.children[0])
      (formsViewer as any).children[0].getGrid = () => g;
    formsViewer.dataFrame = fvTable;
  };

  //resetFormsViewerTable();
  DG.debounce(merge(analysisDf.onFilterChanged, analysisDf.onSelectionChanged), 50).subscribe(() => resetFormsViewerTable());
}

export class DSMMoleculeRenderer extends DG.GridCellRenderer {
  get name() {
    return DSM_MOLECULE_ROW_SEMTYPE;
  }

  get cellType() {
    return DSM_MOLECULE_ROW_SEMTYPE;
  }

  get defaultHeight(): number | null {
    return 300;
  }
  get defaultWidth(): number | null {
    return 300;
  }

  private _MolCache = new Map<string, DG.LruCache<number, string>>();

  private _renderQueue: EvictingArray<() => Promise<void>> = new EvictingArray(20);

  private _renderTimer: any = null;
  private renderDebounced(ctx: CanvasRenderingContext2D, rFunc: () => Promise<void>) {
    // draw a loader while rendering
    ctx.save();
    ctx.fillStyle = 'lightgray';
    ctx.font = '16px Arial';
    ctx.fillText('Fetching...', 10, 20);
    ctx.restore();
    this._renderQueue.push(rFunc);
    if (this._renderTimer)
      clearTimeout(this._renderTimer);
    this._renderTimer = setTimeout(async () => {
      const renderFuncs = this._renderQueue.toArray();
      this._renderQueue.empty();
      for (const f of renderFuncs)
        await f();
    }, 300);
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, value: DG.GridCell, context: any): void {
    // only render if the rendering happens in the tooltip or in the forms, i.e self widgets. this can easily be checked by looking at the rendering size
    const canvas = g.canvas;
    const table = value.cell.dataFrame;
    const dpr = window.devicePixelRatio;
    if (x > 2 * dpr || y > 2 * dpr || w * dpr < (canvas.width - 4) || h * dpr < (canvas.height - 4) || !table || !value.cell.column.temp[getMolFuncTemp] || !value.cell.column.getTag(DSM_MOLECULE_COLUMN_VID_TAG))
      return;

    const val = value.cell.value;
    const mapKey = value.cell.column.getTag(DSM_MOLECULE_COLUMN_VID_TAG);
    if (!this._MolCache.has(mapKey))
      this._MolCache.set(mapKey, new DG.LruCache<number, string>(100));
    const molCache = this._MolCache.get(mapKey)!;
    const mol = molCache.get(val);
    // const mol = molCache.getOrCreate(val, (k) => value.cell.column.temp[getMolFuncTemp](k));
    if (mol !== undefined && mol !== null)
      grok.chem.canvasMol(x, y, w, h, canvas, mol, null);
    else {
      this.renderDebounced(g, async () => {
        if (!document.contains(canvas))
          return;
        const molFetched = await value.cell.column.temp[getMolFuncTemp](val);
        molCache.set(val, molFetched ?? '');
        grok.chem.canvasMol(x, y, w, h, canvas, molFetched, null);
      });
    }
  }
}

class EvictingArray<T> {
  private _capacity: number;
  private _array: T[] = [];
  constructor(capacity: number) {
    this._capacity = capacity;
  }
  push(item: T) {
    if (this._array.length >= this._capacity)
      this._array.shift();
    this._array.push(item);
  }
  toArray(): T[] {
    return this._array;
  }

  empty(): void {
    this._array = [];
  }
}

class ChemSpaceView extends DG.View {
  public grid: DG.Grid | null = null;

  static createWithGrid(grid: DG.Grid | null): ChemSpaceView {
    const view = super.create();
    (view as ChemSpaceView).grid = grid;
    return view as ChemSpaceView;
  }
}
