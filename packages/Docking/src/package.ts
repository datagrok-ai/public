/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '@datagrok-libraries/bio/src/types/ngl'; // To enable import from the NGL module declared in bio lib
import {GridSize, IAutoDockService} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {AutoDockApp, AutoDockDataType} from './apps/auto-dock-app';
import {_runAutodock, AutoDockService, _runAutodock2} from './utils/auto-dock-service';
import {_package, TARGET_PATH, CACHED_DOCKING, CACHED_MOLSTAR, AFFINITY_COL, POSE_COL} from './utils/constants';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

// -- Services & Helpers --

//name: getAutoDockService
//output: object result
export async function getAutoDockService(): Promise<IAutoDockService> {
  const resSvc: IAutoDockService = await AutoDockService.getSvc();
  return resSvc;
}

// -- Test apps --

//name: autoDockApp
export async function autoDockApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('AutoDock app...');
  try {
    const app = new AutoDockApp('autoDockApp');
    await app.init();
  } finally {
    pi.close();
  }
}

// -- AutoDock --

//name: runAutodock
//input: file receptor
//input: file ligand
//input: int x
//input: int y
//input: int z
//output: dataframe result
export async function runAutodock(
  receptor: DG.FileInfo, ligand: DG.FileInfo, x: number, y: number, z: number,
): Promise<DG.DataFrame | null> {
  const npts = new GridSize(x, y, z);
  return await _runAutodock(receptor, ligand, npts);
}

//name: runAutodock2
// //input: dataframe df
// //input: column molCol { semType: Molecule }
// //input: file receptorFi { optional: true }
export async function runAutodock2(): Promise<void> {
  const [csv, receptorPdb] = await Promise.all([
    grok.functions.call(`${_package.name}:readAsText`,
      {file: 'CHEMBL2366517/ic50.mol.csv'}),
    grok.functions.call(`${_package.name}:readAsText`, {file: 'CHEMBL2366517/1bdq.pdb'}),
  ]);
  if (!csv || !receptorPdb)
    throw new Error('Empty input data');

  const df = DG.DataFrame.fromCsv(csv);
  const molCol: DG.Column<string> = df.getCol('molecule');

  return await _runAutodock2(molCol, receptorPdb);
}

//name: runAutodock3
//input: dataframe name
//input: column ligandCol
export async function runAutodock3(df: DG.DataFrame, ligandCol: DG.Column): Promise<void> {

}

//name: runAutodock4
//input: file receptor { caption: 'Receptor structure file' }
//input: file gridParams { nullable: true, caption: 'Grid parameters file' }
//input: dataframe ligandTable { caption: 'Ligand table' }
//input: column ligandColumn { caption: 'Ligand molecule column', semType: Molecule }
//editor: Docking:getRunAutodockFuncEditor
export async function runAutodock4(
  receptor: DG.FileInfo, gridParams: DG.FileInfo | null, ligandTable: DG.DataFrame, ligandColumn: DG.Column
): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('AutoDock load data ...');
  try {
    // AutoDock works with .pdb or .pdbqt files, both are text formats.
    const receptorData: BiostructureData = {
      binary: false,
      data: (new TextDecoder).decode(receptor.data) ?? (await grok.dapi.files.readAsText(receptor)),
      ext: receptor.extension,
      options: {name: receptor.name,},
    };

    const data: AutoDockDataType = {
      ligandDf: ligandTable,
      ligandMolColName: ligandColumn.name,
      receptor: receptorData,
    };

    const app = new AutoDockApp();
    await app.init(data);
  } finally {
    pi.close();
  }
}

//name: getConfigFiles
//output: list<string> configFiles
export async function getConfigFiles(): Promise<string[]> {
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
  return targetsFiles.filter(file => file.isDirectory).map(file => file.name);
}

//top-menu: Chem | Autodock...
//name: Autodock
//tags: HitTriageFunction
//description: Autodock plugin UI
//input: dataframe table [Input data table]
//input: column ligands {type:categorical; semType: Molecule} [Small molecules to dock]
//input: string target {choices: Docking: getConfigFiles} [Folder with config and macromolecule]
//input: int confirmations = 10 [Number of confirmations]
export async function runAutodock5(table: DG.DataFrame, ligands: DG.Column, target: string, confirmations: number): Promise<void> {
  const isGpfFile = (file: DG.FileInfo): boolean => file.extension === 'gpf';
  const configFile = (await grok.dapi.files.list(`${TARGET_PATH}/${target}`, true)).find(isGpfFile)!;
  const receptor = (await grok.dapi.files.list(`${TARGET_PATH}/${target}`)).find((file) => file.extension === 'pdbqt')!;
  const desirableHeight = 100;
  const pi = DG.TaskBarProgressIndicator.create('AutoDock load data ...');
  try {
    // AutoDock works with .pdb or .pdbqt files, both are text formats.
    const receptorData: BiostructureData = {
      binary: false,
      data: (await grok.dapi.files.readAsText(receptor.fullPath)),
      ext: receptor.extension,
      options: {name: receptor.name,},
    };

    const data: AutoDockDataType = {
      ligandDf: table,
      ligandMolColName: ligands.name,
      receptor: receptorData,
      gpfFile: (await grok.dapi.files.readAsText(configFile.fullPath)),
      confirmationNum: confirmations,
    };

    const app = new AutoDockApp();
    const autodockResults = CACHED_DOCKING.has(data) ? CACHED_DOCKING.get(data) : await app.init(data);
    if (!autodockResults)
      return;

    const processedResults = processAutodockResults(autodockResults);
    table.join(processedResults, [], [], null, null, DG.JOIN_TYPE.INNER, true);
    /**Temporary fix (renderer works incorrectly) */
    const grid = grok.shell.addTableView(table).grid;
    
    addColorCoding(grid.columns.byName(AFFINITY_COL)!);
    grid.sort([AFFINITY_COL]);
    grid.onCellPrepare(function (gc) {
      if (gc.grid.col(POSE_COL) !== null && gc.grid.props.rowHeight !== desirableHeight)
          grid.setOptions({ 'rowHeight': desirableHeight });
    });

  } finally {
    pi.close();
  }
}

function addColorCoding(column: DG.GridColumn) {
  column.isTextColorCoded = true;
  column.column!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
  column.column!.tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.green}, ${DG.Color.red}]`;
}

function processAutodockResults(table: DG.DataFrame): DG.DataFrame {
  const affinityDescription = 'Estimated Free Energy of Binding.\
    Lower values correspond to stronger binding.';
  const processedTable = DG.DataFrame.fromColumns([table.col(POSE_COL)!, table.col(AFFINITY_COL)!]);
  processedTable.col(AFFINITY_COL)!.setTag(DG.TAGS.DESCRIPTION, affinityDescription);
  return processedTable;
}

//name: Docking
//tags: panel, chem, widgets
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
export async function autodockWidget(molecule: DG.SemanticValue): Promise<DG.Widget<any> | null> {
  let result = new DG.Widget(ui.divH([]));
  const update = () => {
    ui.remove(result.root);
    result.root.append(ui.loader());
    getAutodockSingle(molecule).then((dockingResults) => {
      ui.remove(result.root);
      result = dockingResults!;
    });
  }
  update();
  return await getAutodockSingle(molecule);
}

//name: getAutodockSingle
export async function getAutodockSingle(molecule: DG.SemanticValue): Promise<DG.Widget<any> | null> {
  const currentTable = grok.shell.tv.table;
  //@ts-ignore
  const key: AutoDockDataType = CACHED_DOCKING.K.find((key: AutoDockDataType) => key.ligandDf === currentTable);
  if (!key) {
    grok.shell.warning('Run Chem | Docking first');
    return null; 
  }

  const autodockResults = CACHED_DOCKING.get(key);
  const matchDf: DG.DataFrame = autodockResults!.rows.match(`${molecule.cell.column.name} = ${molecule.value}`).toDataFrame();
  const widgetKey = `${key.receptor}/${key.ligandMolColName}`
  const cachedWidget = CACHED_MOLSTAR.has(widgetKey);

  if (cachedWidget)
    return CACHED_MOLSTAR.get(widgetKey)!;

  const widget = new DG.Widget(ui.div([]));
  const targetViewer = await currentTable!.plot.fromType('Biostructure', {
    dataJson: BiostructureDataJson.fromData(key.receptor),
    ligandColumnName: key.ligandMolColName,
    zoom: true,
  });
  const result = ui.div();
  const map: { [_: string]: any } = {};
  for (let i = 4; i < matchDf.columns.length; ++i) {
    const columnName = matchDf.columns.names()[i];
    const columnValue = matchDf.get(columnName, 0);
    map[columnName] = columnValue;
  }
  result.appendChild(ui.tableFromMap(map));
  widget.root.append(targetViewer.root);
  widget.root.append(result);
  CACHED_MOLSTAR.set(widgetKey, widget);

  return widget;
}

// -- Demo --
