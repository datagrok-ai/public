/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '@datagrok-libraries/bio/src/types/ngl'; // To enable import from the NGL module declared in bio lib
import {IAutoDockService} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {AutoDockApp, AutoDockDataType} from './apps/auto-dock-app';
import {_runAutodock, AutoDockService, _runAutodock2, ensureNoDockingError} from './utils/auto-dock-service';
import {TARGET_PATH, BINDING_ENERGY_COL, POSE_COL, BINDING_ENERGY_COL_UNUSED, POSE_COL_UNUSED, ERROR_COL_NAME, ERROR_MESSAGE, AUTODOCK_PROPERTY_DESCRIPTIONS} from './utils/constants';
import { _demoDocking } from './demo/demo';
import { DockingViewApp } from './demo/docking-app';
import { addColorCoding, formatColumns, getFromPdbs, getReceptorData, processAutodockResults, prop } from './utils/utils';
export * from './package.g';
export const _package = new DG.Package();


export class PackageFunctions{
  @grok.decorators.func()
  static info() {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func({
    'outputs': [
      {
        'name': 'result',
        'type': 'object',
      }
    ]
  })
  static async getAutoDockService(): Promise<IAutoDockService> {
    const resSvc: IAutoDockService = await AutoDockService.getSvc();
    return resSvc;
  }

  @grok.decorators.func()
  static async autoDockApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('AutoDock app...');
    try {
      const app = new AutoDockApp('autoDockApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async getConfigFiles(): Promise<string[]> {
    const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
    const directoriesWithGpf = await Promise.all(
      targetsFiles.filter(file => file.isDirectory).map(async dir => {
        const filesInDir = await grok.dapi.files.list(dir.fullPath, true);
        return filesInDir.some(file => file.path.endsWith('.gpf')) ? dir.name : null;
      })
    );
    return directoriesWithGpf.filter((dir): dir is string => Boolean(dir));
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async dockLigandCached(
    jsonForm: string,
    containerId: string): Promise<string> {

    const params: RequestInit = {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: jsonForm,
    };

    const path = `/autodock/dock_ligand`;
    const dockingResult = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, path, params);
    const dockingResultText = await dockingResult.json();
    ensureNoDockingError(dockingResultText);
    return dockingResultText;
  }

  @grok.decorators.func({
    name: 'getAutodockResults',
    outputs: [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async getAutodockResults(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) ligands: DG.Column,
    target: string,
    poses: number
  ): Promise<DG.DataFrame> {
    const data = await prepareAutoDockData(target, table, ligands.name, poses);
    if (!data)
      return DG.DataFrame.create();

    const app = new AutoDockApp();
    const autodockResults = await app.init(data);
    if (!autodockResults)
      return DG.DataFrame.create();

    formatColumns(autodockResults);
    const processedResults = processAutodockResults(autodockResults, table);
    await grok.data.detectSemanticTypes(processedResults);
    return processedResults;
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Docking | AutoDock...',
    'name': 'AutoDock',
    'description': 'Autodock plugin UI',
    meta: {role: 'hitTriageFunction'},
  })
  static async runAutodock(
    @grok.decorators.param({options: {description: '\'Input data table\''}}) table: DG.DataFrame,
    @grok.decorators.param({options: {type: 'categorical', semType: 'Molecule', description: '\'Small molecules to dock\''}}) ligands: DG.Column,
    @grok.decorators.param({options: {choices: 'Docking:getConfigFiles', description: '\'Folder with config and macromolecule\''}}) target: string,
    @grok.decorators.param({type: 'int', options: {initialValue: '10', description: '\'Number of output conformations for each small molecule\''}}) poses: number): Promise<void> {

    const desirableHeight = 100;
    const desirableWidth = 100;
    const pi = DG.TaskBarProgressIndicator.create('AutoDock load data ...');
    try {
      await grok.functions.call('Docking:getAutodockResults', { table, ligands, target, poses });

      if (!BINDING_ENERGY_COL_UNUSED || !POSE_COL_UNUSED)
        return;

      const {grid} = grok.shell.getTableView(table.name);

      addColorCoding(grid.columns.byName(BINDING_ENERGY_COL_UNUSED)!);
      grid.sort([BINDING_ENERGY_COL_UNUSED]);

      await grok.data.detectSemanticTypes(table);

      grid.onCellRender.subscribe((args: any) => {
        grid.setOptions({ 'rowHeight': desirableHeight });
        grid.col(POSE_COL_UNUSED)!.width = desirableWidth;
        grid.col(BINDING_ENERGY_COL_UNUSED)!.width = desirableWidth + 50;

        const { cell, g, bounds } = args;
        const value = cell.cell.value;
        const isPoseCell = cell.isTableCell && cell.cell.column.name === POSE_COL_UNUSED;
        const isErrorValue = typeof value === 'string' && value.toLowerCase().includes(ERROR_COL_NAME);

        if (isPoseCell && isErrorValue) {
          g.fillStyle = 'black';
          g.fillText(ERROR_MESSAGE, bounds.x + bounds.width / 2, bounds.y + bounds.height / 2);
          args.preventDefault();
        }
      });
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static isApplicableAutodock(molecule: string): boolean {
    return molecule.includes('binding energy');
  }

  @grok.decorators.panel({
    'name': 'AutoDock',
    'condition': 'Docking:isApplicableAutodock(molecule)',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async autodockWidget(
    @grok.decorators.param({'options':{'semType':'Molecule3D'}}) molecule: DG.SemanticValue): Promise<DG.Widget<any> | null> {
    return await PackageFunctions.getAutodockSingle(molecule);
  }

  @grok.decorators.func()
  static async getAutodockSingle(
    molecule: DG.SemanticValue, showProperties: boolean = true,
    table?: DG.DataFrame): Promise<DG.Widget<any> | null> {

    const value = molecule.value;
    if (value.toLowerCase().includes(ERROR_COL_NAME))
      return new DG.Widget(ui.divText(value));

    const tableView = grok.shell.tv;
    const currentTable = table ?? tableView.dataFrame;

    const addedToPdb = value.includes(BINDING_ENERGY_COL);
    if (!addedToPdb)
      return new DG.Widget(ui.divText('Docking has not been run'));

    const autodockResults: DG.DataFrame = getFromPdbs(molecule);
    const widget = new DG.Widget(ui.div([]));

    if (table)
      currentTable.currentRowIdx = 0;

    const receptorData = await getReceptorData(value);
    const targetViewer = await currentTable.plot.fromType('Biostructure', {
      dataJson: BiostructureDataJson.fromData(receptorData),
      ligandColumnName: molecule.cell.column.name,
      zoom: true,
    });
    targetViewer.root.classList.add('bsv-container-info-panel');
    widget.root.append(targetViewer.root);
    if (!showProperties) return widget;

    const result = ui.div();
    const map: { [_: string]: any } = {};
    for (let i = 0; i < autodockResults!.columns.length; ++i) {
      const columnName = autodockResults!.columns.names()[i];
      const propertyCol = autodockResults!.col(columnName);
      map[columnName] = prop(molecule, propertyCol!, result, AUTODOCK_PROPERTY_DESCRIPTIONS);
    }
    result.appendChild(ui.tableFromMap(map));
    widget.root.append(result);

    return widget;
  }

  @grok.decorators.func({
    'meta': {
      'demoPath': 'Bioinformatics | Docking'
    },
    'name': 'Demo Docking',
    'description': 'Small molecule docking to a macromolecule with pose visualization'
  })
  static async demoDocking(): Promise<void> {
    await _demoDocking();
  }

  @grok.decorators.panel({
    'name': 'Biology | AutoDock',
    meta: {role: 'widgets'},
  })
  static async autodockPanel(
    @grok.decorators.param({'options':{'semType':'Molecule'}}) smiles: DG.SemanticValue): Promise<DG.Widget> {

    const items = await PackageFunctions.getConfigFiles();
    const target = ui.input.choice('Target', {value: items[0], items: items})
    const poses = ui.input.int('Poses', {value: 10});

    const resultsContainer = ui.div();
    const button = ui.button('Run', async () => {
      resultsContainer.innerHTML = '';

      const loader = ui.loader();
      resultsContainer.appendChild(loader);

      const widget = await runDocking(smiles, target.value!, poses.value!);
      resultsContainer.removeChild(loader);
      if (widget) {
        resultsContainer.appendChild(widget.root);
      }
    });

    const form = ui.form([target, poses]);
    const panels = ui.divV([form, button, resultsContainer]);

    return DG.Widget.fromRoot(panels);
  }

  @grok.decorators.app({
    'meta': {
      'icon': 'images/docking-icon.png',
      'browsePath': 'Bio'
    },
    'name': 'Docking'
  })
  static async dockingView(
    @grok.decorators.param({ 'options':{'meta.url':true,'optional':true}})  path?: string): Promise<DG.ViewBase | null> {
    const parent = grok.functions.getCurrentCall();
    const app = new DockingViewApp(parent);
    await app.init();
    return app.tableView!;
  }
}

export async function runDocking(
  smiles: DG.SemanticValue,
  target: string,
  poses: number
): Promise<DG.Widget | null> {

  const ligandColumnName = 'ligand';
  const column = DG.Column.fromStrings(ligandColumnName, [smiles.value]);
  const table = DG.DataFrame.fromColumns([column]);
  column.semType = DG.SEMTYPE.MOLECULE;
  await grok.data.detectSemanticTypes(table);

  const data = await prepareAutoDockData(target, table, ligandColumnName, poses);
  if (!data)
    return null;

  const app = new AutoDockApp();
  const autodockResults = await app.init(data);

  if (autodockResults) {
    const pose = autodockResults.cell(0, POSE_COL);
    return await PackageFunctions.getAutodockSingle(DG.SemanticValue.fromTableCell(pose!), false, autodockResults);
  }

  return null;
}

export async function prepareAutoDockData(
  target: string,
  table: DG.DataFrame,
  ligandColumn: string,
  poses: number
): Promise<AutoDockDataType | null> {
  const isGpfFile = (file: DG.FileInfo): boolean => file.extension === 'gpf';
  const configFile = (await grok.dapi.files.list(`${TARGET_PATH}/${target}`, true)).find(isGpfFile)!;
  const receptor = (await grok.dapi.files.list(`${TARGET_PATH}/${target}`))
    .find(file => ['pdbqt', 'pdb'].includes(file.extension)) || null;

  if (!configFile || !receptor) {
    grok.shell.warning('Missing .gpf or .pdbqt/.pdb file in the target folder.');
    return null;
  }

  const receptorData: BiostructureData = {
    binary: false,
    data: await grok.dapi.files.readAsText(receptor.fullPath),
    ext: receptor.extension,
    options: { name: receptor.name },
  };

  return {
    ligandDf: table,
    ligandMolColName: ligandColumn,
    receptor: receptorData,
    gpfFile: await grok.dapi.files.readAsText(configFile.fullPath),
    posesNum: poses,
    ligandDfString: table.columns.byName(ligandColumn).toString(),
  };
}
