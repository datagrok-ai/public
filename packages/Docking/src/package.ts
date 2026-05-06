/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/docking.css';

import '@datagrok-libraries/bio/src/types/ngl'; // To enable import from the NGL module declared in bio lib
import {IAutoDockService} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {AutoDockApp, AutoDockDataType} from './apps/auto-dock-app';
import {_runAutodock, AutoDockService, _runAutodock2, ensureNoDockingError} from './utils/auto-dock-service';
import {TARGET_PATH, BINDING_ENERGY_COL, POSE_COL, BINDING_ENERGY_COL_UNUSED, POSE_COL_UNUSED, ERROR_COL_NAME, ERROR_MESSAGE, AUTODOCK_PROPERTY_DESCRIPTIONS} from './utils/constants';
import { _demoDocking } from './demo/demo';
import { DockingViewApp } from './demo/docking-app';
import { addColorCoding, buildComparisonTable, formatColumns, getRemarksFromPdb, getRemarksFromPdbs, getReceptorData, processAutodockResults, prop } from './utils/utils';
export * from './package.g';
export const _package = new DG.Package();

const PROLIF_SKIP_RESNAMES = new Set([
  'HOH', 'WAT', 'H2O', 'D2O', 'DOD',
  'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU', 'NI',
  'SO4', 'PO4', 'NO3', 'ACT', 'CO3',
  'GOL', 'EDO', 'PEG', 'PG4', 'DMS', 'TRS', 'IMD', 'BME',
]);

// Returns each unique non-water HETATM ligand instance as "RESNAME CHAIN RESID".
function detectNonWaterHetatmInstances(pdbText: string): string[] {
  if (!pdbText) return [];
  const seen = new Map<string, number>();
  for (const line of pdbText.split('\n')) {
    if (line.startsWith('HETATM') && line.length >= 26) {
      const rn = line.slice(17, 20).trim();
      const chain = (line[21] || '').trim() || 'A';
      const resid = line.slice(22, 26).trim();
      if (rn && !PROLIF_SKIP_RESNAMES.has(rn)) {
        const key = `${rn} ${chain} ${resid}`;
        seen.set(key, (seen.get(key) || 0) + 1);
      }
    }
  }
  return Array.from(seen.entries()).sort((a, b) => b[1] - a[1]).map(([k]) => k);
}

// See BiostructureViewer/src/package.ts for the full comment.
function injectCompactCss(html: string): string {
  const css = `
<style>
  * { box-sizing: border-box; }
  html { margin: 0 !important; padding: 0 !important; }
  body { margin: 0 !important; padding: 0 !important; background: white; }
  body > * { margin: 0 !important; }
  div[id^="mynetwork"], .vis-network, #mynetwork {
    height: 520px !important;
    width: 100% !important;
    padding: 0 !important;
    margin: 0 !important;
  }
  [class*="legend"], div[class*="-legend"] {
    padding: 4px 8px !important;
    margin: 0 !important;
  }
</style>
<script>
  document.addEventListener('DOMContentLoaded', function() {
    var sent = false;
    var send = function() {
      if (sent) return;
      sent = true;
      var h = Math.max(
        document.documentElement.scrollHeight || 0,
        document.body ? (document.body.scrollHeight || 0) : 0
      );
      parent.postMessage({ type: 'prolif-ready', height: h }, '*');
    };
    var tryHook = function() {
      if (typeof network !== 'undefined' && network && network.on) {
        network.on('stabilizationIterationsDone', function() {
          try { network.fit({animation: false}); } catch (e) {}
          setTimeout(send, 60);
        });
        setTimeout(send, 250);
      } else {
        setTimeout(tryHook, 50);
      }
    };
    tryHook();
  });
</script>`;
  if (html.includes('<head>'))
    return html.replace('<head>', '<head>' + css);
  if (html.includes('<body>'))
    return html.replace('<body>', '<body>' + css);
  return css + html;
}

function makeProlifWidget(params: {
  protein: string;
  ligand?: string;
  ligand_resname?: string;
}): DG.Widget {
  const host = ui.div([], 'd4-empty-parent');

  const detectionSource = (params.ligand && params.ligand.trim()) || params.protein;
  const ligands = detectNonWaterHetatmInstances(detectionSource);

  const compute = (resname: string, body: HTMLElement) => {
    ui.empty(body);
    const loader = ui.loader();
    body.append(loader);
    (async () => {
      try {
        const html = await grok.functions.call(
          'BiostructureViewer:ProteinLigandInteractionDiagram',
          {
            protein: params.protein,
            ligand: params.ligand ?? '',
            ligand_resname: resname,
          },
        ) as string;

        const iframe = document.createElement('iframe');
        iframe.srcdoc = injectCompactCss(html);
        iframe.style.cssText =
          'width:100%; height:600px; border:0; display:block; opacity:0; transition:opacity 0.2s;';
        iframe.setAttribute('sandbox', 'allow-scripts');

        const reveal = (h?: number) => {
          if (typeof h === 'number' && h > 0)
            iframe.style.height = `${Math.max(300, Math.min(h + 4, 900))}px`;
          if (loader.isConnected) loader.remove();
          iframe.style.opacity = '1';
        };
        const onMsg = (e: MessageEvent) => {
          if (e.source !== iframe.contentWindow) return;
          const data = e.data;
          if (data && typeof data === 'object' && data.type === 'prolif-ready') {
            window.removeEventListener('message', onMsg);
            reveal(typeof data.height === 'number' ? data.height : undefined);
          } else if (data === 'prolif-ready') {
            window.removeEventListener('message', onMsg);
            reveal();
          }
        };
        window.addEventListener('message', onMsg);
        setTimeout(() => { window.removeEventListener('message', onMsg); reveal(); }, 8000);

        body.append(iframe);
      } catch (err) {
        ui.empty(body);
        body.append(ui.divText(
          `Could not compute interactions: ${err instanceof Error ? err.message : String(err)}`,
        ));
      }
    })();
  };

  if (ligands.length === 0 && !params.ligand_resname) {
    host.append(ui.divText('No non-water HETATM ligand found in this structure.'));
  } else if (params.ligand_resname || ligands.length === 1) {
    compute(params.ligand_resname || ligands[0], host);
  } else {
    const body = ui.div();
    const picker = ui.input.choice('Ligand', {
      value: null,
      items: ligands,
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v) compute(v, body);
      },
    });
    host.append(ui.divV([
      ui.divText(
        `${ligands.length} ligands found. Select one to compute interactions:`,
        {style: {marginBottom: '6px', color: 'var(--grey-5)'}},
      ),
      picker.root,
      body,
    ]));
  }

  return new DG.Widget(host);
}


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
    meta: {vectorFunc: 'true'},
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
    // `processAutodockResults` returns an empty DataFrame when AutoDock's
    // output is missing required columns — skip downstream wiring.
    if (processedResults.columns.length === 0)
      return processedResults;
    await grok.data.detectSemanticTypes(processedResults);

    // NOTE: the SMILES↔pose atom-picker link is NOT written here. Users
    // activate it explicitly via the "Link SMILES column" checkbox in the
    // AutoDock context panel (`mol3dAtomPickerLinkWidget`). Auto-linking
    // was rejected to avoid ambiguous pairings when the table has more
    // than one SMILES column and to keep the interactive-highlight opt-in
    // rather than implicit.
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
    @grok.decorators.param({options: {semType: 'Molecule', description: '\'Small molecules to dock\''}}) ligands: DG.Column,
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

  @grok.decorators.panel({
    name: 'Protein-Ligand Interactions',
    condition: 'Docking:isApplicableAutodock(molecule)',
  })
  static async dockingInteractionsWidget(
    @grok.decorators.param({options: {semType: 'Molecule3D'}}) molecule: DG.SemanticValue
  ): Promise<DG.Widget> {
    const pose = molecule.value as string;
    const receptorData = await getReceptorData(pose);
    const receptor = typeof receptorData.data === 'string'
      ? receptorData.data
      : new TextDecoder().decode(receptorData.data as Uint8Array);
    return makeProlifWidget({protein: receptor, ligand: pose});
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

    const autodockResults: DG.DataFrame = getRemarksFromPdbs(molecule);
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
    const buildPropertyMap = () => {
      const map: { [_: string]: any } = {};
      for (let i = 0; i < autodockResults!.columns.length; ++i) {
        const columnName = autodockResults!.columns.names()[i];
        const propertyCol = autodockResults!.col(columnName);
        map[columnName] = prop(molecule, propertyCol!, result, AUTODOCK_PROPERTY_DESCRIPTIONS);
      }
      return ui.tableFromMap(map);
    };
    result.appendChild(buildPropertyMap());
    widget.root.append(result);

    const currentRowIdx = molecule.cell.rowIndex;
    const ligandCol = molecule.cell.column;
    const currentValues = getRemarksFromPdb(ligandCol.get(currentRowIdx));

    const sub = DG.debounce(currentTable.onMouseOverRowChanged, 50).subscribe(() => {
      const hovIdx = currentTable.mouseOverRowIdx;
      ui.empty(result);
      if (hovIdx >= 0 && hovIdx !== currentRowIdx) {
        const hovPdb = ligandCol.get(hovIdx);
        if (hovPdb && hovPdb.includes(BINDING_ENERGY_COL)) {
          result.appendChild(buildComparisonTable(currentValues, getRemarksFromPdb(hovPdb)));
          return;
        }
      }
      result.appendChild(buildPropertyMap());
    });
    widget.subs.push(sub);

    // Append the atom-picker "Link SMILES column" checkbox at the bottom
    // of the AutoDock panel. Rendered via BiostructureViewer's exposed
    // function rather than as a separate cell panel, so the UI stays as
    // a single "AutoDock" section in the context panel instead of two.
    // Fails silently if BSV isn't installed / the function isn't found.
    try {
      const linkWidget = await grok.functions.call(
        'BiostructureViewer:mol3dAtomPickerLinkWidget',
        {mol3DCol: molecule.cell.column},
      ) as DG.Widget | null;
      if (linkWidget?.root)
        widget.root.append(linkWidget.root);
    } catch (err: unknown) {
      _package.logger.debug(
        `mol3dAtomPickerLinkWidget unavailable: ${err instanceof Error ? err.message : String(err)}`);
    }

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
