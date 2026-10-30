import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {getRdKitModule, drawMoleculeToCanvas, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {HIGHLIGHT_BY_SCAFFOLD_TAG} from '../constants';
import {IColoredScaffold} from '../rendering/rdkit-cell-renderer';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {byEngineModule} from '../engine/chem-engine';


let alertsDf: DG.DataFrame | null = null;
/** The alerts' query molecules, one map per module (RDKit's, and Crux's in Crux mode: `byEngineModule`). */
const _smartsMaps: Map<RDModule, Map<string, RDMol | null>> = new Map();
let rdKitModule: RDModule | null = null;
const NO_HIGHLIGHT = 0;

/** The rows of the alerts `molecule` matches: on Crux in Crux mode, on RDKit for a molecule Crux cannot read. */
export async function getStructuralAlerts(molecule: string): Promise<number[]> {
  if (alertsDf == null)
    alertsDf = await grok.data.loadTable(getRdKitWebRoot() + 'files/alert-collection.csv');
  if (!grok.chem.isMolBlock(molecule) && molecule?.length > MAX_SMILES_LENGTH)
    throw new Error(`SMILES string longer than ${MAX_SMILES_LENGTH} characters not supported`);
  return byEngineModule('structuralAlerts', getRdKitModule(), (module, crux) => {
    const mol = module.get_mol(molecule);
    if (!mol && crux)
      return null;
    try {
      const queries = queryMols(module);
      if (crux && [...queries.values()].includes(null))
        throw new Error('not supported by Crux: an alert\'s SMARTS');
      const smartsCol = alertsDf!.getCol('smarts');
      const alerts: number[] = [];
      for (let i = 0; i < smartsCol.length; i++) {
        const subMol = queries.get(smartsCol.get(i));
        if (subMol && mol!.get_substruct_matches(subMol) !== '{}')
          alerts.push(i);
      }
      return alerts;
    } finally {
      mol?.delete();
    }
  })!;
}

/** The alerts' query molecules on `module`, read once (a SMARTS it cannot read: none). */
function queryMols(module: RDModule): Map<string, RDMol | null> {
  let map = _smartsMaps.get(module);
  if (!map) {
    map = new Map();
    for (const smarts of alertsDf!.getCol('smarts').toList()) {
      try {
        map.set(smarts, module.get_qmol(smarts));
      } catch {
        map.set(smarts, null);
      }
    }
    _smartsMaps.set(module, map);
  }
  return map;
}

export async function structuralAlertsWidget(molecule: string): Promise<DG.Widget> {
  const colors = [NO_HIGHLIGHT].concat(DG.Color.categoricalPalette.slice(0, 10));
  rdKitModule ??= getRdKitModule();
  let alerts: any[] = [];
  try {
    alerts = await getStructuralAlerts(molecule);
  } catch (e) {
    console.warn(e);
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  if (alerts.length === 0)
    return new DG.Widget(ui.divText('No alerts'));

  const width = 200;
  const height = 100;
  const descriptionCol = alertsDf!.getCol('description');
  const smartsCol = alertsDf!.getCol('smarts');
  const calcForWholeButton = ui.button('Calculate for whole dataset', async () => {
    const alertsFunc = DG.Func.find({package: 'Chem', name: 'structuralAlertsTopMenu'})[0];
    const alertsFuncCall = alertsFunc.prepare();
    ui.dialog('Structural alerts')
      .add(await alertsFuncCall.getEditor())
      .onOK(() => {
        const args: any = {};
        Object.entries(alertsFuncCall.inputs).forEach(([key, val]) => args[key] = val);
        alertsFunc.apply(args);
      })
      .show({center: true});
  });
  calcForWholeButton.style.justifyContent = 'flex-start';
  const list = ui.div(alerts.map((i) => {
    const description = ui.divText(descriptionCol.get(i));
    const imageHost = ui.canvas(width, height);
    //in case molecule is smiles setting correct coordinates to save molecule orientation
    if (!DG.chem.isMolBlock(molecule))
      molecule = _convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdKitModule!);
    drawMoleculeToCanvas(0, 0, width, height, imageHost, molecule, smartsCol.get(i));
    const moreBtn = ui.iconFA(
      'ellipsis-v',
      (e: MouseEvent) => {
        e.stopImmediatePropagation();
        const menu = DG.Menu.popup();
        menu.group('Highlight fragment')
          .items(colors.map((color) => ui.tools.click(getColoredDiv(color), () => {
            const substr = smartsCol.get(i);
            const col = grok.shell.tv.dataFrame.currentCol;
            const array: IColoredScaffold[] = col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG) ?
              JSON.parse(col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG)) : [];
            const substrIdx = array.findIndex((it) => it.molecule === substr);
            if (substrIdx !== -1) {
              if (color !== NO_HIGHLIGHT)
                array[substrIdx].color = DG.Color.toHtml(color)!;
              else
                array.splice(substrIdx, 1);
            } else {
              if (color !== NO_HIGHLIGHT)
                array.push({molecule: smartsCol.get(i), color: DG.Color.toHtml(color)!});
            }
            col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(array));
            grok.shell.tv.dataFrame?.fireValuesChanged();
          })), () => { });
        /*           menu.item('Filter by alert', () => {
            filterByAlert(grok.shell.tv.dataFrame.currentCol, smartsCol.get(i));
          }) */
        menu.show();
      },
      'More',
    );
    $(moreBtn).addClass('chem-mol-view-icon pep-more-icon');
    const host = ui.div([description,
      ui.divV([moreBtn, imageHost], 'chem-mol-box struct-alerts-mol-box')], 'd4-flex-col');
    host.style.margin = '5px';
    return host;
  }), {classes: 'd4-flex-wrap chem-search-panel-wrapper'});

  return new DG.Widget(ui.divV([calcForWholeButton, ui.box(list, {style: {height: '100%'}})]));
}

function getColoredDiv(color: number): HTMLDivElement {
  return color === NO_HIGHLIGHT ?
    ui.div('None', {style: {width: '100%', minHeight: '20px', marginLeft: '2px'}}) :
    ui.div('', {style:
      {width: '100%', minHeight: '20px', marginRight: '6px', backgroundColor: DG.Color.toHtml(color)}});
}

function filterByAlert(molCol: DG.Column, alert: string): void {
  const tv = grok.shell.tv;
  //@ts-ignore
  const filterState =
    tv.getFiltersGroup({createDefaultFilters: false}).getStates(molCol.name, 'Chem:scaffoldTreeFilter');
  let newScaffoldTree = [{scaffold: alert}];
  if (filterState.length) {
    // @ts-ignore
    newScaffoldTree = JSON.parse(filterState[0].savedTree);
    if (newScaffoldTree.findIndex((it: any) => it.scaffold === alert) !== -1)
      grok.shell.warning(`Structure already added to filter`);
    else
      newScaffoldTree = newScaffoldTree.concat([{scaffold: alert}]);
  }
  tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
    type: 'Chem:scaffoldTreeFilter',
    column: molCol.name,
    columnName: molCol.name,
    savedTree: JSON.stringify(newScaffoldTree),
  });
}
