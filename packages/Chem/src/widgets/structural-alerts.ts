import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import { drawMoleculeToCanvas } from '../chem_common';

let structuralAlertsRdKitModule: any = null;
let table: DG.DataFrame | null = null;
let smartsMap: Map<string, any> = new Map();

export async function setStructuralAlertsRdKitModule(module: any) {
  structuralAlertsRdKitModule = module;
}

export async function loadAlertsCollection() {
  const path = _package.webRoot + 'data-samples/alert_collection.csv';
  table = await grok.data.loadTable(path);
  for (let i = 0; i < table.rowCount; i++) {
    const currentSmarts = table.get('smarts', i);
    smartsMap.set(currentSmarts, structuralAlertsRdKitModule.get_qmol(currentSmarts));
  }
}

export function structuralAlertsWidget(smiles: string) {
  const alerts = getStructuralAlerts(smiles);

  const width = 200;
  const height = 100;

  const list = ui.div(alerts.map((i) => {
    const description = ui.divText(table!.get('description', i));
    const imageHost = ui.canvas(width, height);
    drawMoleculeToCanvas(0, 0, width, height, imageHost, smiles, table!.get('smarts', i));
    return ui.div([description, imageHost], 'd4-flex-col');
  }), 'd4-flex-wrap');
  if (!alerts.length) {
    list.innerText = 'No Alerts';
  }
  return new DG.Widget(list);
}

function getStructuralAlerts(smiles: string) {
  const alerts: number[] = [];
  const mol = structuralAlertsRdKitModule.get_mol(smiles);
  for (let i = 0; i < table!.rowCount; i++) {
    const subMol = smartsMap.get(table!.get('smarts', i));
    const matches = mol.get_substruct_matches(subMol);
    if (matches !== '{}') {
      alerts.push(i);
    }
  }
  return alerts;
}
