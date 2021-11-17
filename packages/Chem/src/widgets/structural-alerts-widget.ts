// This file may not be used in
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// The file is imported from a WebWorker. Don't use Datagrok imports
import { drawMoleculeToCanvas } from '../chem_common_rdkit';
import { getStructuralAlerts } from './structural-alerts';

let structuralAlertsRdKitModule: any = null;
let _webRoot: string | null = null;
let table: DG.DataFrame | null = null;
let smartsMap: Map<string, any> = new Map();

export function setStructuralAlertsRdKitModule(module: any, webRoot: string) {
  structuralAlertsRdKitModule = module;
  _webRoot = webRoot;
}

export function structuralAlertsWidget(smiles: string) {
  const alerts = getStructuralAlerts(smiles);
  const width = 200;
  const height = 100;
  const list = ui.div(alerts.map((i) => {
    const description = ui.divText(table!.get('description', i));
    const imageHost = ui.canvas(width, height);
    drawMoleculeToCanvas(0, 0, width, height, imageHost, smiles, table!.get('smarts', i));
    const host = ui.div([description, imageHost], 'd4-flex-col');
    host.style.margin = '5px';
    return host;
  }), 'd4-flex-wrap');
  if (!alerts.length) {
    list.innerText = 'No Alerts';
  }
  return new DG.Widget(list);
}