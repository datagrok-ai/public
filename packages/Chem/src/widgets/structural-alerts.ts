import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import { drawMoleculeToCanvas } from '../chem_common';

let structuralAlertsRdKitModule: any = null;

export async function setStructuralAlertsRdKitModule(module: any) {
  structuralAlertsRdKitModule = module;
}

export async function structuralAlertsWidget(smiles: string) {
  const path = _package.webRoot + 'data-samples/alert_collection.csv';
  const table = await grok.data.loadTable(path);
  const alerts = await getStructuralAlerts(table, smiles);

  const ruleRowMap: {[index: number]: DG.Row} = {};
  for (const row of table.rows) {
    //@ts-ignore
    ruleRowMap[row['rule_id']] = row;
  }

  const width = 200;
  const height = 100;

  const list = ui.div(alerts.map((i) => {
    //@ts-ignore
    const description = ui.divText(ruleRowMap[i]['description']);
    // const alertImage = drawMolHighlights(smiles, ruleRowMap[i]['smarts']);
    const imageHost = ui.canvas(width, height);
    //@ts-ignore
    drawMoleculeToCanvas(0, 0, width, height, imageHost, smiles, ruleRowMap[i]['smarts']);
    return ui.div([description, imageHost], 'd4-flex-col');
  }), 'd4-flex-wrap');
  if (!alerts.length) {
    list.innerText = 'No Alerts';
  }
  return new DG.Widget(list);
}

async function getStructuralAlerts(table: DG.DataFrame, smiles: string) {
  const alerts: number[] = [];
  const mol = structuralAlertsRdKitModule.get_mol(smiles);
  for (let i = 0; i < table.rowCount; i++) {
    const subMol = structuralAlertsRdKitModule.get_qmol(table.get('smarts', i));
    const matches = mol.get_substruct_matches(subMol);
    if (matches !== '{}') {
      alerts.push(table.get('rule_id', i));
    }
  }
  return alerts;
}
