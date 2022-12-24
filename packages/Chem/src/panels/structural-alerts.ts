// This file may not be used in
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export async function checkForStructuralAlerts(col: DG.Column<string>): Promise<void> {
  if (col.length > 500) {
    grok.shell.info('Could not check structural alerts: too much data (max. 500 rows)');
    return;
  }
  const df = col.dataFrame;
  const alertsDf = await grok.data.loadTable(getRdKitWebRoot() + 'files/alert-collection.csv');
  const ruleSetCol = alertsDf.getCol('rule_set_name');
  const smartsCol = alertsDf.getCol('smarts');
  const ruleIdCol = alertsDf.getCol('rule_id');
  const alertsDfLength = alertsDf.rowCount;
  const rdkitModule = getRdKitModule();

  const dialog = ui.dialog('Structural Alerts');
  for (const ruleSet of ruleSetCol.categories)
    dialog.add(ui.boolInput(ruleSet, true));

  dialog.onOK(() => {
    const progress = DG.TaskBarProgressIndicator.create('Checking for structural alerts...');
    const ruleSetList = dialog.inputs.filter((input) => input.value).map((input) => input.caption);
    if (ruleSetList.length == 0)
      return;
    runStructuralAlertsDetection(df, ruleSetList, col, ruleSetCol, ruleIdCol, smartsMap, rdkitModule);
    progress.close();
  });

  dialog.show();

  // Caching the molecules
  const smartsMap = new Map<string, RDMol>();
  for (let i = 0; i < alertsDfLength; i++)
    smartsMap.set(ruleIdCol.get(i), rdkitModule.get_qmol(smartsCol.get(i)));
}

export function runStructuralAlertsDetection(df: DG.DataFrame, ruleSetList: string[], col: DG.Column<string>,
  ruleSetCol: DG.Column<string>, ruleIdCol: DG.Column<string>, smartsMap: Map<string, RDMol>,
  rdkitModule: RDModule): DG.DataFrame {
  ruleSetList.forEach((ruleSetName) => df.columns.addNewBool(ruleSetName));
  const skipRuleSetList: string[] = [];
  const originalDfLength = df.rowCount;
  const alertsDfLength = ruleSetCol.length;

  for (let i = 0; i < originalDfLength; i++) {
    const mol = rdkitModule.get_mol(col.get(i)!);

    for (let j = 0; j < alertsDfLength; j++) {
      const currentRuleSet = ruleSetCol.get(j)!;
      if (!ruleSetList.includes(currentRuleSet) || skipRuleSetList.includes(currentRuleSet))
        continue;

      const matches = mol.get_substruct_match(smartsMap.get(ruleIdCol.get(j)!)!);
      if (matches != '{}') {
        skipRuleSetList.push(currentRuleSet);
        df.set(currentRuleSet, i, true);
      }
    }
    skipRuleSetList.length = 0;
    mol.delete();
  }

  return df;
}
