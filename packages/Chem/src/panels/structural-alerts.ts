// This file may not be used in
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {RDMol} from '../rdkit-api';

export async function checkForStructuralAlerts(col: DG.Column<string>): Promise<void> {
  const df = col.dataFrame;
  const originalDfLength = df.rowCount;
  const alertsDf = await grok.data.loadTable(getRdKitWebRoot() + 'files/alert-collection.csv');
  const ruleSetCol = alertsDf.getCol('rule_set_name');
  const smartsCol = alertsDf.getCol('smarts');
  const ruleIdCol = alertsDf.getCol('rule_id');
  const alertsDfLength = alertsDf.rowCount;
  const rdkitModule = getRdKitModule();
  
  // Caching the molecules
  const smartsMap = new Map<string, RDMol>();
  for (let i = 0; i < alertsDfLength; i++)
    smartsMap.set(ruleIdCol.get(i), rdkitModule.get_qmol(smartsCol.get(i)));

  const dialog = ui.dialog('Structural Alerts');
  for (const ruleSet of ruleSetCol.categories)
    dialog.add(ui.boolInput(ruleSet, true));

  dialog.onOK(() => {
    const progress = DG.TaskBarProgressIndicator.create('Checking for structural alerts...');
    const ruleSetList = dialog.inputs.filter((input) => input.value).map((input) => input.caption);
    if (ruleSetList.length == 0)
      return;

    ruleSetList.forEach((ruleSetName) => df.columns.addNewBool(ruleSetName))
    const skipRuleSetList: string[] = [];

    for (let i = 0; i < originalDfLength; i++) {
      const mol = rdkitModule.get_mol(col.get(i)!);

      for (let j = 0; j < alertsDfLength; j++) {
        const currentRuleSet = ruleSetCol.get(j);
        if (!ruleSetList.includes(currentRuleSet) || skipRuleSetList.includes(currentRuleSet))
          continue;

        const matches = mol.get_substruct_match(smartsMap.get(ruleIdCol.get(j))!);
        if (matches != '{}') {
          skipRuleSetList.push(currentRuleSet);
          df.set(currentRuleSet, i, true);
        }
      }
      skipRuleSetList.length = 0;
    }

    progress.close();
  });

  dialog.show();
}