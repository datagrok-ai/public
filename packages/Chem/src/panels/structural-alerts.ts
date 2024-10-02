import * as DG from 'datagrok-api/dg';
import {RdKitService} from '../rdkit-service/rdkit-service';

export type RuleId = 'PAINS' | 'BMS' | 'SureChEMBL' | 'MLSMR' | 'Dundee' | 'Inpharmatica' | 'LINT' | 'Glaxo';
export type RuleSet = {[rule in RuleId]: boolean};

export async function runStructuralAlertsDetection(moleculeCol: DG.Column<string>, ruleSet: RuleSet,
  alertsDf: DG.DataFrame,
  rdkitService: RdKitService): Promise<DG.DataFrame> {
  const ruleSetCol: DG.Column<string> = alertsDf.getCol('rule_set_name');
  const ruleSetColCategories = ruleSetCol.categories;
  const ruleSetColData = ruleSetCol.getRawData();
  const smartsCol: DG.Column<string> = alertsDf.getCol('smarts');
  const smartsColCategories = smartsCol.categories;
  const smartsColData = smartsCol.getRawData();

  const smarts: {[rule in RuleId]?: string[]} = {};
  for (let i = 0; i < smartsColData.length; i++) {
    const ruleName = ruleSetColCategories[ruleSetColData[i]] as RuleId;
    if (ruleSet[ruleName]) {
      smarts[ruleName] ??= [];
      smarts[ruleName]!.push(smartsColCategories[smartsColData[i]]);
    }
  }

  const result = await rdkitService.getStructuralAlerts(smarts, moleculeCol.toList());

  // Build the result dataframe
  const resultDf = DG.DataFrame.create(moleculeCol.length);
  for (const [ruleName, boolArray] of result)
    resultDf.columns.addNewBool(resultDf.columns.getUnusedName(ruleName)).init((i) => boolArray[i]);

  return resultDf;
}
