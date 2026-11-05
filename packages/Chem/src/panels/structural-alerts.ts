import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {alertsByEngine} from '../engine/worker-jobs';

export type RuleId = 'PAINS' | 'BMS' | 'SureChEMBL' | 'MLSMR' | 'Dundee' | 'Inpharmatica' | 'LINT' | 'Glaxo';
export type RuleSet = {[rule in RuleId]: boolean};
export const STRUCT_ALERTS_RULES_NAMES = ['PAINS', 'BMS', 'SureChEMBL', 'MLSMR', 'Dundee', 'Inpharmatica', 'LINT', 'Glaxo'];

export async function runStructuralAlertsDetection(moleculeCol: DG.Column<string>, ruleSet: RuleSet,
  alertsDf: DG.DataFrame,
  matcher: typeof alertsByEngine): Promise<DG.DataFrame> {
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

  const result = await matcher.getStructuralAlerts(smarts, moleculeCol.toList());

  // Build the result dataframe
  const resultDf = DG.DataFrame.create(moleculeCol.length);
  for (const [ruleName, boolArray] of result)
    resultDf.columns.addNewBool(resultDf.columns.getUnusedName(ruleName)).init((i) => boolArray[i]);

  return resultDf;
}

export async function getStructuralAlertsByRules(molecules: DG.Column, ruleSet: RuleSet): Promise<DG.DataFrame | undefined> {
  const alertsDf = await grok.data.loadTable(chemCommonRdKit.getRdKitWebRoot() + 'files/alert-collection.csv');
  const progress = DG.TaskBarProgressIndicator.create('Detecting structural alerts...');
  try {
    // by the "Chemical engine" setting: the RDKit pool, or in Crux mode the Crux pool (RDKit's behind it)
    const resultDf = await runStructuralAlertsDetection(molecules, ruleSet, alertsDf, alertsByEngine);
    return resultDf;
  } catch (e) {
    grok.shell.error('Structural alerts detection failed');
    grok.log.error(`Structural alerts detection failed: ${e}`);
  } finally {
    progress.close();
  }
}