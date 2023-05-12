import * as DG from 'datagrok-api/dg';

import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export type RuleId = 'PAINS' | 'BMS' | 'SureChEMBL' | 'MLSMR' | 'Dandee' | 'Inpharmatica' | 'LINT' | 'Glaxo';
export type RuleSet = {[rule in RuleId]: boolean};

export function runStructuralAlertsDetection(moleculeCol: DG.Column<string>, ruleSet: RuleSet, alertsDf: DG.DataFrame,
  rdkitModule: RDModule): DG.DataFrame {
  const moleculeColCategories = moleculeCol.categories;
  const moleculeColRawData = moleculeCol.getRawData();
  
  // Get the necessary columns from the alerts dataframe
  const alertsDfLength = alertsDf.rowCount;
  const ruleSetCol: DG.Column<string> = alertsDf.getCol('rule_set_name');
  const ruleSetColCategories = ruleSetCol.categories;
  const ruleSetColData = ruleSetCol.getRawData();
  const smartsCol: DG.Column<string> = alertsDf.getCol('smarts');
  const smartsColCategories = smartsCol.categories;
  const smartsColData = smartsCol.getRawData();

  // Cache SMARTS
  const smartsList: RDMol[] = new Array(alertsDfLength);
  for (let alertIdx = 0; alertIdx < alertsDfLength; alertIdx++) {
    const currentRuleSet = ruleSetColCategories[ruleSetColData[alertIdx]] as RuleId;
    if (ruleSet[currentRuleSet])
      smartsList[alertIdx] = rdkitModule.get_qmol(smartsColCategories[smartsColData[alertIdx]]);
  }

  // Prepare the result storage
  const skipRuleSetList: string[] = [];
  const moleculeCount = moleculeCol.length;
  const resultValues: {[ruleId in RuleId]?: BitArray} = {};
  for (const [ruleName, isIncluded] of Object.entries(ruleSet)) {
    if (isIncluded)
      resultValues[ruleName as RuleId] = new BitArray(moleculeCount, false);
  }

  // Run the structural alerts detection
  for (let molIdx = 0; molIdx < moleculeCount; molIdx++) {
    const molStr = moleculeColCategories[moleculeColRawData[molIdx]];
    if (molStr === '')
      continue;

    const mol = rdkitModule.get_mol(molStr);

    for (let alertIdx = 0; alertIdx < alertsDfLength; alertIdx++) {
      const currentRuleSet = ruleSetColCategories[ruleSetColData[alertIdx]] as RuleId;
      if (!ruleSet[currentRuleSet] || skipRuleSetList.includes(currentRuleSet))
        continue;

      const matches = mol.get_substruct_match(smartsList[alertIdx]);
      if (matches !== '{}') {
        skipRuleSetList.push(currentRuleSet);
        resultValues[currentRuleSet]!.setTrue(molIdx);
      }
    }
    skipRuleSetList.length = 0;
    mol.delete();
  }

  for (const smarts of smartsList)
    smarts.delete();

  // Build the result dataframe
  const resultDf = DG.DataFrame.create(moleculeCount);
  for (const [ruleName, bitArray] of Object.entries(resultValues))
    resultDf.columns.addNewBool(ruleName).init((i) => bitArray.getBit(i));

  return resultDf;
}
