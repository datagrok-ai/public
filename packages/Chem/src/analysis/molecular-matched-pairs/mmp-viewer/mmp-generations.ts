import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MMP_NAMES, columnsDescriptions} from './mmp-constants';
import {getRdKitService} from '../../../utils/chem-common-rdkit';
import {MMPA} from '../mmp-analysis/mmpa';

export async function getGenerations(mmpa: MMPA, allPairsGrid: DG.Grid):
  Promise<[DG.Grid, DG.Grid]> {
  const rulesColumns = allPairsGrid.dataFrame.columns;

  const rulesFrom = rulesColumns.byName(MMP_NAMES.FROM).getRawData();
  const rulesTo = rulesColumns.byName(MMP_NAMES.TO).getRawData();
  const rulesFromCats = rulesColumns.byName(MMP_NAMES.FROM).categories;
  const rulesToCats = rulesColumns.byName(MMP_NAMES.TO).categories;

  const genRes = await mmpa.calculateGenerations(rulesFrom, rulesTo, rulesFromCats, rulesToCats);

  const generation = await (await getRdKitService()).mmpLinkFragments(genRes.cores, genRes.to);
  const cols = [];
  cols.push(createColWithDescription('string', 'Structure', genRes.allStructures, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('double', `Initial value`, Array.from(genRes.allInitActivities)));
  cols.push(createColWithDescription('string', `Activity`, genRes.activityName));
  cols.push(createColWithDescription('string', `Core`, genRes.cores, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', `From`, genRes.from, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', `To`, genRes.to, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('double', `Prediction`, Array.from(genRes.prediction)));
  cols.push(createColWithDescription('string', `Generation`, generation, DG.SEMTYPE.MOLECULE));
  const grid = DG.DataFrame.fromColumns(cols).plot.grid();
  createMolExistsCol(mmpa.frags.smiles, generation, grid);

  const gridCorr = getCorGrid(mmpa);

  return [grid, gridCorr];
}

export function createColWithDescription(colType: any, colName: string, list: any[], semType?: DG.SemType): DG.Column {
  const col = DG.Column.fromList(colType, colName, list);
  if (columnsDescriptions[colName])
    col.setTag('description', columnsDescriptions[colName]);
  if (semType)
    col.semType = semType;
  return col;
}

export async function createMolExistsCol(molecules: string[], generation: string[], grid: DG.Grid): Promise<void> {
  const moleculesSet = new Set(molecules);
  const boolCol = DG.Column.bool('Existing', generation.length).init((i) => moleculesSet.has(generation[i]));
  grid.dataFrame.columns.add(boolCol);
  grid.col(boolCol.name)!.editable = false;
  grid.invalidate();
}

function getCorGrid(mmpa: MMPA): DG.Grid {
  const activitiesCount = mmpa.initData.activitiesCount;
  const length = mmpa.allCasesNumber*activitiesCount;
  const initialMol = new Array<string>(length);
  const initialFragment = new Array<string>(length);
  const endFragment = new Array<string>(length);
  const activity = new Array<string>(length);
  const observed = new Array<Number>(length);
  const predicted = new Array<Number>(length);

  let counter = 0;
  const rules = mmpa.rules.rules;
  for (let i = 0; i < rules.length; i ++) {
    for (let j = 0; j <rules[i].pairs.length; j++) {
      for (let k = 0; k < activitiesCount; k++) {
        initialMol[counter] = rules[i].pairs[j].core;
        initialFragment[counter] = mmpa.rules.smilesFrags[rules[i].smilesRule1];
        endFragment[counter] = mmpa.rules.smilesFrags[rules[i].smilesRule2];
        activity[counter] = mmpa.initData.activitiesNames[k];
        observed[counter] = mmpa.initData.activities[k][rules[i].pairs[j].secondStructure];
        predicted[counter] = mmpa.initData.activities[k][rules[i].pairs[j].firstStructure] +
          mmpa.rulesBased.meanDiffs[k][i];
        counter++;
      }
    }
  }

  const colsCorr = [];
  colsCorr.push(createColWithDescription('string', 'Core', initialMol, DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', 'From', initialFragment, DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', 'To', endFragment, DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', 'Activity', activity));
  colsCorr.push(createColWithDescription('double', 'Observed', observed));
  colsCorr.push(createColWithDescription('double', 'Predicted', predicted));
  return DG.DataFrame.fromColumns(colsCorr).plot.grid();
}
