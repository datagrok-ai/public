import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MMP_NAMES, columnsDescriptions} from './mmp-constants';
import {getRdKitService} from '../../../utils/chem-common-rdkit';
import {MMPA} from '../mmp-analysis/mmpa';

export async function getGenerations(mmpa: MMPA, allPairsGrid: DG.Grid):
  Promise<DG.Grid> {
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


  const initVals = new Float32Array(generation.length);
  const col1 = grid.dataFrame.columns.byName('Existing');
  const col2 = grid.dataFrame.columns.byName('Initial value');
  const col3 = grid.dataFrame.columns.byName('Prediction');
  for (let i = 0; i < generation.length; i++) {
    if (col1.get(i))
      initVals[i] = col3.get(i) - col2.get(i)/4;
    else
      initVals[i] = 0;
  }
  const corrCol = DG.Column.float('~sss', generation.length).init((i) => initVals[i] == 0 ? null : initVals[i]);
  grid.dataFrame.columns.add(corrCol);
  return grid;
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
