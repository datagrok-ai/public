import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GENERATIONS_GRID_HEADER_DESCRIPTIONS, MMP_NAMES} from './mmp-constants';
import {getRdKitService} from '../../../utils/chem-common-rdkit';
import {MMPA} from '../mmp-analysis/mmpa';
import {resizeGridColsSize} from '../../../utils/ui-utils';

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
  cols.push(createColWithDescription('string', MMP_NAMES.STRUCTURE, genRes.allStructures, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', MMP_NAMES.CORE, genRes.cores, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', MMP_NAMES.FROM, genRes.from, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', MMP_NAMES.TO, genRes.to, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', MMP_NAMES.NEW_MOLECULE, generation, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', MMP_NAMES.PROPERTY_TYPE, genRes.activityName, GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  cols.push(createColWithDescription('double', MMP_NAMES.OBSERVED, Array.from(genRes.allInitActivities), GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  cols.push(createColWithDescription('double', MMP_NAMES.PREDICTED, Array.from(genRes.prediction), GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  cols.push(createColWithDescription('double', MMP_NAMES.DELTA_ACTIVITY,
    Array.from(genRes.prediction).map((val: number, idx: number) => val - genRes.allInitActivities[idx]), GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  const grid = DG.DataFrame.fromColumns(cols).plot.grid();

  grid.onCellTooltip(function(cell: DG.GridCell, x: number, y: number) {
    if (cell.isColHeader && cell.tableColumn && cell.tableColumn.semType === DG.SEMTYPE.MOLECULE) {
      let tooltip = '';
      if (GENERATIONS_GRID_HEADER_DESCRIPTIONS[cell.tableColumn.name])
        tooltip = GENERATIONS_GRID_HEADER_DESCRIPTIONS[cell.tableColumn.name];
      ui.tooltip.show(ui.divText(tooltip), x, y);
      return true;
    } else
      return false;
  });

  const gridSub = grid.onAfterDrawContent.subscribe(() => {
    resizeGridColsSize(grid, [MMP_NAMES.STRUCTURE, MMP_NAMES.CORE, MMP_NAMES.FROM, MMP_NAMES.TO], 150, 70);
    gridSub.unsubscribe();
  });
  createMolExistsCol(mmpa.initData.molecules, generation, grid);

  const gridCorr = getCorGrid(mmpa);

  return [grid, gridCorr];
}

export function createColWithDescription(colType: any, colName: string, list: any, descriptions: {[key: string]: string}, type: string,
  semType?: DG.SemType, mean?: boolean): DG.Column {
  const col = type === 'int32' ? DG.Column.fromInt32Array(colName, list) : type === 'float32' ?
    DG.Column.fromFloat32Array(colName, list) : DG.Column.fromList(colType, colName, list);
  if (!semType) {
    if (descriptions[colName])
      col.setTag('description', descriptions[colName]);
    else if (colName.includes('\u0394'))
      col.setTag('description', colName.replace('\u0394', mean ? 'Mean difference in' : 'Difference in'));
  } else
    col.semType = semType;
  return col;
}

export async function createMolExistsCol(molecules: string[], generation: string[], grid: DG.Grid): Promise<void> {
  const moleculesSet = new Set(molecules);
  const boolCol = DG.Column.bool(MMP_NAMES.PREDICTION, generation.length).init((i) => !moleculesSet.has(generation[i]));
  boolCol.setTag('description', GENERATIONS_GRID_HEADER_DESCRIPTIONS[MMP_NAMES.PREDICTION]);
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
        initialMol[counter] = mmpa.frags.idToName[rules[i].pairs[j].core];
        initialFragment[counter] = mmpa.frags.idToName[mmpa.rules.smilesFrags[rules[i].sr1]];
        endFragment[counter] = mmpa.frags.idToName[mmpa.rules.smilesFrags[rules[i].sr2]];
        activity[counter] = mmpa.initData.activitiesNames[k];
        observed[counter] = mmpa.initData.activities[k][rules[i].pairs[j].ss];
        predicted[counter] = mmpa.initData.activities[k][rules[i].pairs[j].fs] +
          mmpa.rulesBased.meanDiffs[k][i];
        counter++;
      }
    }
  }

  const colsCorr = [];
  colsCorr.push(createColWithDescription('string', MMP_NAMES.CORE, initialMol, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', MMP_NAMES.FROM, initialFragment, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', MMP_NAMES.TO, endFragment, GENERATIONS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE));
  colsCorr.push(createColWithDescription('string', MMP_NAMES.ACTIVITY, activity, GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  colsCorr.push(createColWithDescription('double', 'Observed', observed, GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  colsCorr.push(createColWithDescription('double', 'Predicted', predicted, GENERATIONS_GRID_HEADER_DESCRIPTIONS, ''));
  return DG.DataFrame.fromColumns(colsCorr).plot.grid();
}
