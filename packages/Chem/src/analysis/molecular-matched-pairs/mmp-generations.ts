import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MMP_COLNAME_FROM, MMP_COLNAME_TO, columnsDescriptions} from './mmp-constants';
import {MmpInput} from './mmp-constants';
import {IMmpFragmentsResult} from '../../rdkit-service/rdkit-service-worker-substructure';
import {getRdKitService} from '../../utils/chem-common-rdkit';

export async function getGenerations(mmpInput: MmpInput, moleculesArray: string[],
  fragsOut: IMmpFragmentsResult, meanDiffs: Float32Array[],
  allPairsGrid: DG.Grid, activityMeanNames: Array<string>):
  Promise<DG.Grid> {
  const rulesColumns = allPairsGrid.dataFrame.columns;

  const rulesFrom = rulesColumns.byName(MMP_COLNAME_FROM).getRawData();
  const rulesTo = rulesColumns.byName(MMP_COLNAME_TO).getRawData();
  const rulesFromCats = rulesColumns.byName(MMP_COLNAME_FROM).categories;
  const rulesToCats = rulesColumns.byName(MMP_COLNAME_TO).categories;
  const activityN = activityMeanNames.length;

  //const structures = mmpInput.molecules.toList();
  //const structuresN = structures.length;
  // const structures = mmpInput.molecules.getRawData();
  // const structuresCats = mmpInput.molecules.categories;
  const structuresN = mmpInput.molecules.length;
  const cores = new Array<string>(activityN *structuresN);
  const from = new Array<string>(activityN * structuresN);
  const to = new Array<string>(activityN * structuresN);
  const prediction = new Float32Array(activityN * structuresN).fill(0);
  const allStructures = Array(structuresN * activityN);
  const allInitActivities = new Float32Array(structuresN * activityN);
  const activityName: Array<string> = Array(structuresN * activityN);
  const activities = new Array<Float32Array>(activityN);
  const activityNames = Array<string>(activityN);
  for (let i = 0; i < activityN; i++) {
    const name = activityMeanNames[i].replace('Mean Difference ', '');
    activities[i] = mmpInput.activities.byName(name).getRawData() as Float32Array;
    activityNames[i] = name;
  }

  for (let i = 0; i < structuresN; i ++) {
    for (let j = 0; j < activityN; j++) {
      allStructures[j * structuresN + i] = moleculesArray[i];//mmpInput.molecules.get(i);
      allInitActivities[j * structuresN + i] = activities[j][i];
      activityName[j * structuresN + i] = activityNames[j];
    }

    for (let j = 0; j < fragsOut.frags[i].length; j++) {
      const core = fragsOut.frags[i][j][0];
      const subst = fragsOut.frags[i][j][1];
      if (core != '') {
        for (let k = 0; k < rulesFrom.length; k++) {
          if (subst === rulesFromCats[rulesFrom[k]] ) {
            for (let kk = 0; kk < activityN; kk++) {
              const activity = activities[kk][i] + meanDiffs[kk][k];
              const add = kk * structuresN;
              if (activity > prediction[add + i]) {
                prediction[add + i] = activity;
                cores[add + i] = core;
                from[add + i] = subst;
                to[add + i] = rulesToCats[rulesTo[k]];
              }
            }
          }
        }
      }
    }
  }

  const generation = await (await getRdKitService()).mmpLinkFragments(cores, to);
  const cols = [];
  cols.push(createColWithDescription('string', 'Structure', allStructures, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('double', `Initial value`, Array.from(allInitActivities)));
  cols.push(createColWithDescription('string', `Activity`, activityName));
  cols.push(createColWithDescription('string', `Core`, cores, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', `From`, from, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('string', `To`, to, DG.SEMTYPE.MOLECULE));
  cols.push(createColWithDescription('double', `Prediction`, Array.from(prediction)));
  cols.push(createColWithDescription('string', `Generation`, generation, DG.SEMTYPE.MOLECULE));
  const grid = DG.DataFrame.fromColumns(cols).plot.grid();
  createMolExistsCol(fragsOut.smiles, generation, grid);
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
