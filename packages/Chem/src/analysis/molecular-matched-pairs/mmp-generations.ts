import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MMP_COLNAME_FROM, MMP_COLNAME_TO, columnsDescriptions} from './mmp-constants';
import {MmpInput} from './mmp-constants';
import {IMmpFragmentsResult} from '../../rdkit-service/rdkit-service-worker-substructure';
import {getRdKitService} from '../../utils/chem-common-rdkit';

function encodeFragments2(fragsOut: IMmpFragmentsResult):
  [[number, number][][], Record<string, number>, Array<string>] {
  const res: [number, number][][] = new Array(fragsOut.frags.length)
    .fill(null).map((_, i) => new Array(fragsOut.frags[i].length).fill(null).map((_) => [0, 0]));

  let fragmentCounter = 0;
  const fragmentMap: Record<string, number> = {};
  for (let i = 0; i < fragsOut.frags.length; i++) {
    for (let j = 0; j < fragsOut.frags[i].length; j++) {
      if (!fragmentMap[fragsOut.frags[i][j][0]]) {
        fragmentMap[fragsOut.frags[i][j][0]] = fragmentCounter;
        fragmentCounter++;
      }
      if (!fragmentMap[fragsOut.frags[i][j][1]]) {
        fragmentMap[fragsOut.frags[i][j][1]] = fragmentCounter;
        fragmentCounter++;
      }
      res[i][j][0] = fragmentMap[fragsOut.frags[i][j][0]];
      res[i][j][1] = fragmentMap[fragsOut.frags[i][j][1]];
    }
  }

  const maxFragmentIndex = fragmentCounter;
  const fragIdToFragName = new Array<string>(maxFragmentIndex);
  Object.entries(fragmentMap).forEach(([key, val]) => {
    fragIdToFragName[val] = key;
  });

  // const fragSizes = new Uint32Array(maxFragmentIndex).map((_, i) => fragIdToFragName[i].length);


  return [res, fragmentMap, fragIdToFragName];
}

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

  const [encodedFrags, fragmentMap, fragIdToFragName] = encodeFragments2(fragsOut);
  const emptyCore = fragmentMap[''];
  const encodedFrom = new Int32Array(rulesFrom.length);
  const encodedTo = new Int32Array(rulesFrom.length);
  for (let i = 0; i < rulesFrom.length; i++) {
    encodedFrom[i] = fragmentMap[rulesFromCats[rulesFrom[i]]];
    encodedTo[i] = fragmentMap[rulesToCats[rulesTo[i]]];
  }

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

  console.profile('MMP');
  for (let i = 0; i < structuresN; i ++) {
    for (let j = 0; j < activityN; j++) {
      const idx = j * structuresN + i;
      allStructures[idx] = moleculesArray[i];
      allInitActivities[idx] = activities[j][i];
      activityName[idx] = activityNames[j];
    }

    for (let j = 0; j < fragsOut.frags[i].length; j++) {
      const core = encodedFrags[i][j][0];
      const subst = encodedFrags[i][j][1];
      if (core != emptyCore) {
        for (let k = 0; k < rulesFrom.length; k++) {
          if (subst === encodedFrom[k]) {
            for (let kk = 0; kk < activityN; kk++) {
              const activity = activities[kk][i] + meanDiffs[kk][k];
              const idx = kk * structuresN + i;
              if (activity > prediction[idx + i]) {
                prediction[idx] = activity;
                cores[idx] = fragIdToFragName[core];
                from[idx] = fragIdToFragName[subst];
                to[idx] = fragIdToFragName[k];
              }
            }
          }
        }
      }
    }
  }

  console.profileEnd('MMP');

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
