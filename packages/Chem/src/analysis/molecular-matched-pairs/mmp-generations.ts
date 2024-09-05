import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MMP_NAMES, columnsDescriptions} from './mmp-viewer/mmp-constants';
import {MmpInput} from './mmp-viewer/mmp-viewer';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {generationsGPU} from '@datagrok-libraries/math/src/webGPU/mmp/webGPU-generations';
import {MMP_CONSTRICTIONS, MMP_ERRORS} from './mmp-analysis/mmpa-misc';
import {MMPA} from './mmp-analysis/mmpa';

export async function getGenerations(mmpInput: MmpInput, moleculesArray: string[],
  mmpa: MMPA,
  allPairsGrid: DG.Grid, activityMeanNames: Array<string>, gpu: boolean):
  Promise<DG.Grid> {
  const rulesColumns = allPairsGrid.dataFrame.columns;

  const rulesFrom = rulesColumns.byName(MMP_NAMES.FROM).getRawData();
  const rulesTo = rulesColumns.byName(MMP_NAMES.TO).getRawData();
  const rulesFromCats = rulesColumns.byName(MMP_NAMES.FROM).categories;
  const rulesToCats = rulesColumns.byName(MMP_NAMES.TO).categories;


  const activityN = activityMeanNames.length;

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


  await calculateGenerations(structuresN, activityN, moleculesArray, allStructures, allInitActivities,
    activityName, activities, activityNames, mmpa.frags.frags, mmpa.rulesBased.meanDiffs,
    prediction, cores, from, to,
    rulesFrom, rulesTo, rulesFromCats, rulesToCats, gpu);

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
  createMolExistsCol(mmpa.frags.smiles, generation, grid);
  return grid;
}

async function calculateGenerations(structuresN: number, activityN: number, moleculesArray: string[],
  allStructures: string[], allInitActivities: Float32Array, activityName: string[], activities: Float32Array[],
  activityNames: string[], frags: [string, string][][], meanDiffs: Float32Array[], prediction: Float32Array,
  cores: string[], from: string[], to: string[], rulesFrom: ArrayLike<number>, rulesTo: ArrayLike<number>,
  rulesFromCats: string[], rulesToCats: string[], gpu: boolean, strictCPU: boolean = false) {
  try {
    if (structuresN < 10 || !gpu || strictCPU) {
      if (structuresN > MMP_CONSTRICTIONS.CPU)
        throw new Error(MMP_ERRORS.FRAGMENTS_CPU);

      await generationsCPU(structuresN, activityN, moleculesArray, allStructures, allInitActivities,
        activityName, activities, activityNames, frags, meanDiffs, prediction, cores, from, to,
        rulesFrom, rulesTo, rulesFromCats, rulesToCats);
    } else {
      await generationsGPU(structuresN, activityN, moleculesArray, allStructures, allInitActivities,
        activityName, activities, activityNames, frags, meanDiffs, prediction, cores, from, to,
        rulesFrom, rulesTo, rulesFromCats, rulesToCats);
    }
  } catch (e: any) {
    const eMsg: string = e instanceof Error ? e.message : e.toString();
    if (eMsg === MMP_ERRORS.FRAGMENTS_CPU) {
      grok.shell.warning(MMP_ERRORS.GPU_ABORTED);
      grok.shell.error(MMP_ERRORS.FRAGMENTS_CPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    }
    if (gpu) {
      await calculateGenerations(structuresN, activityN, moleculesArray, allStructures, allInitActivities,
        activityName, activities, activityNames, frags, meanDiffs, prediction, cores, from, to,
        rulesFrom, rulesTo, rulesFromCats, rulesToCats, gpu, true);
    } else {
      grok.shell.error(MMP_ERRORS.GENERATIONS);
      throw new Error(MMP_ERRORS.GENERATIONS);
    }
  }
}

// eslint-disable-next-line max-params
async function generationsCPU(structuresN: number, activityN: number, moleculesArray: string[],
  allStructures: string[], allInitActivities: Float32Array, activityName: string[], activities: Float32Array[],
  activityNames: string[], frags: [string, string][][], meanDiffs: Float32Array[], prediction: Float32Array,
  cores: string[], from: string[], to: string[], rulesFrom: ArrayLike<number>, rulesTo: ArrayLike<number>,
  rulesFromCats: string[], rulesToCats: string[]) {
  for (let i = 0; i < structuresN; i ++) {
    for (let j = 0; j < activityN; j++) {
      allStructures[j * structuresN + i] = moleculesArray[i];//mmpInput.molecules.get(i);
      allInitActivities[j * structuresN + i] = activities[j][i];
      activityName[j * structuresN + i] = activityNames[j];
    }

    for (let j = 0; j < frags[i].length; j++) {
      const core = frags[i][j][0];
      const subst = frags[i][j][1];
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
