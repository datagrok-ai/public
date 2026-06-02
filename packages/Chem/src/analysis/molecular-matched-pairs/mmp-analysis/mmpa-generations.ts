import * as grok from 'datagrok-api/grok';
import {generationsGPU} from '@datagrok-libraries/math/src/webGPU/mmp/webGPU-generations';
import {MMP_CONSTRICTIONS, MMP_ERRORS, MmpFragments} from './mmpa-misc';

const _tsLog = (msg: string): void => console.log(`[${new Date().toISOString()}] ${msg}`);

export async function calculateGenerations(structuresN: number, activityN: number, moleculesArray: string[],
  allStructures: string[], allInitActivities: Float32Array, activityName: string[], activities: Float32Array[],
  activityNames: string[], frags: MmpFragments, meanDiffs: Float32Array[], prediction: Float32Array,
  cores: string[], from: string[], to: string[], rulesFrom: ArrayLike<number>, rulesTo: ArrayLike<number>,
  rulesFromCats: string[], rulesToCats: string[], gpu: boolean, strictCPU: boolean = false) {
  const path = (structuresN < 10 || !gpu || strictCPU) ? 'CPU' : 'GPU';
  _tsLog(`[calculateGenerations] entering, path=${path}, structuresN=${structuresN}, gpu=${gpu}, strictCPU=${strictCPU}`);
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
    _tsLog(`[calculateGenerations] ${path} path completed`);
  } catch (e: any) {
    _tsLog(`[calculateGenerations] ${path} path threw: ${e instanceof Error ? e.message : e}`);
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

async function generationsCPU(structuresN: number, activityN: number, moleculesArray: string[],
  allStructures: string[], allInitActivities: Float32Array, activityName: string[], activities: Float32Array[],
  activityNames: string[], frags: MmpFragments, meanDiffs: Float32Array[], prediction: Float32Array,
  cores: string[], from: string[], to: string[], rulesFrom: ArrayLike<number>, rulesTo: ArrayLike<number>,
  rulesFromCats: string[], rulesToCats: string[]) {
  for (let i = 0; i < structuresN; i ++) {
    for (let j = 0; j < activityN; j++) {
      allStructures[j * structuresN + i] = moleculesArray[i];//mmpInput.molecules.get(i);
      allInitActivities[j * structuresN + i] = activities[j][i];
      activityName[j * structuresN + i] = activityNames[j];
    }

    for (let j = 0; j < frags.fragCodes[i].length; j++) {
      const core = frags.idToName[frags.fragCodes[i][j][0]];
      const subst = frags.idToName[frags.fragCodes[i][j][1]];
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
