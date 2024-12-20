import * as DG from 'datagrok-api/dg';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {chemSpace} from '../../chem-space';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

import {MMP_CONSTRICTIONS, MMP_ERRORS,
  MmpRules, MmpInitData, MmpAllCasesBasedData, MmpRulesBasedData, MmpGeneration,
  MmpFragments} from './mmpa-misc';
import {getMmpFrags, getMmpRules} from './mmpa-fragments';
import {getPlainData} from './mmpa-differences';
import {calculateGenerations} from './mmpa-generations';
import {SortData} from '../mmp-viewer/mmp-viewer';

export class MMPA {
  initData: MmpInitData;
  gpu: boolean;

  frags: MmpFragments;
  rules: MmpRules;
  allCasesNumber: number;

  rulesBased: MmpRulesBasedData;
  allCasesBased: MmpAllCasesBasedData;

  chemSpaceResult: Float32Array [] | null = null;
  generationResult: MmpGeneration | null = null;

  constructor(initData: MmpInitData, frags: MmpFragments, rules: MmpRules, allCasesNumber: number,
    rulesBased: MmpRulesBasedData, allCasesBased: MmpAllCasesBasedData,
    chemSpaceResult: Float32Array [] | null, generationResult: MmpGeneration | null, gpu: boolean) {
    this.initData = initData;
    this.gpu = gpu;

    this.frags = frags;
    this.rules = rules;
    this.allCasesNumber = allCasesNumber;

    this.rulesBased = rulesBased;
    this.allCasesBased = allCasesBased;

    this.chemSpaceResult = chemSpaceResult;
    this.generationResult = generationResult;
  }

  static async init(molName: string, molecules: string[], fragmentCutoff: number,
    activities: Float32Array[], activitiesNames: string[], fragSortingInfo: SortData): Promise<MMPA> {
    const initData: MmpInitData =
    {molName, molecules, activities, activitiesNames, activitiesCount: activitiesNames.length};

    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    if (!gpu && molecules.length > MMP_CONSTRICTIONS.CPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    else if (molecules.length > MMP_CONSTRICTIONS.GPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);

    const [frags, canonical] = await getMmpFrags(molecules);
    initData.molecules = canonical;
    const [rules, allCasesNumber] = await getMmpRules(frags, fragmentCutoff, gpu);

    const [rulesBased, allCasesBased] = getPlainData(rules, frags, initData, allCasesNumber, fragSortingInfo);

    return new MMPA(initData, frags, rules, allCasesNumber, rulesBased, allCasesBased, null, null, gpu);
  }

  static async fromData(molName: string, data: string, molecules: string[],
    activities: Float32Array[], activitiesNames: string[], fragSortingInfo: SortData): Promise<MMPA> {
    const initData: MmpInitData =
    {molName, molecules, activities, activitiesNames, activitiesCount: activitiesNames.length};

    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    const totalParsed = JSON.parse(data);

    const frags: MmpFragments = totalParsed['fragments'];
    const rules: MmpRules = totalParsed['rules'];
    const allCasesNumber: number = totalParsed['cases'];

    const [rulesBased, allCasesBased] = getPlainData(rules, frags, initData, allCasesNumber, fragSortingInfo);

    const chemSpaceResult: Float32Array [] | null = totalParsed['chemSpaceResult'];
    const generationResult: MmpGeneration | null = totalParsed['generationResult'];

    if (chemSpaceResult) {
      chemSpaceResult[0] = new Float32Array(Object.values(chemSpaceResult[0]));
      chemSpaceResult[1] = new Float32Array(Object.values(chemSpaceResult[1]));
    }

    if (generationResult) {
      generationResult.allInitActivities = new Float32Array(Object.values(generationResult.allInitActivities));
      generationResult.prediction = new Float32Array(Object.values(generationResult.prediction));
    }

    return new MMPA(initData, frags, rules, allCasesNumber, rulesBased, allCasesBased,
      chemSpaceResult, generationResult, gpu);
  }

  toJSON() {
    return JSON.stringify({
      'fragments': this.frags,
      'rules': this.rules,
      'cases': this.allCasesNumber,
      'chemSpaceResult': this.chemSpaceResult,
      'generationResult': this.generationResult,
    }, function(key, val) {
      return val && !isNaN(Number(val)) && val.toFixed ? Number(val.toFixed(5)) : val;
    });
  }

  async chemSpace(chemSpaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
    if (this.chemSpaceResult) {
      const cols: DG.Column[] = chemSpaceParams.embedAxesNames.map((name: string, index: number) =>
        DG.Column.fromFloat32Array(name, this.chemSpaceResult![index]));
      return {coordinates: new DG.ColumnList(cols)};
    } else {
      const res = await chemSpace(chemSpaceParams);

      const embeddings = res.coordinates;
      const cols: DG.Column [] = [];
      for (const col of embeddings)
        cols.push(col);

      this.chemSpaceResult = new Array<Float32Array>(2);
      this.chemSpaceResult[0] = cols[0].asDoubleList();
      this.chemSpaceResult[1] = cols[1].asDoubleList();

      return res;
    }
  }

  async calculateGenerations(rulesFrom: ArrayLike<number>, rulesTo: ArrayLike<number>,
    rulesFromCats: string[], rulesToCats: string[]): Promise<MmpGeneration> {
    if (this.generationResult)
      return this.generationResult;
    else {
      const activityN = this.initData.activitiesCount;
      const structuresN = this.initData.molecules.length;

      const allStructures = Array(structuresN * activityN);
      const allInitActivities = new Float32Array(structuresN * activityN);
      const cores = new Array<string>(activityN *structuresN);
      const from = new Array<string>(activityN * structuresN);
      const to = new Array<string>(activityN * structuresN);
      const prediction = new Float32Array(activityN * structuresN).fill(0);
      const activityName: Array<string> = Array(structuresN * activityN);

      await calculateGenerations(structuresN, activityN, this.initData.molecules, allStructures, allInitActivities,
        activityName, this.initData.activities, this.initData.activitiesNames,
        this.frags, this.rulesBased.meanDiffs,
        prediction, cores, from, to,
        rulesFrom, rulesTo, rulesFromCats, rulesToCats, this.gpu);

      this.generationResult = {allStructures, allInitActivities, activityName, cores, from, to, prediction};
      return this.generationResult;
    }
  }
}
