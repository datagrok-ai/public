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

const _tsLog = (msg: string): void => console.log(`[${new Date().toISOString()}] ${msg}`);

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
    activities: Float32Array[], activitiesNames: string[], diffTypes: string[],
    fragSortingInfo: SortData): Promise<MMPA> {
    _tsLog(`[MMPA.init] entering, molecules=${molecules.length}, fragmentCutoff=${fragmentCutoff}, ` +
      `activities=${activitiesNames.length}`);
    const initData: MmpInitData =
    {molName, molecules, activities, activitiesNames, activitiesCount: activitiesNames.length, diffTypes};

    _tsLog('[MMPA.init] calling getGPUDevice');
    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;
    _tsLog(`[MMPA.init] getGPUDevice returned, gpu=${gpu}`);

    if (!gpu && molecules.length > MMP_CONSTRICTIONS.CPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    else if (molecules.length > MMP_CONSTRICTIONS.GPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);

    _tsLog('[MMPA.init] calling getMmpFrags');
    const [frags, canonical] = await getMmpFrags(molecules);
    _tsLog(`[MMPA.init] getMmpFrags returned, frags.idToName.length=${frags?.idToName?.length}, ` +
      `canonical.length=${canonical?.length}`);
    initData.molecules = canonical;
    _tsLog('[MMPA.init] calling getMmpRules');
    const [rules, allCasesNumber] = await getMmpRules(frags, fragmentCutoff, gpu);
    _tsLog(`[MMPA.init] getMmpRules returned, rules.rules.length=${rules?.rules?.length}, ` +
      `smilesFrags=${rules?.smilesFrags?.length}, allCasesNumber=${allCasesNumber}`);

    _tsLog('[MMPA.init] calling getPlainData');
    const [rulesBased, allCasesBased] = getPlainData(rules, frags, initData, allCasesNumber, fragSortingInfo);
    _tsLog('[MMPA.init] getPlainData returned, constructing MMPA');

    return new MMPA(initData, frags, rules, allCasesNumber, rulesBased, allCasesBased, null, null, gpu);
  }

  static async fromData(molName: string, data: string, molecules: string[],
    activities: Float32Array[], activitiesNames: string[], diffTypes: string[],
    fragSortingInfo: SortData): Promise<MMPA> {
    _tsLog(`[MMPA.fromData] entering, data.length=${data?.length}, molecules=${molecules.length}`);
    const initData: MmpInitData =
    {molName, molecules, activities, activitiesNames, activitiesCount: activitiesNames.length, diffTypes};

    _tsLog('[MMPA.fromData] calling getGPUDevice');
    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;
    _tsLog(`[MMPA.fromData] getGPUDevice returned, gpu=${gpu}, parsing JSON`);

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

    _tsLog(`[MMPA.fromData] done, rules.rules.length=${rules?.rules?.length}, ` +
      `chemSpaceResult=${!!chemSpaceResult}, generationResult=${!!generationResult}`);
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
    _tsLog(`[MMPA.chemSpace] entering, cached=${!!this.chemSpaceResult}`);
    if (this.chemSpaceResult) {
      const cols: DG.Column[] = chemSpaceParams.embedAxesNames.map((name: string, index: number) =>
        DG.Column.fromFloat32Array(name, this.chemSpaceResult![index]));
      _tsLog('[MMPA.chemSpace] returning cached embeddings');
      return {coordinates: new DG.ColumnList(cols)};
    } else {
      _tsLog('[MMPA.chemSpace] calling chemSpace() (dim. reduction)');
      const res = await chemSpace(chemSpaceParams);
      _tsLog('[MMPA.chemSpace] chemSpace() returned');

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
    _tsLog(`[MMPA.calculateGenerations] entering, cached=${!!this.generationResult}, gpu=${this.gpu}`);
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

      _tsLog(`[MMPA.calculateGenerations] calling calculateGenerations(), structuresN=${structuresN}, ` +
        `activityN=${activityN}`);
      await calculateGenerations(structuresN, activityN, this.initData.molecules, allStructures, allInitActivities,
        activityName, this.initData.activities, this.initData.activitiesNames,
        this.frags, this.rulesBased.meanDiffs,
        prediction, cores, from, to,
        rulesFrom, rulesTo, rulesFromCats, rulesToCats, this.gpu);
      _tsLog('[MMPA.calculateGenerations] calculateGenerations() returned');

      this.generationResult = {allStructures, allInitActivities, activityName, cores, from, to, prediction};
      return this.generationResult;
    }
  }
}
