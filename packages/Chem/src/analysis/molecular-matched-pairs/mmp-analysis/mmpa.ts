import * as DG from 'datagrok-api/dg';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {IMmpFragmentsResult} from '../../../rdkit-service/rdkit-service-worker-substructure';
import {MMP_CONSTRICTIONS, MMP_ERRORS, MmpRules} from './mmpa-misc';
import {chemSpace} from '../../chem-space';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

import {getMmpFrags, getMmpRules} from './mmpa-fragments';
import {getPlainData, MmpAllCasesBasedData, MmpRulesBasedData} from './mmpa-differences';


export class MMPA {
  frags: IMmpFragmentsResult;
  rules: MmpRules;
  allCasesNumber: number;

  rulesBased: MmpRulesBasedData;
  allCasesBased: MmpAllCasesBasedData;

  chemSpaceResult: Float32Array [] | null = null;

  constructor(frags: IMmpFragmentsResult, rules: MmpRules, allCasesNumber: number,
    rulesBased: MmpRulesBasedData, allCasesBased: MmpAllCasesBasedData,
    chemSpaceResult: Float32Array [] | null) {
    this.frags = frags;
    this.rules = rules;
    this.allCasesNumber = allCasesNumber;

    this.rulesBased = rulesBased;
    this.allCasesBased = allCasesBased;

    this.chemSpaceResult = chemSpaceResult;
  }

  static async init(molecules: string[], fragmentCutoff: number, activities: Float32Array[]): Promise<MMPA> {
    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    if (!gpu && molecules.length > MMP_CONSTRICTIONS.CPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    else if (molecules.length > MMP_CONSTRICTIONS.GPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);

    const frags = await getMmpFrags(molecules);
    const [rules, allCasesNumber] = await getMmpRules(frags, fragmentCutoff, gpu);

    const [rulesBased, allCasesBased] = getPlainData(rules, molecules, activities, allCasesNumber);

    return new MMPA(frags, rules, allCasesNumber, rulesBased, allCasesBased, null);
  }

  static async fromData(data: string, molecules: string[], activities: Float32Array[]): Promise<MMPA> {
    const totalParsed = JSON.parse(data);

    const frags: IMmpFragmentsResult = totalParsed['fragments'];
    const rules: MmpRules = totalParsed['rules'];
    const allCasesNumber: number = totalParsed['cases'];

    const [rulesBased, allCasesBased] = getPlainData(rules, molecules, activities, allCasesNumber);

    const chemSpaceResult: Float32Array [] | null = totalParsed['chemSpaceResult'];

    return new MMPA(frags, rules, allCasesNumber, rulesBased, allCasesBased, chemSpaceResult);
  }

  toJSON() {
    return JSON.stringify({
      'fragments': this.frags,
      'rules': this.rules,
      'cases': this.allCasesNumber,
      'chemSpaceResult': this.chemSpaceResult,
    });
  }

  async chemSpace(chemSpaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
    if (this.chemSpaceResult) {
      const cols: DG.Column[] = chemSpaceParams.embedAxesNames.map((name: string, index: number) =>
        DG.Column.fromFloat32Array(name, this.chemSpaceResult![index]));
      return {coordinates: new DG.ColumnList(cols)};
    } else {
      const res = await chemSpace(chemSpaceParams);

      //const embeddings = res.coordinates;

      // this.chemSpaceResult = new Array<Float32Array>(2);
      // this.chemSpaceResult[0] = embeddings.byIndex(0).asDoubleList();
      // this.chemSpaceResult[1] = embeddings.byIndex(1).asDoubleList();

      return res;
    }
  }

  async calculateGenerations() {
    
  }
}
