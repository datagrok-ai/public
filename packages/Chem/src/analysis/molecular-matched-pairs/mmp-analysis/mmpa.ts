import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {IMmpFragmentsResult} from '../../../rdkit-service/rdkit-service-worker-substructure';
import {getMmpFrags, getMmpRules} from './mmpa-fragments';
import {MMP_CONSTRICTIONS, MMP_ERRORS, MmpRules} from './mmpa-misc';


export class MMPA {
  frags: IMmpFragmentsResult;
  rules: MmpRules;
  allCasesNumber: number;

  constructor(frags: IMmpFragmentsResult, rules: MmpRules, allCasesNumber: number) {
    this.frags = frags;
    this.rules = rules;
    this.allCasesNumber = allCasesNumber;
  }

  static async init(molecules: string[], fragmentCutoff: number): Promise<MMPA> {
    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    if (!gpu && molecules.length > MMP_CONSTRICTIONS.CPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    else if (molecules.length > MMP_CONSTRICTIONS.GPU)
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);

    const frags = await getMmpFrags(molecules);
    const [rules, allCasesNumber] = await getMmpRules(frags, fragmentCutoff, gpu);

    return new MMPA(frags, rules, allCasesNumber);
  }

  static async fromData(data: string): Promise<MMPA> {
    const totalParsed = JSON.parse(data);

    const frags: IMmpFragmentsResult = totalParsed['fragments'];
    const rules: MmpRules = totalParsed['rules'];
    const allCasesNumber: number = totalParsed['cases'];

    return new MMPA(frags, rules, allCasesNumber);
  }

  toJSON() {
    return JSON.stringify({'fragments': this.frags, 'rules': this.rules, 'cases': this.allCasesNumber});
  }
}
