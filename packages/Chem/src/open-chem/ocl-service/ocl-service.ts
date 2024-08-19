import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {MolNotationType, OCLServiceCall} from './consts';

export class OCLService {
  _threadCount: number = Math.max((navigator.hardwareConcurrency || 6) - 2, 1);
  _workers: Worker[] = new Array(this._threadCount).fill(null)
    .map(() => new Worker(new URL('workers/ocl-worker', import.meta.url)));
  constructor() {
  }

  private async _doParallel(
    molCol:DG.Column, argList: any[], op: OCLServiceCall,
  ): Promise<{[key: string]: Array<number>}> {
    const molList = molCol.toList();
    // detect if the column is smiles or molblock notation
    const notationType = DG.chem.isMolBlock(molList[0]) ? MolNotationType.MOLBLOCK : MolNotationType.SMILES;

    const colLenght = molList.length;
    const chunkSize = Math.ceil(colLenght / this._threadCount);
    const result: {[key: string]: Array<number>} = {};
    const promises: Promise<{res:{[key: string]: Array<number>}, errors: string[]}>[] = new Array(this._threadCount);
    for (let i = 0; i < this._threadCount; i++) {
      const chunk = molList.slice(i * chunkSize, (i === (this._threadCount - 1)) ? colLenght : (i + 1) * chunkSize);
      this._workers[i].postMessage({op, data: chunk, argList, notationType});
      promises[i] = new Promise((resolve) => {
        this._workers[i].onmessage = ({data: {res, errors}}) => {
          resolve({res, errors});
        };
      });
    }
    const scatterRes = await Promise.all(promises);
    const totalErrors = scatterRes.map(({errors}) => errors.length)
      .reduce((prev, cur) => prev + cur, 0);
    const errorPercentages = totalErrors / colLenght;
    if (errorPercentages > 0.1)
      grok.shell.warning(`Operation resulted in error for ${totalErrors} molecules`);

    scatterRes.forEach((res) => {
      Object.keys(res.res).forEach((propName) => {
        if (!result[propName])
          result[propName] = [];

        result[propName].push(...res.res[propName]);
      });
    });
    return result;
  }


  async getChemProperties(molCol: DG.Column, propList: string[]): Promise<{[key: string]: Array<number>}> {
    return this._doParallel(molCol, propList, OCLServiceCall.CHEM_PROPERTIES);
  }

  async getChemToxicity(molCol: DG.Column, riskTypes: number[]): Promise<{[key: string]: Array<number>}> {
    return this._doParallel(molCol, riskTypes, OCLServiceCall.TOXICITY);
  }

  async getDrugLikeness(molCol: DG.Column): Promise<{[key: string]: Array<number>}> {
    return this._doParallel(molCol, [], OCLServiceCall.DRUG_LIKENESS);
  }

  public async terminate() {
    this._workers.forEach((worker) => worker.terminate());
  }
}
