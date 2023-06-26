import {OCLServiceCall} from './consts';
import {CHEM_PROP_MAP} from './calculations';
import * as OCL from 'openchemlib/full';
onmessage = ({data: {op, data, argList}}) => {
  switch (op) {
  case OCLServiceCall.CHEM_PROPERTIES:
    postMessage(getChemProperties(data, argList));
    break;
  default:
    postMessage(getToxRisks(data, argList));
    break;
  }
};


function getChemProperties(molList: Array<string>, propList: string[]) {
  const res: {[key: string]: Array<number>} = {};
  const errors: string[] = [];
  propList.forEach((propName) => {
    res[propName] = new Array(molList.length).fill(NaN);
  });
  molList.forEach((smiles, i) => {
    try {
      const mol = OCL.Molecule.fromSmiles(smiles);
      propList.forEach((p) => {
        res[p][i] = CHEM_PROP_MAP[p].valueFunc(mol);
      });
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

function getToxRisks(molList: Array<string>, riskTypes: number[]) {
  const res: {[key: string]: Array<number>} = {};
  const errors: string[] = [];
  riskTypes.forEach((propName) => {
    res[propName] = new Array(molList.length).fill(0);
  });
  const toxicityPredictor = new OCL.ToxicityPredictor();
  molList.forEach((smiles, i) => {
    try {
      const mol = OCL.Molecule.fromSmiles(smiles);
      riskTypes.forEach((p) => {
        res[p][i] = toxicityPredictor.assessRisk(mol, p);
      });
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

