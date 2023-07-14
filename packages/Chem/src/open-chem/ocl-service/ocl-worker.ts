import {MolNotationType, OCLServiceCall} from './consts';
import {CHEM_PROP_MAP} from './calculations';
import * as OCL from 'openchemlib/full';

type OCLWorkerReturnType = {res: {[key: string]: Array<number>}, errors: string[]};

onmessage = ({data: {op, data, argList, notationType}}) => {
  switch (op) {
  case OCLServiceCall.CHEM_PROPERTIES:
    postMessage(getChemProperties(data, argList, notationType));
    break;
  case OCLServiceCall.DRUG_LIKENESS:
    postMessage(getDrugLikeliness(data, notationType));
  default:
    postMessage(getToxRisks(data, argList, notationType));
    break;
  }
};

function isSmiles(notationType: string): boolean {
  return notationType === MolNotationType.SMILES;
}

function getChemProperties(molList: Array<string>, propList: string[], notationType: string): OCLWorkerReturnType {
  const res: {[key: string]: Array<number>} = {};
  const errors: string[] = [];
  for (const propName of propList)
    res[propName] = new Array(molList.length).fill(null);
  molList.forEach((smiles, i) => {
    try {
      const mol = isSmiles(notationType) ? OCL.Molecule.fromSmiles(smiles) : OCL.Molecule.fromMolfile(smiles);
      propList.forEach((p) => {
        res[p][i] = CHEM_PROP_MAP[p].valueFunc(mol);
      });
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

function getToxRisks(molList: Array<string>, riskTypes: number[], notationType: string): OCLWorkerReturnType {
  const res: {[key: string]: Array<number>} = {};
  const errors: string[] = [];
  for (const propName of riskTypes)
    res[propName] = new Array(molList.length).fill(0);
  const toxicityPredictor = new OCL.ToxicityPredictor();
  molList.forEach((smiles, i) => {
    try {
      const mol = isSmiles(notationType) ? OCL.Molecule.fromSmiles(smiles) : OCL.Molecule.fromMolfile(smiles);
      riskTypes.forEach((p) => {
        res[p][i] = toxicityPredictor.assessRisk(mol, p);
      });
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

function getDrugLikeliness(molList: Array<string>, notationType: string): OCLWorkerReturnType {
  const colName = 'Drug likeness score';
  const res: {[key: string]: Array<number>} = {[colName]: new Array(molList.length).fill(0)};
  const errors: string[] = [];
  const dlp = new OCL.DruglikenessPredictor();
  molList.forEach((smiles, i) => {
    try {
      const mol = isSmiles(notationType) ? OCL.Molecule.fromSmiles(smiles) : OCL.Molecule.fromMolfile(smiles);
      const score = dlp.assessDruglikeness(mol);
      res[colName][i] = score;
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

