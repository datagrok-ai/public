import {MolNotationType, OCLServiceCall} from '../consts';
import * as OCL from 'openchemlib/full';

type OCLWorkerReturnType = {res: {[key: string]: Array<number>}, errors: string[]};

onmessage = ({data: {op, data, argList, notationType}}) => {
  switch (op) {
  case OCLServiceCall.CHEM_PROPERTIES:
    postMessage(getChemProperties(data, argList, notationType));
    break;
  case OCLServiceCall.DRUG_LIKENESS:
    postMessage(getDrugLikeliness(data, notationType));
  case OCLServiceCall.Molfile_To_V3K:
    postMessage(molfilesToV3K(data));
  case OCLServiceCall.RECALCULATE_COORDINATES:
    postMessage(recalculateCoordinates(data, notationType));
    break;
  default:
    postMessage(getToxRisks(data, argList, notationType));
    break;
  }
};

function isSmiles(notationType: string): boolean {
  return notationType === MolNotationType.SMILES;
}

function oclMol(mol: string, notationType: string): OCL.Molecule {
  if (isSmiles(notationType) && mol.length > 5000)
    throw new Error('Invalid molecule string');
  return isSmiles(notationType) ? OCL.Molecule.fromSmiles(mol) : OCL.Molecule.fromMolfile(mol);
}

function recalculateCoordinates(molList: Array<string>, notationType: string) {
  const res: Array<string> = new Array(molList?.length ?? 0).fill('');
  const errors: string[] = [];
  molList.forEach((molfile, i) => {
    if (!molfile) {
      res[i] = '';
      return;
    }
    try {
      const mol = oclMol(molfile, notationType);
      mol.inventCoordinates();
      const v3 = mol.toMolfileV3();
      res[i] = v3.replace('STERAC1', 'STEABS');
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
      res[i] = '';
    }
  });
  return {res, errors};
}

function molfilesToV3K(molList: Array<string>) {
  const res: Array<string> = new Array(molList?.length ?? 0).fill('');
  const errors: string[] = [];
  molList.forEach((molfile, i) => {
    if (!molfile) {
      res[i] = '';
      return;
    }
    try {
      const mol = OCL.Molecule.fromMolfile(molfile);
      const v3 = mol.toMolfileV3();
      res[i] = v3.replace('STERAC1', 'STEABS');
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
      res[i] = '';
    }
  });
  return {res, errors};
}

function getChemProperties(molList: Array<string>, propList: string[], notationType: string): OCLWorkerReturnType {
  const res: {[key: string]: Array<number>} = {};
  const errors: string[] = [];
  for (const propName of propList)
    res[propName] = new Array(molList.length).fill(null);
  molList.forEach((smiles, i) => {
    try {
      const mol = oclMol(smiles, notationType);
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
      const mol = oclMol(smiles, notationType);
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
      const mol = oclMol(smiles, notationType);
      const score = dlp.assessDruglikeness(mol);
      res[colName][i] = score;
    } catch (e) {
      errors.push(e instanceof Error ? e.message : e as string);
    }
  });
  return {res, errors};
}

// duplication required for webpack to separate the worker OCL context from main thread context
const CHEM_PROP_MAP: {[k: string]: IChemProperty} = {
  'MW': {name: 'MW', type: 'float', valueFunc: (m: OCL.Molecule) => m.getMolecularFormula().absoluteWeight},
  'HBA': {name: 'HBA', type: 'int', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).acceptorCount},
  'HBD': {name: 'HBD', type: 'int', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).donorCount},
  'LogP': {name: 'LogP', type: 'float', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logP},
  'LogS': {name: 'LogS', type: 'float', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logS},
  'PSA': {
    name: 'PSA', type: 'float', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).polarSurfaceArea},
  'Rotatable bonds': {name: 'Rotatable bonds', type: 'int',
    valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).rotatableBondCount},
  'Stereo centers': {name: 'Stereo centers', type: 'int',
    valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).stereoCenterCount},
  'Molecule charge': {name: 'Molecule charge', type: 'int', valueFunc: (m: OCL.Molecule) => getMoleculeCharge(m)},
};

function getMoleculeCharge(mol: OCL.Molecule): number {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx)
    moleculeCharge += mol.getAtomCharge(atomIndx);

  return moleculeCharge;
}
interface IChemProperty {
  name: string;
  valueFunc: (mol: OCL.Molecule) => any;
  type: IChemPropertyType;
}
// we need this instead of DG.TYPE.FLOAT or int, because anything used by workers can not import DG, as it is external
type IChemPropertyType = 'float' | 'int';
