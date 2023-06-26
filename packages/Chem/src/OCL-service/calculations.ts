import * as OCL from 'openchemlib/full';
import {IChemProperty} from './types';

export const CHEM_PROP_MAP: {[k: string]: IChemProperty} = {
  'MW': {name: 'MW', valueFunc: (m: OCL.Molecule) => m.getMolecularFormula().absoluteWeight},
  'HBA': {name: 'HBA', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).acceptorCount},
  'HBD': {name: 'HBD', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).donorCount},
  'LogP': {name: 'LogP', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logP},
  'LogS': {name: 'LogS', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logS},
  'PSA': {
    name: 'PSA', valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).polarSurfaceArea},
  'Rotatable bonds': {name: 'Rotatable bonds',
    valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).rotatableBondCount},
  'Stereo centers': {name: 'Stereo centers',
    valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).stereoCenterCount},
  'Molecule charge': {name: 'Molecule charge', valueFunc: (m: OCL.Molecule) => getMoleculeCharge(m)},
};

export function getMoleculeCharge(mol: OCL.Molecule): number {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx)
    moleculeCharge += mol.getAtomCharge(atomIndx);

  return moleculeCharge;
}
