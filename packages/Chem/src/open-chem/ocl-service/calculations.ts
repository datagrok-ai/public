import * as OCL from 'openchemlib/full';
import {IChemProperty} from './types';

export const CHEM_PROP_MAP: {[k: string]: IChemProperty} = {
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

export function getMoleculeCharge(mol: OCL.Molecule): number {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx)
    moleculeCharge += mol.getAtomCharge(atomIndx);

  return moleculeCharge;
}
