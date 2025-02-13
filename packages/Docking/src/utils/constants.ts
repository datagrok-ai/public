import * as DG from 'datagrok-api/dg';
import { DockingPackage } from '../package-utils';

export const _package = new DockingPackage();
export const TARGET_PATH = 'System:AppData/Docking/targets';
export const CACHED_RESULTS: DG.LruCache<string, DG.DataFrame> = new DG.LruCache<string, DG.DataFrame>();
export let BINDING_ENERGY_COL = 'binding energy';
export let BINDING_ENERGY_COL_UNUSED = '';
export let POSE_COL = 'pose';
export let POSE_COL_UNUSED = '';
export const ERROR_COL_NAME = 'error';
export const ERROR_MESSAGE = 'Error log';

export function setPose(value: string) {
  POSE_COL_UNUSED = value;
}

export function setAffinity(value: string) {
  BINDING_ENERGY_COL_UNUSED = value;
}

export const AUTODOCK_PROPERTY_DESCRIPTIONS: {[colName: string]: string} = {
  'intermolecular (1)': 'Final Intermolecular Energy',
  'electrostatic': 'Electrostatic Energy',
  'ligand fixed': 'Moving Ligand-Fixed Receptor',
  'ligand moving': 'Moving Ligand-Moving Receptor',
  'total internal (2)': 'Final Total Internal Energy',
  'torsional free (3)': 'Torsional Free Energy',
  'unbound systems (4)': 'Unbound System\s Energy' 
}

export const BOLTZ_CONFIG_PATH = 'System:AppData/Chem/boltz';

export const BOLTZ_PROPERTY_DESCRIPTIONS: { [colName: string]: string } = {
  'confidence_score': 'Overall prediction quality score',
  'ptm': 'Global fold similarity measure',
  'iptm': 'Accuracy of chain interactions',
  'ligand_iptm': 'Confidence in ligand binding',
  'protein_iptm': 'Confidence in protein interactions',
  'complex_plddt': 'Average per-residue confidence score',
  'complex_iplddt': 'Confidence in chain interfaces',
  'complex_pde': 'Uncertainty in chain positioning',
  'complex_ipde': 'Uncertainty in interface docking',
  'chains_ptm': 'Confidence per individual chain',
  'pair_chains_iptm': 'Interaction accuracy between chains'
};
  
type SequenceModification = {
  position: number;  // Index of the residue, starting from 1
   ccd: string;       // CCD code of the modified residue
};
  
type ProteinEntity = {
  id: string | string[];  // Chain ID or multiple Chain IDs for identical entities
  sequence: string;       // Sequence (only for protein, dna, rna)
  msa?: string;            // MSA Path (only for protein)
  modifications?: SequenceModification[];  // List of modifications
};
  
type LigandEntity = {
  id: string | string[];  // Chain ID or multiple Chain IDs for identical entities
  smiles: string;         // SMILES (only for ligands)
  ccd: string;            // CCD (only for ligands)
};
  
type SequenceEntity = {
  protein?: ProteinEntity;  // Protein entity with sequence and msa
  ligand?: LigandEntity;    // Ligand entity with smiles and ccd
};
  
export type Config = {
  version: number;
  sequences: SequenceEntity[];  // Array of sequence entities, each with an entity type
  constraints?: Constraint[];   // Optional constraints (can be empty)
};
  
type Constraint = {
  bond?: {
    atom1: [string, number, string];  // [CHAIN_ID, RES_IDX, ATOM_NAME]
    atom2: [string, number, string];  // [CHAIN_ID, RES_IDX, ATOM_NAME]
  };
  pocket?: {
    binder: string;                    // CHAIN_ID for the binder
    contacts: [string, number][];      // List of [CHAIN_ID, RES_IDX] for contacts
  };
};