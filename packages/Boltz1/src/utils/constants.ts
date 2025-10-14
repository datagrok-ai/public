import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const BOLTZ_CONFIG_PATH = 'System:AppData/Boltz1/boltz';

export const CACHED_RESULTS: DG.LruCache<string, DG.DataFrame> = new DG.LruCache<string, DG.DataFrame>();

export const BOLTZ_PROPERTY_DESCRIPTIONS: { [colName: string]: string } = {
  'confidence_score': 'Aggregated score used to sort the predictions, corresponds to 0.8 * complex_plddt + 0.2 * iptm (ptm for single chains)',
  'ptm': 'Predicted TM score for the complex',
  'iptm': 'Predicted TM score when aggregating at the interfaces',
  'ligand_iptm': 'ipTM but only aggregating at protein-ligand interfaces',
  'protein_iptm': 'ipTM but only aggregating at protein-protein interfaces',
  'complex_plddt': 'Average pLDDT score for the complex',
  'complex_iplddt': 'Average pLDDT score when upweighting interface tokens',
  'complex_pde': 'Average PDE score for the complex',
  'complex_ipde': 'Average PDE score when aggregating at interfaces ',
  'chains_ptm': 'Predicted (interface) TM score between each pair of chains',
  'pair_chains_iptm': '# Predicted (interface) TM score between each pair of chains'
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

export interface BoltzResponse {
  success: boolean;
  error: string | null;
  result: string | null;
}