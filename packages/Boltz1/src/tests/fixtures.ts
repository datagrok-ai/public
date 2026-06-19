import type {StructureAndBindingStartParams} from 'boltz-api/resources/predictions/structure-and-binding';
import type {AdmeStartParams} from 'boltz-api/resources/predictions/adme';
import type {DesignStartParams as SmallMoleculeDesignStartParams} from 'boltz-api/resources/small-molecule/design';
import type {
  LibraryScreenStartParams as SmallMoleculeLibraryScreenStartParams,
} from 'boltz-api/resources/small-molecule/library-screen';
import type {DesignStartParams as ProteinDesignStartParams} from 'boltz-api/resources/protein/design';
import type {
  LibraryScreenStartParams as ProteinLibraryScreenStartParams,
} from 'boltz-api/resources/protein/library-screen';

import {ADME_MODEL, STRUCTURE_AND_BINDING_MODEL} from '../utils/boltz-api-constants';

export const PROTEIN_SEQUENCE = 'MKTIIALSYIFCLVFA';
export const BINDER_SEQUENCE = 'MKTAYIVKSHFSRQ';

export const ASPIRIN_SMILES = 'CC(=O)OC1=CC=CC=C1C(=O)O';
export const IBUPROFEN_SMILES = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';
export const CAFFEINE_SMILES = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C';
export const PHENOL_SMILES = 'C1=CC=C(C=C1)O';
export const TOLUENE_SMILES = 'CC1=CC=CC=C1';

// Shared target for the small-molecule design and library-screen examples in the docs.
const smallMoleculeTarget: SmallMoleculeDesignStartParams['target'] = {
  entities: [{type: 'protein', value: PROTEIN_SEQUENCE, chain_ids: ['A']}],
  pocket_residues: {A: [2, 3, 4, 7, 8, 9]},
  reference_ligands: ['CC(=O)Oc1ccccc1C(=O)O'],
  constraints: [
    {type: 'pocket', binder_chain_id: 'L', contact_residues: {A: [2, 3, 4, 7, 8, 9]}, max_distance_angstrom: 6.0},
  ],
};

// Full filter set from the docs, used by design (the generator produces molecules that satisfy it).
const moleculeFilters: SmallMoleculeDesignStartParams['molecule_filters'] = {
  boltz_smarts_catalog_filter_level: 'recommended',
  custom_filters: [
    {type: 'lipinski_filter', max_mw: 500, max_logp: 5, max_hbd: 5, max_hba: 10, allow_single_violation: false},
    {
      type: 'rdkit_descriptor_filter',
      mol_wt: {min: 150, max: 500}, mol_logp: {max: 5}, tpsa: {max: 140},
      num_h_donors: {max: 5}, num_h_acceptors: {max: 10}, num_rotatable_bonds: {max: 10},
      num_heteroatoms: {max: 12}, num_aromatic_rings: {min: 1, max: 4}, num_rings: {max: 6},
      fraction_csp3: {min: 0.2},
    },
    {type: 'smarts_custom_filter', patterns: ['[N+](=O)[O-]', 'C(=O)Cl']},
    {type: 'smarts_catalog_filter', catalog: 'PAINS'},
    {type: 'smiles_regex_filter', patterns: ['P', 'S(=O)(=O)Cl']},
  ],
};

// Lenient filters for screening user-supplied molecules: exercise the filter path without rejecting small ligands.
const screenMoleculeFilters: SmallMoleculeLibraryScreenStartParams['molecule_filters'] = {
  boltz_smarts_catalog_filter_level: 'disabled',
  custom_filters: [
    {type: 'lipinski_filter', max_mw: 900, max_logp: 8, max_hbd: 10, max_hba: 20, allow_single_violation: true},
  ],
};

export const structureAndBindingRequest: StructureAndBindingStartParams = {
  model: STRUCTURE_AND_BINDING_MODEL,
  input: {
    entities: [
      {type: 'protein', value: PROTEIN_SEQUENCE, chain_ids: ['A']},
      {type: 'ligand_ccd', value: 'AIN', chain_ids: ['B']},
    ],
    binding: {type: 'ligand_protein_binding', binder_chain_id: 'B'},
    constraints: [
      {type: 'pocket', binder_chain_id: 'B', contact_residues: {A: [10, 11, 12]},
        max_distance_angstrom: 6.0, force: false},
      {type: 'contact',
        token1: {type: 'polymer_contact', chain_id: 'A', residue_index: 3},
        token2: {type: 'polymer_contact', chain_id: 'A', residue_index: 11},
        max_distance_angstrom: 8.0, force: false},
    ],
    bonds: [
      {atom1: {type: 'polymer_atom', chain_id: 'A', residue_index: 11, atom_name: 'SG'},
        atom2: {type: 'polymer_atom', chain_id: 'A', residue_index: 1, atom_name: 'NZ'}},
    ],
    model_options: {recycling_steps: 3, sampling_steps: 200, step_scale: 1.638},
    templates: [
      {template_structure: {type: 'url', url: 'https://files.rcsb.org/download/1CRN.cif'},
        template_chains: [{input_chain_id: 'A', template_chain_id: 'A'}],
        force_threshold_angstroms: 5.0},
    ],
    num_samples: 3,
  },
};

export const admeRequest: AdmeStartParams = {
  model: ADME_MODEL,
  input: {
    molecules: [
      {smiles: ASPIRIN_SMILES, id: 'aspirin'},
      {smiles: IBUPROFEN_SMILES, id: 'ibuprofen'},
      {smiles: CAFFEINE_SMILES},
    ],
  },
};

export const smallMoleculeDesignRequest: SmallMoleculeDesignStartParams = {
  target: smallMoleculeTarget,
  chemical_space: 'enamine_real',
  molecule_filters: moleculeFilters,
  num_molecules: 100,
};

export const smallMoleculeLibraryScreenRequest: SmallMoleculeLibraryScreenStartParams = {
  target: smallMoleculeTarget,
  molecules: [
    {smiles: ASPIRIN_SMILES, id: 'aspirin'},
    {smiles: PHENOL_SMILES, id: 'phenol'},
    {smiles: TOLUENE_SMILES},
  ],
  molecule_filters: screenMoleculeFilters,
};

export const proteinDesignRequest: ProteinDesignStartParams = {
  target: {
    type: 'no_template',
    entities: [
      {type: 'protein', value: PROTEIN_SEQUENCE, chain_ids: ['A']},
      {type: 'ligand_ccd', value: 'ATP', chain_ids: ['L1']},
    ],
    epitope_residues: {A: [13, 14, 15]},
    non_binding_residues: {A: [0, 1, 2]},
    epitope_ligand_chains: ['L1'],
  },
  binder_specification: {
    type: 'no_template',
    modality: 'custom_protein',
    entities: [{type: 'designed_protein', chain_ids: ['B'], value: 'MKTAYI10..20VKSHFSRQ'}],
    rules: {
      excluded_amino_acids: ['C'],
      max_hydrophobic_fraction: 0.5,
      excluded_sequence_motifs: ['NXS'],
    },
  },
  num_proteins: 10,
};

export const proteinLibraryScreenRequest: ProteinLibraryScreenStartParams = {
  target: {
    type: 'no_template',
    entities: [
      {type: 'protein', value: PROTEIN_SEQUENCE, chain_ids: ['A']},
      {type: 'ligand_ccd', value: 'ATP', chain_ids: ['L1']},
    ],
    epitope_residues: {A: [10, 11, 12]},
    non_binding_residues: {A: [0, 1, 2]},
    epitope_ligand_chains: ['L1'],
  },
  proteins: [
    {entities: [{type: 'protein', value: BINDER_SEQUENCE, chain_ids: ['B']}], id: 'binder-001'},
    {entities: [{type: 'protein', value: 'ACDEFGHIKLMNPQRSTVWY', chain_ids: ['B']}], id: 'binder-002'},
  ],
};
