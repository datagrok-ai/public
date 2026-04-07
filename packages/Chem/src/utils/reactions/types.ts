/** Reaction mode: transform one column or combine two columns */
export type ReactionMode = 'transformation' | 'two-component';

/**
 * A single named reaction definition.
 * Variable regions can be marked with ${variableName} inside the reaction SMARTS string.
 */
export type NamedReaction = {
  /** Unique identifier (slug or uuid) */
  id: string;
  /** Display name, e.g. "Suzuki Coupling" */
  name: string;
  /** Reaction SMIRKS/SMARTS string */
  reactionSmarts: string;
  /** Whether reaction operates on one column or combines two */
  mode: ReactionMode;
  /** Grouping category, e.g. "C-C Bond Formation", "Deprotection" */
  category: string;
  /** Human-readable description (shown in tooltips) */
  description?: string;
  /** Literature references (DOIs or URLs) */
  references?: string[];
  /** User-tunable variable placeholders used in the SMARTS string */
  variables?: {
    [key: string]: {
      name: string;
      type: 'string' | 'int' | 'float';
      defaultValue: any;
      description?: string;
    }
  };
  /** Example reactant SMILES for preview */
  exampleReactants?: string[];
  /** Expected product SMILES for preview */
  exampleProducts?: string[];
  /** Author who created this reaction */
  author?: string;
  /** true = saved by user, false = shipped default */
  isUserDefined?: boolean;
  /** Searchable tags */
  tags?: string[];
}

/** Result of running a reaction on a list of molecules */
export type ReactionResult = {
  /** SMILES array of products (same length as input) */
  products: string[];
  /** Number of molecules that caused RDKit errors */
  errors: number;
  /** Number of molecules where reaction produced no products */
  noProducts: number;
  /** Number of empty/null input molecules */
  emptyMolecules: number;
  /** The resolved SMARTS that was actually executed */
  reactionSmarts: string;
}

/** Result of a two-component reaction */
export type TwoComponentReactionResult = {
  /** Reactant 1 SMILES for each pair */
  reactants1: string[];
  /** Reactant 2 SMILES for each pair */
  reactants2: string[];
  /** Product SMILES for each pair */
  products: string[];
  /** Number of pairs that caused errors */
  errors: number;
  /** Number of pairs that produced no products */
  noProducts: number;
}

/** @deprecated Use NamedReaction instead */
export type CommonTransformationReaction = {
  name: string;
  reactionSmirks: string;
  category?: string;
  variables?: {
    [key: string]: {
      name: string;
      type: 'string' | 'int' | 'float';
      defaultValue: any;
    }
  };
}
