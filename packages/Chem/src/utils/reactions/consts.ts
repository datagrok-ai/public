/* eslint-disable max-len */
import {NamedReaction} from './types';

export const UNCATEGORIZED_REACTION_NAME = 'Uncategorized';

/** Canonical SMILES fragments for salts, water and common byproducts */
export const SALT_AND_WATER_FRAGMENTS = [
  'O', '[H]O[H]', 'O[H]', '[H]O', '[H][O][H]',
  '[Cl][H]', '[H][Cl]', 'Cl[H]', 'Cl', '[H]Cl',
  '[Br][H]', '[H][Br]', 'Br[H]', 'Br', '[H]Br',
  '[I][H]', '[H][I]', 'I[H]', 'I', '[H]I',
  '[F][H]', '[H][F]', 'F[H]', 'F', '[H]F',
  'O=S(=O)(O)O', 'S(=O)(=O)(O)(O)',
];

/** When true, use the RDKitReactionRenderer class (offscreen canvas + LRU cache) for reaction previews.
 *  When false, use simple direct RDKit rendering: get_rxn() → draw_to_canvas() → delete(). */
export const USE_RDKIT_REACTION_RENDERER = true;

/** Storage paths */
export const REACTIONS_STORAGE_PATH = 'System:AppData/Chem/reactions';
export const DEFAULT_REACTIONS_FILE = 'reactions/default-reactions.json';
export const USER_REACTIONS_SUFFIX = '-reactions.json';

/** Get all unique categories from a list of reactions. */
export function getReactionCategories(reactions: NamedReaction[]): string[] {
  const categories = new Set<string>();
  for (const r of reactions)
    categories.add(r.category);
  return Array.from(categories).sort();
}
