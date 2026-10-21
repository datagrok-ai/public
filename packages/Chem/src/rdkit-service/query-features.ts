import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/** RDKit's descriptions of the queries that read what only sanitization perceives in a molecule: ring sizes (r<n>, of
 * atoms and bonds, and their bounds) and bridgeheads, which read RDKit's rings, hybridization (^), radicals. */
const SANITIZED_PROPERTY = /RingSize|IsBridgehead|Hybridization|NumRadicalElectrons/;
/** A ring count (R<n>, of atoms and bonds; `less_`, `greater_` and `range_` for its bounds). */
const RING_COUNT = /InNRings$/;

interface QueryNode {descr?: string, val?: number, children?: QueryNode[], subquery?: JsonMolecule}
interface JsonMolecule {
  atoms?: {nRad?: number}[],
  extensions?: {name: string, atomQueries?: QueryNode[], bondQueries?: QueryNode[]}[],
}

function asksForSanitized(node: QueryNode, defaultRadicals: number): boolean {
  const descr = node.descr ?? '';
  if (SANITIZED_PROPERTY.test(descr))
    return true;
  // R<n> for n > 0 and its bounds count rings; R (-1) and R0 ask for ring membership only, as do r, x and @, which the
  // fast read's rings give as sanitization does
  if (RING_COUNT.test(descr) && (descr.includes('_') || (node.val ?? 0) > 0))
    return true;
  return (node.children ?? []).some((child) => asksForSanitized(child, defaultRadicals)) ||
    (!!node.subquery && moleculeAsksForSanitized(node.subquery, defaultRadicals));
}

function moleculeAsksForSanitized(mol: JsonMolecule, defaultRadicals: number): boolean {
  if ((mol.atoms ?? []).some((atom) => (atom.nRad ?? defaultRadicals) > 0))
    return true;
  const queries = mol.extensions?.find((e) => e.name === 'rdkitQueries');
  return [...queries?.atomQueries ?? [], ...queries?.bondQueries ?? []]
    .some((node) => asksForSanitized(node, defaultRadicals));
}

/**
 * Whether matching `query` needs the molecules as RDKit sanitizes them, not as Chem's RDKit search reads its candidates
 * for speed (as written, unsanitized: rings from a depth-first search, no hybridization, no radicals). Ring counts
 * (R<n>) and ring sizes (r<n>) count RDKit's rings (the symmetrized smallest set of smallest rings), which only
 * sanitization finds; hybridization (^) and radicals only sanitization sets. Read from the query molecule (its RDKit
 * JSON: the atoms and their queries, recursive SMARTS included); a query whose description can't be read needs them.
 */
export function queryNeedsSanitizedMolecules(query: RDMol): boolean {
  try {
    const json = JSON.parse(query.get_json());
    const defaultRadicals: number = json.defaults?.atom?.nRad ?? 0;
    return (json.molecules ?? []).some((mol: JsonMolecule) => moleculeAsksForSanitized(mol, defaultRadicals));
  } catch {
    return true;
  }
}
