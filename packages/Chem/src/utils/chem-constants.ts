/**
 * Leaf module for primitive chemistry constants shared across the package — including
 * web-worker code. Keep it import-free (NO `datagrok-api` imports) so it stays safe to
 * bundle into worker entry points.
 */

/**
 * Maximum SMILES string length the package will attempt to parse. Longer SMILES are
 * skipped — parsing very long SMILES is slow and can crash RDKit/OCL. Molblocks
 * (which contain newlines) are exempt from this cap.
 */
export const MAX_SMILES_LENGTH = 5000;
