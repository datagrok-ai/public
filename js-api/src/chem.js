/**
 * Cheminformatics support
 * @module chem
 * */

/**
 * @typedef {string} SimilarityMetric
 **/

import {BitSet, DataFrame} from "./dataframe";
import {SIMILARITY_METRIC} from "./const";

/** Cheminformatics-related routines */

/**
 * Returns molecules similar to the reference one.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/similarity-search}
 * @async
 * @param {Column} column - Molecule column to search in.
 * @param {string} molecule - Reference molecule in SMILES format.
 * @param {SimilarityMetric} metric - Metric to use.
 * @param {number} limit - Maximum number of results to return.
 * @param {number} minScore - Minimum similarity score for a molecule to be included.
 * @returns {Promise<DataFrame>}
 * */
export function similaritySearch(column, molecule, metric = SIMILARITY_METRIC.TANIMOTO, limit = 10, minScore = 0.7) {
    return new Promise((resolve, reject) => grok_Chem_SimilaritySearch(column.d, molecule, metric,
        limit, minScore, (t) => resolve(new DataFrame(t))));
}

/**
 * Returns the specified number of most diverse molecules in the column.
 * See example: {@link https://datagrok.ai/help/domains/chem/diversity-search}
 * @async
 * @param {Column} column - Column with molecules in which to search.
 * @param {SimilarityMetric} metric - Metric to use.
 * @param {number} limit - Number of molecules to return.
 * @returns {Promise<DataFrame>}
 * */
export function diversitySearch(column, metric = SIMILARITY_METRIC.TANIMOTO, limit = 10) {
    return new Promise((resolve, reject) => grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols)));
}

/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 * See example: {@link substructure-search}
 * @async
 * @param {Column} column - Column with molecules to search.
 * @param {string} pattern - Pattern, either SMARTS or SMILES.
 * @param {boolean} isSmarts - Whether the pattern is SMARTS.
 * @returns {Promise<BitSet>}
 * */
export function substructureSearch(column, pattern, isSmarts = true) {
    return new Promise((resolve, reject) => grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs) => resolve(new BitSet(bs))));
}

/**
 * Performs R-group analysis.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 * @async
 * @param {DataFrame} table - Table.
 * @param {string} column - Column name with SMILES to analyze.
 * @param {string} core - Core in the SMILES format.
 * @returns {Promise<DataFrame>}
 * */
export function rGroup(table, column, core) {
    return new Promise((resolve, reject) => grok_Chem_RGroup(table.d, column, core, () => resolve(table)));
}

/**
 * Finds Most Common Substructure in the specified column.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
 * @async
 * @param {Column} column - Column with SMILES to analyze.
 * @returns {Promise<string>}
 * */
export function mcs(column) {
    return new Promise((resolve, reject) => grok_Chem_MCS(column.d, (mcs) => resolve(mcs)));
}

/**
 * Calculates specified descriptors for the molecular column.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 *
 * @async
 * @param {DataFrame} table - Table.
 * @param {string} column - Column name with SMILES to calculate descriptors for.
 * @param {string[]} descriptors - RDKit descriptors to calculate.
 * @returns {Promise<DataFrame>}
 * */
export function descriptors(table, column, descriptors) {
    return new Promise((resolve, reject) => grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table)));
}

/**
 * Returns available descriptors tree.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 *
 * @async
 * @returns {Promise<Object>}
 * */
export function descriptorsTree() {
    return new Promise((resolve, reject) => grok_Chem_DescriptorsTree((tree) => resolve(JSON.parse(tree))));
}

/**
 * Renders a molecule to SVG
 * https://public.datagrok.ai/js/samples/domains/chem/mol-rendering
 * @param {string} smiles
 * @param {number} width
 * @param {number} height
 * @returns {HTMLDivElement}
 * */
export function svgMol(smiles, width = 300, height = 200) {
    let root = document.createElement('div');
    import('openchemlib/full.js').then((OCL) => {
        let m = OCL.Molecule.fromSmiles(smiles);
        root.innerHTML = m.toSVG(width, height);
    });
    return root;
}

/**
 * Molecule sketcher.
 * @param {function} handler Molecule on change handler, SMILES
 * @param {string} smiles Initial molecule
 * @returns {HTMLElement}
 * */
export function sketcher(handler, smiles = '') { return grok_Chem_Sketcher(handler, smiles); }
