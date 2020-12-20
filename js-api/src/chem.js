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
 * @param {Column} column - Molecule column to search in
 * @param {string} pattern - Reference molecule in one of formats supported by RDKit:
 *     smiles, cxsmiles, molblock, v3Kmolblock
 * @param {boolean} settings.sorted -
 *     if set, returns a two-column dataframe with molecule strings and scores,
 *     sorted in descending order by the score
 * @returns {Promise<DataFrame>, if sorted; Promise<Column>, otherwise}
 * */
export async function similarityScoring(column, pattern = null, settings = { sorted: false }) {
    
    let foo = await grok.functions.eval('Chem:similarityScoring');
    let call = await foo.prepare({
        'molStringsColumn': column,
        'molString': pattern == null ? "" : pattern,
        'sorted': settings.sorted
    });
    await call.call();
    if (pattern != null) {
        let result = call.getParamValue('result');
        return settings.sorted ? result : result.columns.byIndex(0);
    }
    
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
    return new Promise((resolve, reject) => grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols), (e) => reject(e)));
}

/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 * See example: {@link substructure-search}
 * @async
 * @param {Column} column - Column with molecules to search
 * @param {string} pattern - Pattern, either one of which RDKit supports
 * @returns {Promise<BitSet>}
 * */
export async function substructureSearch(column, pattern = null, settings = null) {
    
    let foo = await grok.functions.eval('Chem:substructureSearch');
    let call = await foo.prepare({
        'molStringsColumn': column,
        'molString': pattern == null ? "" : pattern,
        'substructLibrary':
          (settings && settings.hasOwnProperty('substructLibrary') && settings.substructLibrary === false) ?
            false : true
    });
    await call.call();
    // unpacking our BitSet object from a synthetic column
    return call.getParamValue('result').get(0);
    
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
    return new Promise((resolve, reject) => grok_Chem_RGroup(table.d, column, core, () => resolve(table), (e) => reject(e)));
}

/**
 * Finds Most Common Substructure in the specified column.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
 * @async
 * @param {Column} column - Column with SMILES to analyze.
 * @returns {Promise<string>}
 * */
export function mcs(column) {
    return new Promise((resolve, reject) => grok_Chem_MCS(column.d, (mcs) => resolve(mcs), (e) => reject(e)));
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
    return new Promise((resolve, reject) => grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table), (e) => reject(e)));
}

/**
 * Returns available descriptors tree.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 *
 * @async
 * @returns {Promise<Object>}
 * */
export function descriptorsTree() {
    return new Promise((resolve, reject) => grok_Chem_DescriptorsTree((tree) => resolve(JSON.parse(tree)), (e) => reject(e)));
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
 * Sketches Molecule sketcher.
 * @param {function} onChangedCallback - a function that accepts (smiles, molfile)
 * @param {string} smiles Initial molecule
 * @returns {HTMLElement}
 * */
export function sketcher(onChangedCallback, smiles = '') { return grok_Chem_Sketcher(onChangedCallback, smiles); }
