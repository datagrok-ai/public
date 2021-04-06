/**
 * Cheminformatics support
 * @module chem
 * */
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
  function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
  return new (P || (P = Promise))(function (resolve, reject) {
    function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
    function rejected(value) { try { step(generator['throw'](value)); } catch (e) { reject(e); } }
    function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
    step((generator = generator.apply(thisArg, _arguments || [])).next());
  });
};
import { BitSet, DataFrame } from './dataframe';
import { SIMILARITY_METRIC } from './const';
let api = window;
/** Cheminformatics-related routines */
/**
 * Returns molecules similar to the reference one.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/similarity-search}
 * @async
 * @deprecated
 * @param {Column} column - Molecule column to search in
 * @param {string} molecule - Reference molecule in one of formats supported by RDKit:
 *     smiles, cxsmiles, molblock, v3Kmolblock
 * @param {Object} settings
 * @param {boolean} settings.sorted -
 *     if set, returns a two-column dataframe with molecule strings and scores,
 *     sorted in descending order by the score
 * @returns {Promise<DataFrame>, if sorted; Promise<Column>, otherwise}
 * */
export function similarityScoring(column, molecule = '', settings = { sorted: false }) {
  return __awaiter(this, void 0, void 0, function* () {
    let foo = yield grok.functions.eval('Chem:similarityScoring');
    let call = yield foo.prepare({
      'molStringsColumn': column,
      'molString': molecule,
      'sorted': settings.sorted
    });
    yield call.call();
    if (molecule.length != 0) {
      let result = call.getParamValue('result');
      return settings.sorted ? result : result.columns.byIndex(0);
    }
  });
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
  return new Promise((resolve, reject) => api.grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols), (e) => reject(e)));
}
/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 * See example: {@link substructure-search}
 * @async
 * @deprecated
 * @param {Column} column - Column with molecules to search
 * @param {string} pattern - Pattern, either one of which RDKit supports
 * @param settings
 * @returns {Promise<BitSet>}
 * */
export function substructureSearch(column, pattern = '', settings) {
  return __awaiter(this, void 0, void 0, function* () {
    let foo = yield grok.functions.eval('Chem:substructureSearch');
    let call = yield foo.prepare({
      'molStringsColumn': column,
      'molString': pattern,
      'substructLibrary': !(settings.hasOwnProperty('substructLibrary') && !settings.substructLibrary)
    });
    yield call.call();
    // unpacking our BitSet object from a synthetic column
    return call.getParamValue('result').get(0);
  });
}
/**
 * Computes similarity scores for molecules in the input vector based on a preferred similarity score.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/similarity-scoring-scores}
 * @async
 * @param {Column} column - Column with molecules to search in
 * @param {string} molecule - Reference molecule in one of formats supported by RDKit:
 *   smiles, cxsmiles, molblock, v3Kmolblock, and inchi
 * @param {Object} settings - Properties for the similarity function (type, parameters, etc.)
 * @returns {Promise<Column>} - Column of corresponding similarity scores
 * */
export function getSimilarities(column, molecule = '', settings = {}) {
  return __awaiter(this, void 0, void 0, function* () {
    let foo = yield grok.functions.eval('Chem:getSimilarities');
    let call = yield foo.prepare({
      'molStringsColumn': column,
      'molString': molecule
    });
    yield call.call();
    // TODO: figure out what's the state in returning columns from package functions
    return (molecule.length != 0) ? call.getParamValue('result').columns.byIndex(0) : null;
  });
}
/**
 * Computes similarity scores for molecules in the input vector based on a preferred similarity score.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/similarity-scoring-sorted}
 * @async
 * @param {Column} column - Column with molecules to search in
 * @param {string} molecule - Reference molecule in one of formats supported by RDKit:
 *   smiles, cxsmiles, molblock, v3Kmolblock, and inchi
 * @param {Object} settings - Properties for the similarity function
 * @param {int} limit - Would return top limit molecules based on the score
 * @param {int} cutoff - Would drop molecules which score is lower than cutoff
 * @returns {Promise<DataFrame>} - DataFrame with 3 columns:
 *   - molecule: original molecules string representation from the input column
 *   - score: similarity scores within the range from 0.0 to 1.0;
 *            DataFrame is sorted descending by this column
 *   - index: indices of the molecules in the original input column
 * */
export function findSimilar(column, molecule = '', settings = { limit: Number.MAX_VALUE, cutoff: 0.0 }) {
  return __awaiter(this, void 0, void 0, function* () {
    let foo = yield grok.functions.eval('Chem:findSimilar');
    let call = yield foo.prepare({
      'molStringsColumn': column,
      'molString': molecule,
      'limit': settings.limit,
      'cutoff': settings.cutoff
    });
    yield call.call();
    return (molecule.length != 0) ? call.getParamValue('result') : null;
  });
}
/**
 * Returns molecules similar to the reference one.
 * @async
 * @param {Column} column - Molecule column to search in.
 * @param {string} molecule - Reference molecule in SMILES format.
 * @param {SimilarityMetric} metric - Metric to use.
 * @param {number} limit - Maximum number of results to return.
 * @param {number} minScore - Minimum similarity score for a molecule to be included.
 * @returns {Promise<DataFrame>}
 * */
export function findSimilarServer(column, molecule, metric = SIMILARITY_METRIC.TANIMOTO, limit = 10, minScore = 0.7) {
  return new Promise((resolve, reject) => api.grok_Chem_SimilaritySearch(column.d, molecule, metric, limit, minScore, (t) => resolve(new DataFrame(t))));
}
/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 * @async
 * @param {Column} column - Column with molecules to search.
 * @param {string} pattern - Pattern, either SMARTS or SMILES.
 * @param {boolean} isSmarts - Whether the pattern is SMARTS.
 * @returns {Promise<BitSet>}
 * */
export function searchSubstructureServer(column, pattern, isSmarts = true) {
  return new Promise((resolve, reject) => api.grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs) => resolve(new BitSet(bs))));
}
/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/substructure-search}
 * @async
 * @param {Column} column - Column with molecules to search
 * @param {string} molecule - Substructure being sought, either one of which RDKit supports:
 *   smiles, cxsmiles, molblock, v3Kmolblock, and inchi
 * @param settings
 * @returns {Promise<BitSet>}
 * */
export function searchSubstructure(column, molecule = '', settings) {
  return __awaiter(this, void 0, void 0, function* () {
    return substructureSearch(column, molecule, settings);
  });
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
  return new Promise((resolve, reject) => api.grok_Chem_RGroup(table.d, column, core, () => resolve(table), (e) => reject(e)));
}
/**
 * Finds Most Common Substructure in the specified column.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
 * @async
 * @param {Column} column - Column with SMILES to analyze.
 * @returns {Promise<string>}
 * */
export function mcs(column) {
  return new Promise((resolve, reject) => api.grok_Chem_MCS(column.d, (mcs) => resolve(mcs), (e) => reject(e)));
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
  return new Promise((resolve, reject) => api.grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table), (e) => reject(e)));
}
/**
 * Returns available descriptors tree.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 *
 * @async
 * @returns {Promise<Object>}
 * */
export function descriptorsTree() {
  return new Promise((resolve, reject) => api.grok_Chem_DescriptorsTree((tree) => resolve(JSON.parse(tree)), (e) => reject(e)));
}
/**
 * Renders a molecule to SVG
 * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mol-rendering}
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
export function sketcher(onChangedCallback, smiles = '') {
  return api.grok_Chem_Sketcher(onChangedCallback, smiles);
}
