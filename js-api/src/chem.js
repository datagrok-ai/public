/**
 * Cheminformatics support
 * @module chem
 * */
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
  Object.defineProperty(o, 'default', { enumerable: true, value: v });
}) : function(o, v) {
  o['default'] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) if (k !== 'default' && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
  __setModuleDefault(result, mod);
  return result;
};
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
  function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
  return new (P || (P = Promise))(function (resolve, reject) {
    function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
    function rejected(value) { try { step(generator['throw'](value)); } catch (e) { reject(e); } }
    function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
    step((generator = generator.apply(thisArg, _arguments || [])).next());
  });
};
define(['require', 'exports', './dataframe', './const'], function (require, exports, dataframe_1, const_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.sketcher = exports.svgMol = exports.descriptorsTree = exports.descriptors = exports.mcs = exports.rGroup = exports.searchSubstructure = exports.searchSubstructureServer = exports.findSimilarServer = exports.findSimilar = exports.getSimilarities = exports.substructureSearch = exports.diversitySearch = exports.similarityScoring = void 0;
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
  function similarityScoring(column, molecule = '', settings = { sorted: false }) {
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
  exports.similarityScoring = similarityScoring;
  /**
     * Returns the specified number of most diverse molecules in the column.
     * See example: {@link https://datagrok.ai/help/domains/chem/diversity-search}
     * @async
     * @param {Column} column - Column with molecules in which to search.
     * @param {SimilarityMetric} metric - Metric to use.
     * @param {number} limit - Number of molecules to return.
     * @returns {Promise<DataFrame>}
     * */
  function diversitySearch(column, metric = const_1.SIMILARITY_METRIC.TANIMOTO, limit = 10) {
    return new Promise((resolve, reject) => api.grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols), (e) => reject(e)));
  }
  exports.diversitySearch = diversitySearch;
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
  function substructureSearch(column, pattern = '', settings) {
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
  exports.substructureSearch = substructureSearch;
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
  function getSimilarities(column, molecule = '', settings = {}) {
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
  exports.getSimilarities = getSimilarities;
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
  function findSimilar(column, molecule = '', settings = { limit: Number.MAX_VALUE, cutoff: 0.0 }) {
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
  exports.findSimilar = findSimilar;
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
  function findSimilarServer(column, molecule, metric = const_1.SIMILARITY_METRIC.TANIMOTO, limit = 10, minScore = 0.7) {
    return new Promise((resolve, reject) => api.grok_Chem_SimilaritySearch(column.d, molecule, metric, limit, minScore, (t) => resolve(new dataframe_1.DataFrame(t))));
  }
  exports.findSimilarServer = findSimilarServer;
  /**
     * Searches for a molecular pattern in a given column, returning a bitset with hits.
     * @async
     * @param {Column} column - Column with molecules to search.
     * @param {string} pattern - Pattern, either SMARTS or SMILES.
     * @param {boolean} isSmarts - Whether the pattern is SMARTS.
     * @returns {Promise<BitSet>}
     * */
  function searchSubstructureServer(column, pattern, isSmarts = true) {
    return new Promise((resolve, reject) => api.grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs) => resolve(new dataframe_1.BitSet(bs))));
  }
  exports.searchSubstructureServer = searchSubstructureServer;
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
  function searchSubstructure(column, molecule = '', settings) {
    return __awaiter(this, void 0, void 0, function* () {
      return substructureSearch(column, molecule, settings);
    });
  }
  exports.searchSubstructure = searchSubstructure;
  /**
     * Performs R-group analysis.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
     * @async
     * @param {DataFrame} table - Table.
     * @param {string} column - Column name with SMILES to analyze.
     * @param {string} core - Core in the SMILES format.
     * @returns {Promise<DataFrame>}
     * */
  function rGroup(table, column, core) {
    return new Promise((resolve, reject) => api.grok_Chem_RGroup(table.d, column, core, () => resolve(table), (e) => reject(e)));
  }
  exports.rGroup = rGroup;
  /**
     * Finds Most Common Substructure in the specified column.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
     * @async
     * @param {Column} column - Column with SMILES to analyze.
     * @returns {Promise<string>}
     * */
  function mcs(column) {
    return new Promise((resolve, reject) => api.grok_Chem_MCS(column.d, (mcs) => resolve(mcs), (e) => reject(e)));
  }
  exports.mcs = mcs;
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
  function descriptors(table, column, descriptors) {
    return new Promise((resolve, reject) => api.grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table), (e) => reject(e)));
  }
  exports.descriptors = descriptors;
  /**
     * Returns available descriptors tree.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
     *
     * @async
     * @returns {Promise<Object>}
     * */
  function descriptorsTree() {
    return new Promise((resolve, reject) => api.grok_Chem_DescriptorsTree((tree) => resolve(JSON.parse(tree)), (e) => reject(e)));
  }
  exports.descriptorsTree = descriptorsTree;
  /**
     * Renders a molecule to SVG
     * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mol-rendering}
     * @param {string} smiles
     * @param {number} width
     * @param {number} height
     * @returns {HTMLDivElement}
     * */
  function svgMol(smiles, width = 300, height = 200) {
    let root = document.createElement('div');
    new Promise((resolve_1, reject_1) => { require(['openchemlib/full.js'], resolve_1, reject_1); }).then(__importStar).then((OCL) => {
      let m = OCL.Molecule.fromSmiles(smiles);
      root.innerHTML = m.toSVG(width, height);
    });
    return root;
  }
  exports.svgMol = svgMol;
  /**
     * Sketches Molecule sketcher.
     * @param {function} onChangedCallback - a function that accepts (smiles, molfile)
     * @param {string} smiles Initial molecule
     * @returns {HTMLElement}
     * */
  function sketcher(onChangedCallback, smiles = '') {
    return api.grok_Chem_Sketcher(onChangedCallback, smiles);
  }
  exports.sketcher = sketcher;
});
