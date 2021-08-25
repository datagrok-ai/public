/**
 * Cheminformatics support
 * @module chem
 * */

import {BitSet, Column, DataFrame} from './dataframe';
import {SIMILARITY_METRIC, SimilarityMetric, TYPE} from './const';
import {Observable, Subject, Subscription} from "rxjs";
import {Widget} from "./widgets";
import {Func} from "./entities";
import * as ui from "../ui";
import {SemanticValue} from "./grid";

let api = <any>window;
declare let grok: any;

/** Cheminformatics-related routines */
export namespace chem {

  /** A common interface that all sketchers should implement */
  export abstract class SketcherBase extends Widget {
    readonly SMILES: string = 'smiles';
    readonly SMARTS = 'smarts';
    readonly MOL = 'mol';

    onChanged: Subject<any> = new Subject<any>();
    _smiles: string;
    _smarts: string;
    _molFile: string;

    constructor() {
      super(ui.box());

      this._smiles = this.addProperty('smiles', TYPE.STRING);
      this._smarts = this.addProperty('smarts', TYPE.STRING);
      this._molFile = this.addProperty('molFile', TYPE.STRING);
    }

    /** SMILES representation of the molecule */
    get smiles(): string {
      return this._smiles;
    }

    set smiles(s: string) {
      this._smiles = s;
    }

    get supportedExportFormats(): string[] {
      return [];
    }

    async exportStructure(format: string): Promise<string> {
      return 'dummy';
    }

    /** SMARTS query */
    async getSmarts(): Promise<string> {
      return this.exportStructure('smarts');
    }

    setSmarts(s: string) {
      this._smarts = s;
    }

    /** MolFile representation of the molecule */
    get molFile(): string {
      return this._molFile;
    }

    set molFile(s: string) {
      this._molFile = s;
    }

    /** Override to provide custom initialization. At this point, the root is already in the DOM. */
    async init() {
    }
  }


  export class Sketcher extends Widget {

    host: HTMLDivElement = ui.box(null, {style: {width: '500px', height: '500px'}});
    changedSub: Subscription | null = null;
    sketcher: SketcherBase | null = null;
    listeners: Function[] = [];

    getSmiles(): string {
      return this.sketcher!.smiles;
    }

    setSmiles(x: string) {
      this.sketcher!.smiles = x;
    }

    getMolFile(): string {
      return this.sketcher!.molFile;
    }

    setMolFile(x: string) {
      this.sketcher!.molFile = x;
    }

    async getSmarts(): Promise<string> {
      return this.sketcher!.getSmarts();
    }

    setChangeListenerCallback(callback: () => void) {
      this.listeners.push(callback);
      this.sketcher!.onChanged.subscribe((_) => callback());
    }

    constructor() {
      super(ui.div());

      let funcs = Func.find({tags: ['moleculeSketcher']});
      let input = ui.choiceInput('Sketcher', funcs[0].name, funcs.map((f) => f.name), (name: string) => this.setSketcher(name))
      this.root.appendChild(ui.div([input.root, this.host]));

      this.setSketcher(funcs[0].name);
    }

    async setSketcher(name: string) {
      ui.empty(this.host);
      this.changedSub?.unsubscribe();
      let molFile = this.sketcher?.molFile;
      let f = Func.find({name: name})[0];
      this.sketcher = await f.apply();
      this.host.appendChild(this.sketcher!.root);
      await ui.tools.waitForElementInDom(this.root);
      await this.sketcher!.init();
      this.changedSub = this.sketcher!.onChanged.subscribe((_) => {
        this.listeners.forEach((callback) => callback());
        grok.shell.o = SemanticValue.fromValueType(this.sketcher!.smiles, 'Molecule');
      });
      if (molFile != null)
        this.sketcher!.molFile = molFile;
    }
  }


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
  export async function similarityScoring(column: Column, molecule: string = '', settings: { sorted?: boolean } = {sorted: false}) {

    let foo = await grok.functions.eval('Chem:similarityScoring');
    let call = await foo.prepare({
      'molStringsColumn': column,
      'molString': molecule,
      'sorted': settings.sorted
    });
    await call.call();
    if (molecule.length != 0) {
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
  export function diversitySearch(column: Column, metric: SimilarityMetric = SIMILARITY_METRIC.TANIMOTO, limit: number = 10): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_Chem_DiversitySearch(column.d, metric, limit, (mols: any) => resolve(mols), (e: any) => reject(e)));
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
  export async function searchSubstructure(column: Column, pattern: string = '', settings: {
    substructLibrary?: boolean;
  } = {}): Promise<BitSet> {

    let foo = await grok.functions.eval('Chem:searchSubstructure');
    let call = await foo.prepare({
      'molStringsColumn': column,
      'molString': pattern,
      'substructLibrary':
        !(settings?.hasOwnProperty('substructLibrary') && !settings.substructLibrary)
    });
    await call.call();
    // unpacking our BitSet object from a synthetic column
    return call.getParamValue('result').get(0);

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
  export async function getSimilarities(column: Column, molecule: string = '', settings: object = {}): Promise<Column | null> {

    let foo = await grok.functions.eval('Chem:getSimilarities');
    let call = await foo.prepare({
      'molStringsColumn': column,
      'molString': molecule
    });
    await call.call();
    // TODO: figure out what's the state in returning columns from package functions
    return (molecule.length != 0) ? call.getParamValue('result').columns.byIndex(0) : null;

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
  export async function findSimilar(column: Column, molecule: string = '', settings = {limit: Number.MAX_VALUE, cutoff: 0.0}): Promise<DataFrame | null> {

    let foo = await grok.functions.eval('Chem:findSimilar');
    let call = await foo.prepare({
      'molStringsColumn': column,
      'molString': molecule,
      'limit': settings.limit,
      'cutoff': settings.cutoff
    });
    await call.call();
    return (molecule.length != 0) ? call.getParamValue('result') : null;

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
  export function findSimilarServer(column: Column, molecule: string, metric: SimilarityMetric = SIMILARITY_METRIC.TANIMOTO, limit: number = 10, minScore: number = 0.7): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_Chem_SimilaritySearch(column.d, molecule, metric,
      limit, minScore, (t: any) => resolve(new DataFrame(t))));
  }

  /**
   * Searches for a molecular pattern in a given column, returning a bitset with hits.
   * @async
   * @param {Column} column - Column with molecules to search.
   * @param {string} pattern - Pattern, either SMARTS or SMILES.
   * @param {boolean} isSmarts - Whether the pattern is SMARTS.
   * @returns {Promise<BitSet>}
   * */
  export function searchSubstructureServer(column: Column, pattern: string, isSmarts: boolean = true): Promise<BitSet> {
    return new Promise((resolve, reject) => api.grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs: any) => resolve(new BitSet(bs))));
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
  export async function substructureSearch(column: Column, molecule: string = '', settings: {
    substructLibrary: boolean;
  }): Promise<BitSet> {
    return searchSubstructure(column, molecule, settings);
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
  export function rGroup(table: DataFrame, column: string, core: string): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_Chem_RGroup(table.d, column, core, () => resolve(table), (e: any) => reject(e)));
  }

  /**
   * Finds Most Common Substructure in the specified column.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
   * @async
   * @param {Column} column - Column with SMILES to analyze.
   * @returns {Promise<string>}
   * */
  export function mcs(column: Column): Promise<string> {
    return new Promise((resolve, reject) => api.grok_Chem_MCS(column.d, (mcs: any) => resolve(mcs), (e: any) => reject(e)));
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
  export function descriptors(table: DataFrame, column: string, descriptors: string[]): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table), (e: any) => reject(e)));
  }

  /**
   * Returns available descriptors tree.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
   *
   * @async
   * @returns {Promise<Object>}
   * */
  export function descriptorsTree(): Promise<object> {
    return new Promise((resolve, reject) => api.grok_Chem_DescriptorsTree((tree: any) => resolve(JSON.parse(tree)), (e: any) => reject(e)));
  }

  /**
   * Renders a molecule to SVG
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mol-rendering}
   * @param {string} smiles
   * @param {number} width
   * @param {number} height
   * @returns {HTMLDivElement}
   * */
  export function svgMol(smiles: string, width: number = 300, height: number = 200): HTMLDivElement {
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
  export function sketcher(onChangedCallback: Function, smiles: string = ''): HTMLElement {
    return api.grok_Chem_Sketcher(onChangedCallback, smiles);
  }

  export async function createSketcher(): Promise<SketcherBase> {
    let func = Func.find({name: 'createMarvinSketcher'})[0];
    return func.apply();
  }

  export function convert(s: string, sourceFormat: string, targetFormat: string) {
    if (sourceFormat == 'mol' && targetFormat == 'smiles') {
      // @ts-ignore
      let mol = new OCL.Molecule.fromMolfile(s);
      return mol.toSmiles();
    }
  }

  export async function showSketcherDialog() {

    ui.dialog()
      .add(new Sketcher().root)
      .show();
    return;

    let funcs = Func.find({tags: ['moleculeSketcher']});
    let host = ui.box(null, {style: {width: '500px', height: '500px'}});
    let changedSub: Subscription | null = null;
    let sketcher: SketcherBase | null = null;

    async function setSketcher(name: string) {
      ui.empty(host);
      changedSub?.unsubscribe();
      let molFile = sketcher?.molFile;
      let f = Func.find({name: name})[0];
      sketcher = await f.apply();
      host.appendChild(sketcher!.root);
      await sketcher!.init();
      changedSub = sketcher!.onChanged.subscribe((_) => grok.shell.o = SemanticValue.fromValueType(sketcher!.smiles, 'Molecule'));
      if (molFile != null)
        sketcher!.molFile = molFile;
    }

    ui.dialog()
      .add(ui.choiceInput('Sketcher', funcs[0].name, funcs.map((f) => f.name), setSketcher))
      .add(host)
      .show();

    await setSketcher(funcs[0].name);
  }
}