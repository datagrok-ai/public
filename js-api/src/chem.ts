/**
 * Cheminformatics support
 * @module chem
 * */

import {BitSet, Column, DataFrame} from './dataframe';
import {SIMILARITY_METRIC, SimilarityMetric, TYPE} from './const';
import {Observable, Subject, Subscription} from "rxjs";
import {Menu, Widget} from "./widgets";
import {Func} from "./entities";
import * as ui from "../ui";
import {SemanticValue} from "./grid";
import $ from "cash-dom";
import {element} from "../ui";

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


  /**
   * Molecule sketcher that supports multiple dynamically initialized implementations.
   * */
  export class Sketcher extends Widget {

    molInput: HTMLInputElement = ui.element('input');
    host: HTMLDivElement = ui.box(null, 'grok-sketcher');
    /* molInputSubscription: Subscription | null = null; */
    changedSub: Subscription | null = null;
    sketcher: SketcherBase | null = null;
    onChanged: Subject<any> = new Subject<any>();
    listeners: Function[] = [];

    _smiles: string = '';
    _molFile: string = '';
    _smarts: string = '';

    getSmiles(): string {
      return this.sketcher ? this.sketcher.smiles : this._smiles;
    }

    setSmiles(x: string) {
      this._smiles = x;
      if (this.sketcher != null)
        this.sketcher!.smiles = x;
    }

    getMolFile(): string {
      return this.sketcher ? this.sketcher.molFile : this._molFile;
    }

    setMolFile(x: string) {
      this._molFile = x;
      if (this.sketcher != null)
        this.sketcher!.molFile = x;
    }

    get supportedExportFormats(): string[] {
      return this.sketcher ? this.sketcher.supportedExportFormats : [];
    }

    async getSmarts(): Promise<string> {
      return this.sketcher!.getSmarts();
    }

    setChangeListenerCallback(callback: () => void) {
      this.changedSub?.unsubscribe();
      this.listeners.push(callback);
      if (this.sketcher)
        this.changedSub = this.sketcher.onChanged.subscribe((_: any) => callback());
    }

    constructor() {
      super(ui.div());

      let funcs = Func.find({tags: ['moleculeSketcher']});
      if (funcs.length == 0)
        throw 'Sketcher functions not found. Please install OpenChemLib, or MarvinJS package.';

      
      $(this.molInput).attr('placeholder', 'SMILES, Inchi, Inchi keys, ChEMBL id, etc');
      const smilesInputHandler = (e: any) => {
        const newSmilesValue: string = (e?.target as HTMLTextAreaElement).value;
        if (this.getSmiles() !== newSmilesValue) {
          this.setSmiles(newSmilesValue);
        }
        const currentSmiles = this.getSmiles();
        if (currentSmiles !== newSmilesValue) {
          (e?.target as HTMLTextAreaElement).value = currentSmiles ?? '';
        }
      };
      this.molInput.addEventListener('focusout', (e) => {
        smilesInputHandler(e);
      });
      this.molInput.addEventListener('keydown', (e) => {
        if (e.keyCode == 13)
          smilesInputHandler(e);
      });

      let optionsIcon = ui.iconFA('bars', () => {
        Menu
          .popup()
          .item('Add to favorites', () => console.log(this.sketcher!.molFile))
          .separator()
          .items(funcs.map((f) => f.name), (name: string) => this.setSketcher(name))
          .show();
      });
      $(optionsIcon).addClass('d4-input-options');

      let molInputDiv = ui.div([this.molInput, optionsIcon], 'grok-sketcher-input');

      //let sketcherChoice = ui.choiceInput('Sketcher', funcs[0].name, funcs.map((f) => f.name), (name: string) => this.setSketcher(name))
      this.root.appendChild(ui.div([
        molInputDiv,
        //sketcherChoice.root,
        this.host]));

      this.setSketcher(funcs[0].name);
    }

    detach() {
      this.changedSub?.unsubscribe();
      this.sketcher?.detach();
      super.detach();
    }

    async setSketcher(name: string) {
      ui.empty(this.host);
      this.changedSub?.unsubscribe();

      let molFile = this.sketcher?.molFile ?? this._molFile;
      let smiles = this.sketcher?.smiles ?? this._smiles;

      let f = Func.find({name: name})[0];
      this.sketcher = await f.apply();
      /*
      // A backward connection from the sketched molecule to the text
      this.molInputSubscription?.unsubscribe();
      this.molInputSubscription = this.sketcher!.onChanged.subscribe((_) => {
        this.molInput.value = this.getSmiles() ?? ''; });
      */
      this.host.appendChild(this.sketcher!.root);
      await ui.tools.waitForElementInDom(this.root);
      await this.sketcher!.init();
      this.changedSub = this.sketcher!.onChanged.subscribe((_: any) => {
        this.onChanged.next(null);
        for (let callback of this.listeners)
          callback();
        grok.shell.o = SemanticValue.fromValueType(this.sketcher!.smiles, 'Molecule');
      });

      if (molFile != null)
        this.sketcher!.molFile = molFile;
      else if (smiles != null)
        this.sketcher!.smiles = smiles;
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

    const result = await grok.functions.call('Chem:similarityScoring', {
      'molStringsColumn': column,
      'molString': molecule,
      'sorted': settings.sorted
    });
    if (molecule.length != 0) {
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
    return new Promise((resolve, reject) => api.grok_Chem_DiversitySearch(column.dart, metric, limit, (mols: any) => resolve(mols), (e: any) => reject(e)));
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

    return (await grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn': column,
      'molString': pattern,
      'substructLibrary':
        !(settings?.hasOwnProperty('substructLibrary') && !settings.substructLibrary),
      'molStringSmarts': ''
    })).get(0);
    
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

    const result = await grok.functions.call('Chem:getSimilarities', {
      'molStringsColumn': column,
      'molString': molecule
    });
    // TODO: figure out what's the state in returning columns from package functions
    return (molecule.length != 0) ? result.columns.byIndex(0) : null;

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

    const result = await grok.functions.call('Chem:findSimilar', {
      'molStringsColumn': column,
      'molString': molecule,
      'limit': settings.limit,
      'cutoff': settings.cutoff
    });
    return (molecule.length != 0) ? result : null;

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
    return new Promise((resolve, reject) => api.grok_Chem_SimilaritySearch(column.dart, molecule, metric,
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
    return new Promise((resolve, reject) => api.grok_Chem_SubstructureSearch(column.dart, pattern, isSmarts, (bs: any) => resolve(new BitSet(bs))));
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
    return new Promise((resolve, reject) => api.grok_Chem_RGroup(table.dart, column, core, () => resolve(table), (e: any) => reject(e)));
  }

  /**
   * Finds Most Common Substructure in the specified column.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
   * @async
   * @param {Column} column - Column with SMILES to analyze.
   * @returns {Promise<string>}
   * */
  export function mcs(column: Column): Promise<string> {
    return new Promise((resolve, reject) => api.grok_Chem_MCS(column.dart, (mcs: any) => resolve(mcs), (e: any) => reject(e)));
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
    return new Promise((resolve, reject) => api.grok_Chem_Descriptors(table.dart, column, descriptors, () => resolve(table), (e: any) => reject(e)));
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
   * @param {string} smiles - accepts smiles/molfile format
   * @param {number} width
   * @param {number} height
   * @param {object} options - OCL.IMoleculeToSVGOptions
   * @returns {HTMLDivElement}
   * */
  export function svgMol(
      smiles: string, width: number = 300, height: number = 200,
      options?: {[key: string]: boolean | number | string}
    ): HTMLDivElement {
    let root = document.createElement('div');
    import('openchemlib/full.js').then((OCL) => {
      let m = smiles.endsWith("M END") ? OCL.Molecule.fromMolfile(smiles): OCL.Molecule.fromSmiles(smiles);
      root.innerHTML = m.toSVG(width, height, undefined, options);
    });
    return root;
  }
  
    /**
   * Renders a molecule to canvas (using RdKit)
   * TODO: should NOT be async
   * See example: {@link }
   * @param {string} smiles - accepts smiles/molfile format
   * @param {number} x
   * @param {number} y
   * @param {number} w
   * @param {number} h
   * @param {Object} canvas
   * @param {string} molString
   * @param {string} scaffoldMolString
   * */
  export async function canvasMol(
    x: number, y: number, w: number, h: number,
    canvas: Object, molString: string, scaffoldMolString: string | null = null): Promise<void> {
      await grok.functions.call('Chem:canvasMol', {
        'x': x, 'y': y, 'w': w, 'h': h,
        'canvas': canvas, 'molString': molString,
        'scaffoldMolString': scaffoldMolString ?? '' 
      });
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

    /*
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
    */
  }
}