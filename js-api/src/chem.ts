/**
 * Cheminformatics support
 * @module chem
 * */

import {BitSet, Column, DataFrame} from './dataframe';
import {FUNC_TYPES, SEMTYPE, SIMILARITY_METRIC, SimilarityMetric, TYPE, UNITS} from './const';
import {Subject, Subscription} from "rxjs";
import {Menu, Widget} from "./widgets";
import {Func} from "./entities";
import * as ui from "../ui";
import {SemanticValue} from "./grid";
import $ from "cash-dom";
import {Utils} from "./utils";

let api = <any>window;
declare let grok: any;

const STORAGE_NAME = 'sketcher';
const KEY = 'selected';
const DEFAULT_SKETCHER = 'openChemLibSketcher';

let extractors: Func[];  // id => molecule

export function isMolBlock(s: string | null) {
  return s != null && s.includes('M  END');
}

/** Cheminformatics-related routines */
export namespace chem {

  export enum SKETCHER_MODE {
    INPLACE = 'Inplace',
    EXTERNAL = 'External'
  }

  /** A common interface that all sketchers should implement */
  export abstract class SketcherBase extends Widget {
    readonly SMILES: string = 'smiles';
    readonly SMARTS = 'smarts';
    readonly MOL = 'mol';

    onChanged: Subject<any> = new Subject<any>();
    _smiles: string;
    _smarts: string;
    _molFile: string;
    _mode: string;

    constructor() {
      super(ui.box());

      this._smiles = this.addProperty('smiles', TYPE.STRING);
      this._smarts = this.addProperty('smarts', TYPE.STRING);
      this._molFile = this.addProperty('molFile', TYPE.STRING);
      this._mode = this.addProperty('mode', TYPE.STRING, SKETCHER_MODE.EXTERNAL, { choices: [ SKETCHER_MODE.EXTERNAL, SKETCHER_MODE.INPLACE ] });
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

    async exportStructure(_format: string): Promise<string> {
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

    /** Sketcher mode */
    get mode(): string {
      return this._mode;
    }

    set mode(s: string) {
      this._mode = s;
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
    changedSub: Subscription | null = null;
    sketcher: SketcherBase | null = null;
    onChanged: Subject<any> = new Subject<any>();

    /** Whether the currently drawn molecule becomes the current object as you sketch it */
    syncCurrentObject: boolean = true;

    listeners: Function[] = [];

    _smiles: string = '';
    _molFile: string = '';
    _smarts: string = '';
    _mode = SKETCHER_MODE.INPLACE;

    extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});

    getMode(): string {
      return this.sketcher ? this.sketcher.mode : this._mode;
    }

    setMode(x: SKETCHER_MODE) {
      this._mode = x;
      if (this.sketcher != null)
        this.sketcher!.mode = x;
    }

    getSmiles(): string {
      return this.sketcher ? this.sketcher.smiles : this._smiles;
    }

    setSmiles(x: string) {
      this._smiles = x;
      this._molFile = convert(x, 'smiles', 'mol');
      if (this.sketcher != null)
        this.sketcher!.smiles = x;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    getMolFile(): string {
      return this.sketcher ? this.sketcher.molFile : this._molFile;
    }

    setMolFile(x: string) {
      this._molFile = x;
      this._smiles = convert(x, 'mol', 'smiles');
      if (this.sketcher != null)
        this.sketcher!.molFile = x;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    get isEmpty(): boolean {
      return (this.getSmiles() == null || this.getSmiles() == '') &&
        (this.getMolFile() == null || this.getMolFile() == '' || this.getMolFile().split("\n")[3].trimLeft()[0] === '0');
    }

    /** Sets the molecule, supports either SMILES or MOLBLOCK formats */
    setMolecule(molString: string) {
      if (isMolBlock(molString))
        this.setMolFile(molString)
      else
        this.setSmiles(molString);
    }

    get supportedExportFormats(): string[] {
      return this.sketcher ? this.sketcher.supportedExportFormats : [];
    }

    async getSmarts(): Promise<string> {
      return this.sketcher ? await this.sketcher.getSmarts() : '';
    }

    setChangeListenerCallback(callback: () => void) {
      this.changedSub?.unsubscribe();
      this.listeners.push(callback);
      if (this.sketcher)
        this.changedSub = this.sketcher.onChanged.subscribe((_: any) => callback());
    }

    /** Sets SMILES, MOLBLOCK, or any other molecule representation */
    setValue(x: string) {
      const extractor = extractors
        .find((f) => new RegExp(f.options['inputRegexp']).test(x));

      if (extractor != null)
        extractor
          .apply([ new RegExp(extractor.options['inputRegexp']).exec(x)![1] ])
          .then((mol) => this.setMolecule(mol));
      else
        this.setMolecule(x);
    }

    constructor(mode?: SKETCHER_MODE) {
      super(ui.div());
      if (mode)
        this.setMode(mode);
      this.root.append(ui.div([ ui.divText('') ]));
      setTimeout(() => this.createSketcher(), 100);
    }

    /** In case sketcher is opened in filter panel use EXTERNAL mode*/
    setExternalModeForSubstrFilter() {
      if (this.root.closest('.d4-filter'))
        this.setMode(SKETCHER_MODE.EXTERNAL);
    }

    createSketcher() {
      this.setExternalModeForSubstrFilter();
      this.root.innerHTML = '';
      if (this._mode === SKETCHER_MODE.INPLACE)
        this.root.appendChild(this.createInplaceModeSketcher());
      else
        this.root.appendChild(this.createExternalModeSketcher());
    }

    _updateExtSketcherInnerHTML(content: HTMLElement) {
      this.extSketcherDiv.innerHTML = '';
      this.extSketcherDiv.append(content);
    }

    updateExtSketcherContent(extSketcherDiv: HTMLElement) {
      if (!this.isEmpty && extSketcherDiv.parentElement) {
        const width = extSketcherDiv.parentElement!.clientWidth;
        const height = width / 2;
        ui.empty(this.extSketcherDiv);
        let canvas = ui.canvas(width, height);
        canvas.style.height = '100%';
        canvas.style.width = '100%';
        canvasMol(0, 0, width, height, canvas, this.getMolFile())
          .then((_) => {
            ui.empty(this.extSketcherDiv);
            this.extSketcherDiv.append(canvas);
          });
      }

      let sketchLink = ui.button('Sketch', () => this.updateExtSketcherContent(extSketcherDiv));
      sketchLink.style.paddingLeft = '0px';
      sketchLink.style.marginLeft = '0px';
      this._updateExtSketcherInnerHTML(sketchLink);
    };

    createExternalModeSketcher(): HTMLElement {
      this.extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});
      ui.tooltip.bind(this.extSketcherDiv, 'Click to edit filter');

      this.extSketcherDiv.addEventListener('mousedown', () => {
        let savedMolFile = this.getMolFile();
        ui.dialog()
          .add(this.createInplaceModeSketcher())
          .onCancel(() => this.setMolFile(savedMolFile))
          .onOK(() => {
            this.updateExtSketcherContent(this.extSketcherDiv);
            Sketcher.addRecent(savedMolFile);
          })
          .show();
      });

      ui.onSizeChanged(this.extSketcherDiv).subscribe((_) => {
        this.updateExtSketcherContent(this.extSketcherDiv);
      });

      this.updateExtSketcherContent(this.extSketcherDiv);
      return this.extSketcherDiv;
    }

    createInplaceModeSketcher(): HTMLElement {
      const molInputDiv = ui.div();
      grok.dapi.userDataStorage.getValue(STORAGE_NAME, KEY, true).then((sname: string) => {
        let funcs = Func.find({ tags: [ 'moleculeSketcher' ] });
        let fr = funcs.find(e => e.friendlyName == sname || e.name == sname)
          ?? funcs.find(e => e.name == DEFAULT_SKETCHER);

        $(this.molInput).attr('placeholder', 'SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');

      if (extractors == null) {
        const extractorSearchOptions = {
          meta: {
            role: FUNC_TYPES.CONVERTER,
            inputRegexp: null,
          },
          returnSemType: SEMTYPE.MOLECULE
        };

        const load: Promise<any> = api.grok_Func_LoadQueriesScripts();
        load
          .then((_) => { extractors = Func.find(extractorSearchOptions); })
          .catch((_) => extractors = []);
      }

      const applyInput = (e: any) => {
        const newSmilesValue: string = (e?.target as HTMLTextAreaElement).value;

        if (this.getSmiles() !== newSmilesValue)
          this.setValue(newSmilesValue);

        const currentSmiles = this.getSmiles();

        if (currentSmiles !== newSmilesValue)
          (e?.target as HTMLTextAreaElement).value = currentSmiles ?? '';
      };

      this.molInput.addEventListener('keydown', (e) => {
        if (e.key == 'Enter') {
          applyInput(e);
          e.stopImmediatePropagation();
        }
      });

      this.molInput.addEventListener('paste', (e) => {
        const text = e.clipboardData?.getData('text/plain');
        if (text != null && isMolBlock(text)) {
          e.preventDefault();
          this.setValue(text);
        }
      });

      let optionsIcon = ui.iconFA('bars', () => {
        Menu.popup()
          .item('Copy as SMILES', () => navigator.clipboard.writeText(this.sketcher!.smiles))
          .item('Copy as MOLBLOCK', () => navigator.clipboard.writeText(this.sketcher!.molFile))
          .group('Recent')
            .items(Sketcher.getRecent().map((m) => ui.tools.click(svgMol(m, 100, 70), () => this.setMolecule(m))), () => { })
          .endGroup()
          .group('Favorites')
            .item('Add to Favorites', () => Sketcher.addFavorite(this.sketcher!.molFile))
            .separator()
            .items(Sketcher.getFavorites().map((m) => ui.tools.click(svgMol(m, 100, 70), () => this.setMolecule(m))), () => { })
          .endGroup()
          .separator()
          .items(funcs.map((f) => f.friendlyName), (name: string) => this.setSketcher(name),
            { isChecked: (item) => item === sname, toString: item => item })
          .show();
      });
      $(optionsIcon).addClass('d4-input-options');
        molInputDiv.append(ui.div([ this.molInput, optionsIcon ], 'grok-sketcher-input'));
        this.setSketcher(fr!.friendlyName);
      });

      return ui.div([
        molInputDiv,
        this.host]);
    }

    static readonly FAVORITES_KEY = 'chem-molecule-favorites';
    static readonly RECENT_KEY = 'chem-molecule-recent';

    static getFavorites(): string[] {
      return JSON.parse(localStorage.getItem(Sketcher.FAVORITES_KEY) ?? '[]');
    }

    static addFavorite(molecule: string) {
      let s = JSON.stringify([...Sketcher.getFavorites().slice(-9), molecule]);
      localStorage.setItem(Sketcher.FAVORITES_KEY, s);
    }

    static getRecent(): string[] {
      return JSON.parse(localStorage.getItem(Sketcher.RECENT_KEY) ?? '[]');
    }

    static addRecent(molecule: string) {
      if (!Sketcher.getRecent().includes(molecule)) {
        let s = JSON.stringify([...Sketcher.getRecent().slice(-9), molecule]);
        localStorage.setItem(Sketcher.RECENT_KEY, s);
      }
    }

    detach() {
      this.changedSub?.unsubscribe();
      this.sketcher?.detach();
      super.detach();
    }

    async setSketcher(name: string) {
      ui.empty(this.host);
      this.changedSub?.unsubscribe();

      if (this.sketcher?.molFile) this._molFile = this.sketcher?.molFile;

      let funcs = Func.find({tags: ['moleculeSketcher']});
      let f = funcs.find(e => e.friendlyName == name || e.name == name);

      grok.dapi.userDataStorage.postValue(STORAGE_NAME, KEY, f!.friendlyName, true);

      this.sketcher = await f!.apply();
      this.host!.style.minWidth = '500px';
      this.host!.style.minHeight = '400px';
      this.host.appendChild(this.sketcher!.root);
      await ui.tools.waitForElementInDom(this.root);
      await this.sketcher!.init();
      this.changedSub = this.sketcher!.onChanged.subscribe((_: any) => {
        this.onChanged.next(null);
        for (let callback of this.listeners)
          callback();
        if (this.syncCurrentObject)
          grok.shell.o = SemanticValue.fromValueType(this.sketcher!.molFile, SEMTYPE.MOLECULE, UNITS.Molecule.MOLBLOCK);
      });

      if (!Utils.isEmpty(this._molFile))
        this.sketcher!.molFile = this._molFile;
      else if (!Utils.isEmpty(this._smiles))
        this.sketcher!.smiles = this._smiles;
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
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/substructure-search-library}
   * @async
   * @param {Column} column - Column with molecules to search
   * @param {string} pattern - Pattern, either one of which RDKit supports
   * @param settings
   * @returns {Promise<BitSet>}
   * */
  export async function searchSubstructure(column: Column, pattern: string = ''): Promise<BitSet> {
    return (await grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn': column,
      'molString': pattern,
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
   * @param {int} settings.limit - Would return top limit molecules based on the score
   * @param {int} settings.cutoff - Would drop molecules which score is lower than cutoff
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
    return new Promise((resolve, _reject) => api.grok_Chem_SimilaritySearch(column.dart, molecule, metric,
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
    return new Promise((resolve, _reject) => api.grok_Chem_SubstructureSearch(column.dart, pattern, isSmarts, (bs: any) => resolve(new BitSet(bs))));
  }

  /**
   * Searches for a molecular pattern in a given column, returning a bitset with hits.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/substructure-search}
   * @async
   * @deprecated
   * @param {Column} column - Column with molecules to search
   * @param {string} molecule - Substructure being sought, either one of which RDKit supports:
   *   smiles, cxsmiles, molblock, v3Kmolblock, and inchi
   * @param settings
   * @returns {Promise<BitSet>}
   * */
  export async function substructureSearch(column: Column, molecule: string = ''): Promise<BitSet> {
    return searchSubstructure(column, molecule);
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
    // @ts-ignore
    import('openchemlib/full.js').then((OCL) => {
      let m = smiles.includes("M  END") ? OCL.Molecule.fromMolfile(smiles): OCL.Molecule.fromSmiles(smiles);
      root.innerHTML = m.toSVG(width, height, undefined, options);
    });
    return root;
  }
  
    /**
   * Renders a molecule to canvas (using RdKit)
   * TODO: should NOT be async
   * See example: {@link }
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
    } else if (sourceFormat == 'smiles' && targetFormat == 'mol'){
      // @ts-ignore
      let mol = new OCL.Molecule.fromSmiles(s);
      return mol.toMolfile();
    }
  }

  export async function showSketcherDialog() {
    ui.dialog()
      .add(new Sketcher().root)
      .show();
  }
}