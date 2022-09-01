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

let api = <any>window;
declare let grok: any;

const STORAGE_NAME = 'sketcher';
const KEY = 'selected';
const DEFAULT_SKETCHER = 'openChemLibSketcher';
const WHITE_MOLBLOCK = `
  Datagrok empty molecule

0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;

let extractors: Func[];  // id => molecule

export function isMolBlock(s: string | null) {
  return s != null && s.includes('M  END');
}

/** Cheminformatics-related routines */
export namespace chem {

  export const SMILES = 'smiles';
  export const MOLV2000 = 'molv2000';
  export const SMARTS = 'smarts';

  export enum SKETCHER_MODE {
    INPLACE = 'Inplace',
    EXTERNAL = 'External'
  }

  /** A common interface that all sketchers should implement */
  export abstract class SketcherBase extends Widget {

    _sketcher: any;
    onChanged: Subject<any> = new Subject<any>();
    host?: Sketcher;

    constructor() {
      super(ui.box());
    }

    /** SMILES representation of the molecule */
    get smiles(): string {
      return this.host!._smiles;
    }

    set smiles(s: string) { }

    /** MolFile representation of the molecule */
    get molFile(): string {
      return this.host!._molfile;
    }

    set molFile(s: string) { }    

    /** SMARTS query */
    async getSmarts(): Promise<string> {
      return this.host!._smarts;
    }

    set smarts(s: string) { }


    get supportedExportFormats(): string[] {
      return [];
    }

    /** Override to provide custom initialization. At this point, the root is already in the DOM. */
    async init(host: Sketcher) {
      this.host = host;
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
    sketcherCreated = new Subject<boolean>();
    sketcherFunctions: Func[] = [];
    selectedSketcher: Func | undefined = undefined;

    /** Whether the currently drawn molecule becomes the current object as you sketch it */
    syncCurrentObject: boolean = true;

    listeners: Function[] = [];
    _mode = SKETCHER_MODE.INPLACE;
    _smiles = '';
    _molfile = '';
    _smarts = '';
    unitsBeforeInit = '';

    extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});
    inplaceSketcherDiv: HTMLDivElement|null = null;

    getSmiles(): string {
      let returnConvertedSmiles = () => { // in case getter is called before sketcher initialized
        if(this._molfile) {
          this._smiles = chem.convert(this._molfile, 'mol', 'smiles');
          return this._smiles;
        } else {
          return this._smarts; //to do - convert from smarts to smiles
        }
      }
      return this.sketcher && this.sketcher._sketcher ? this.sketcher.smiles : !this._smiles ? returnConvertedSmiles() : this._smiles;
    }

    setSmiles(x: string): void {
      this._smiles = x;
      this.sketcher && this.sketcher._sketcher ? this.sketcher!.smiles = x : this.unitsBeforeInit = SMILES;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    getMolFile(): string {
      console.log('#######get molfile js')
      let returnConvertedMolfile = () => { // in case getter is called before sketcher initialized
        if(this._smiles) {
          this._molfile = chem.convert(this._smiles, 'smiles', 'mol');
          return this._molfile;
        } else {
          return this._smarts; //to do - convert from smarts to molfile
        }
      }
      console.log(this.sketcher && this.sketcher._sketcher ? this.sketcher.molFile : !this._molfile ? returnConvertedMolfile() : this._molfile)
      return this.sketcher && this.sketcher._sketcher ? this.sketcher.molFile : !this._molfile ? returnConvertedMolfile() : this._molfile;
    }

    setMolFile(x: string): void {
      this._molfile = x;
      this.sketcher && this.sketcher._sketcher ? this.sketcher!.molFile = x : this.unitsBeforeInit = MOLV2000;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    async getSmarts(): Promise<string | null> {
      let returnConvertedSmarts = async() => { // in case getter is called before sketcher initialized
        if(this._smiles) {
          const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(this._smiles);
          this._smarts = mol.get_smarts();
          mol?.delete();
          return this._smarts;
        } else {
          const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(this._molfile);
          this._smarts = mol.get_smarts();
          mol?.delete();
          return this._smarts;
        }
      }
      return this.sketcher && this.sketcher._sketcher ? await this.sketcher.getSmarts() : !this._smarts ? await returnConvertedSmarts() : this._smarts;
    }

    setSmarts(x: string): void {
      this._smarts = x;
      this.sketcher ? this.sketcher!.smarts = x : this.unitsBeforeInit = SMARTS;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    get supportedExportFormats(): string[] {
      return this.sketcher ? this.sketcher.supportedExportFormats : [];
    }

    isEmpty(): boolean {
      const molFile = this.getMolFile();
      return (molFile == null || molFile == '' || molFile.split("\n")[3].trimStart()[0] === '0');
    }

    /** Sets the molecule, supports either SMILES or MOLBLOCK formats */
    async setMolecule(molString: string, substructure: boolean = false) {
      if(substructure)
        this.setSmarts(molString)
      else if (isMolBlock(molString))
        this.setMolFile(molString)
      else {
        const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molString);
        if (!mol.has_coords())
          mol.set_new_coords();
        mol.normalize_depiction();
        mol.straighten_depiction();
        this._molfile = mol.get_molblock();
        this.setMolFile(this._molfile);
        mol?.delete();
      }
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
        this._mode = mode;
      this.root.append(ui.div([ ui.divText('') ]));
      this.sketcherCreated.subscribe(() => {
        const molecule = this.unitsBeforeInit === SMILES ? this._smiles : this.unitsBeforeInit === MOLV2000 ? this._molfile : this._smarts;
        this.setMolecule(molecule, this.unitsBeforeInit === SMARTS);
      });
      setTimeout(() => this.createSketcher(), 100);
    }

    /** In case sketcher is opened in filter panel use EXTERNAL mode*/
    setExternalModeForSubstrFilter() {
      if (this.root.closest('.d4-filter'))
        this._mode = SKETCHER_MODE.EXTERNAL;
    }

    async createSketcher() {
      const lastSelecttedSketcher = await grok.dapi.userDataStorage.getValue(STORAGE_NAME, KEY, true);
      this.sketcherFunctions = Func.find({ tags: [ 'moleculeSketcher' ] });
      this.selectedSketcher = this.sketcherFunctions.find(e => e.name == lastSelecttedSketcher) ?? this.sketcherFunctions.find(e => e.name == DEFAULT_SKETCHER);
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

    async updateExtSketcherContent(extSketcherDiv: HTMLElement) {
      await ui.tools.waitForElementInDom(extSketcherDiv);
      const width = extSketcherDiv.parentElement!.clientWidth < 100 ? 100 : extSketcherDiv.parentElement!.clientWidth;
      const height = width / 2;
      if (!(this.isEmpty()) && extSketcherDiv.parentElement) {
        ui.empty(this.extSketcherDiv);
        let canvas = ui.canvas(width, height);
        canvas.style.height = '100%';
        canvas.style.width = '100%';
        canvasMol(0, 0, width, height, canvas, this.getMolFile()!)
          .then((_) => {
            ui.empty(this.extSketcherDiv);
            this.extSketcherDiv.append(canvas);
          });
      }

      const sketchLinkStyle = {style: {
        width: `${width}px`, 
        height: `${height/2}px`,
        textAlign: 'center',
        verticalAlign: 'middle',
        lineHeight: `${height/2}px`,
        border: '1px solid #dbdcdf'
      }}
      let sketchLink = ui.divText('Click to edit', sketchLinkStyle);
      sketchLink.onclick = () => this.updateExtSketcherContent(extSketcherDiv);
      sketchLink.style.paddingLeft = '0px';
      sketchLink.style.marginLeft = '0px';
      this._updateExtSketcherInnerHTML(sketchLink);
    };

    createExternalModeSketcher(): HTMLElement {
      this.extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});
      ui.tooltip.bind(this.extSketcherDiv, 'Click to edit');

      this.extSketcherDiv.addEventListener('mousedown', () => {

        let savedMolFile = this.getMolFile();
        savedMolFile = savedMolFile == '' ?  WHITE_MOLBLOCK : savedMolFile;

        let dlg = ui.dialog();
        dlg.add(this.createInplaceModeSketcher(savedMolFile!))
          .onOK(() => {
            this.updateExtSketcherContent(this.extSketcherDiv);
            Sketcher.addRecent(savedMolFile!);
          })
          .onCancel(() => {
            this.setMolFile(savedMolFile!);
          })
          .show();
      });

      ui.onSizeChanged(this.extSketcherDiv).subscribe((_) => {
        this.updateExtSketcherContent(this.extSketcherDiv);
      });

      this.updateExtSketcherContent(this.extSketcherDiv);
      return this.extSketcherDiv;
    }

    createInplaceModeSketcher(molStr?: string): HTMLElement {
      console.log(`create inplace sketcher ${molStr}`)
      const molInputDiv = ui.div();
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
          .item('Copy as SMILES', () => navigator.clipboard.writeText(this.getSmiles()))
          .item('Copy as MOLBLOCK', () => navigator.clipboard.writeText(this.getMolFile()))
          .group('Recent')
          .items(Sketcher.getRecent().map((m) => ui.tools.click(svgMol(m, 100, 70), () => this.setMolecule(m))), () => { })
          .endGroup()
          .group('Favorites')
          .item('Add to Favorites', () => Sketcher.addFavorite(this.getMolFile()))
          .separator()
          .items(Sketcher.getFavorites().map((m) => ui.tools.click(svgMol(m, 100, 70), () => this.setMolecule(m))), () => { })
          .endGroup()
          .separator()
          .items(this.sketcherFunctions.map((f) => f.friendlyName), (friendlyName: string) => {
            this.selectedSketcher = this.sketcherFunctions.filter(f => f.friendlyName === friendlyName)[0];
            this.setSketcher(this.getMolFile());
          },
            { isChecked: (item) => item === this.selectedSketcher?.friendlyName, toString: item => item })
          .show();
      });
      $(optionsIcon).addClass('d4-input-options');
      molInputDiv.append(ui.div([this.molInput, optionsIcon], 'grok-sketcher-input'));
      this.setSketcher(molStr);


      this.inplaceSketcherDiv = ui.div([
        molInputDiv,
        this.host]);

      return this.inplaceSketcherDiv;
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

    async setSketcher(molString?: string) {
      ui.empty(this.host);
      ui.setUpdateIndicator(this.host, true);
      this.changedSub?.unsubscribe();

      grok.dapi.userDataStorage.postValue(STORAGE_NAME, KEY, this.selectedSketcher!.name, true);

      this.sketcher = await this.selectedSketcher!.apply();
      this.host!.style.minWidth = '500px';
      this.host!.style.minHeight = '400px';
      this.host.appendChild(this.sketcher!.root);
      await ui.tools.waitForElementInDom(this.root);
      await this.sketcher!.init(this);
      ui.setUpdateIndicator(this.host, false);
      molString ? this.setMolecule(molString) : this.sketcherCreated.next(true);
      this.changedSub = this.sketcher!.onChanged.subscribe((_: any) => {
        this.onChanged.next(null);
        for (let callback of this.listeners)
          callback();
        if (this.syncCurrentObject) {
          const molFile = this.getMolFile();
          grok.shell.o = SemanticValue.fromValueType(molFile, SEMTYPE.MOLECULE, UNITS.Molecule.MOLBLOCK);
        }
      });
    }
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
   * Returns the specified number of most diverse molecules in the column.
   * See example: {@link https://datagrok.ai/help/domains/chem/diversity-search}
   * @async
   * @param {Column} column - Column with molecules to search in
   * @param {Object} settings - Settings
   * @param {int} settings.limit - Would return top limit molecules
   * @returns {Promise<DataFrame>} - DataFrame with 1 column:
   *   - molecule: set of diverse structures
   * */
  export async function diversitySearch(column: Column, settings = {limit: Number.MAX_VALUE}): Promise<DataFrame> {
    const result = await grok.functions.call('Chem:getDiversities', {
      'molStringsColumn': column,
      'limit': settings.limit
    });
    return result;
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
  export async function searchSubstructure(column: Column, pattern: string = '', settings: {
    molBlockFailover?: string;
  } = {}): Promise<BitSet> {

    return (await grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn': column,
      'molString': pattern,
      'molBlockFailover': (settings.hasOwnProperty('molBlockFailover') ? settings.molBlockFailover : '') ?? ''
    })).get(0);
  }
 
  /**
   * Performs R-group analysis.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
   * @async
   * @param {DataFrame} table - Table.
   * @param {string} column - Column name with molecules to analyze.
   * @param {string} core - Core molecule.
   * @returns {Promise<DataFrame>}
   * */
  export async function rGroup(table: DataFrame, column: string, core: string): Promise<DataFrame> {
    return await grok.functions.call('Chem:FindRGroups', {
      column, table, core, prefix: 'R'});
  }

  /**
   * Finds Most Common Substructure in the specified column.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
   * @async
   * @param {Column} column - Column with SMILES to analyze.
   * @returns {Promise<string>}
   * */
   export async function mcs(column: Column): Promise<string> {
    return await grok.functions.call('Chem:searchSubstructure', {
      'molecules': column
    });
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

}