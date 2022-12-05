/**
 * Cheminformatics support
 * @module chem
 * */

import {BitSet, Column, DataFrame} from './dataframe';
import {FUNC_TYPES, SEMTYPE, SIMILARITY_METRIC, SimilarityMetric, UNITS} from './const';
import {Subject, Subscription} from 'rxjs';
import {Menu, Widget} from './widgets';
import {Func} from './entities';
import * as ui from '../ui';
import {SemanticValue} from './grid';
import $ from 'cash-dom';

let api = <any>window;
declare let grok: any;

const DEFAULT_SKETCHER = 'openChemLibSketcher';
export const WHITE_MOLBLOCK = `
  Datagrok empty molecule

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;
export const WHITE_MOLBLOCK_V_3000 = `Datagrok macromolecule handler
0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 0 0 0 0 0
M  V30 END CTAB
M  END
$$$$`

let extractors: Func[];  // id => molecule

export function isMolBlock(s: string | null) {
  return s != null && s.includes('M  END');
}

interface IStoredMolecule {
  molfile: string;
  smiles: string
}

/** Cheminformatics-related routines */
export namespace chem {

  export let SKETCHER_LOCAL_STORAGE = 'sketcher';
  export const STORAGE_NAME = 'sketcher';
  export const KEY = 'selected';
  export const SMILES = 'smiles';
  export const MOLV2000 = 'molv2000';
  export const MOLV3000 = 'molV3000';
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

    get molV3000(): string {
      return this.host!._molfile;
    }

    set molV3000(s: string) { }

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
    extractorsCreated = new Subject<boolean>();
    sketcherDialogOpened = false;

    /** Whether the currently drawn molecule becomes the current object as you sketch it */
    syncCurrentObject: boolean = true;

    listeners: Function[] = [];
    _mode = SKETCHER_MODE.INPLACE;
    _smiles = '';
    _molfile = WHITE_MOLBLOCK;
    _smarts = '';
    unitsBeforeInit = '';
    molFileUnits = MOLV2000;


    extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});
    inplaceSketcherDiv: HTMLDivElement | null = null;

    getSmiles(): string {
      return this.sketcher && this.sketcher._sketcher ? this.sketcher.smiles : this._smiles;
    }

    setSmiles(x: string): void {
      this._smiles = x;
      this.sketcher && this.sketcher._sketcher ? this.sketcher!.smiles = x : this.unitsBeforeInit = SMILES;
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    getMolFile(): string {
      return this.sketcher && this.sketcher._sketcher ?
        this.molFileUnits === MOLV2000 ? this.sketcher.molFile : this.sketcher.molV3000 : this._molfile;
    }

    setMolFile(x: string): void {
      this.molFileUnits = x && x.includes('V3000') ? MOLV3000 : MOLV2000;
      this._molfile = x;
      if (this.sketcher && this.sketcher._sketcher) {
        this.molFileUnits === MOLV2000 ? this.sketcher!.molFile = x : this.sketcher!.molV3000 = x;
      } else {
        this.unitsBeforeInit = this.molFileUnits;
      }
      this.updateExtSketcherContent(this.extSketcherDiv);
    }

    async getSmarts(): Promise<string | null> {
      let returnConvertedSmarts = async () => { // in case getter is called before sketcher initialized
        if (this._smiles) {
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
      };
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
      return Sketcher.isEmptyMolfile(molFile);
    }

    /** Sets the molecule, supports either SMILES or MOLBLOCK formats */
    async setMolecule(molString: string, substructure: boolean = false) {
      if (substructure)
        this.setSmarts(molString);
      else if (isMolBlock(molString))
        this.setMolFile(molString);
      else {
        this._smiles = molString;
        //setting molFile instead of smiles to draw in coordinates similar to dataframe cell
        const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molString, '{"mergeQueryHs":true}');
        if (!mol.has_coords())
          mol.set_new_coords();
        mol.normalize_depiction(1);
        mol.straighten_depiction(false);
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
    async setValue(x: string) {
      const extractor = extractors
        .find((f) => new RegExp(f.options['inputRegexp']).test(x));

      if (extractor != null)
        extractor
          .apply([new RegExp(extractor.options['inputRegexp']).exec(x)![1]])
          .then((mol) => this.setMolecule(mol));
      else
        this.setMolecule(x);
    }

    constructor(mode?: SKETCHER_MODE) {
      super(ui.div());
      if (mode)
        this._mode = mode;
      this.root.append(ui.div([ui.divText('')]));
      this.sketcherCreated.subscribe(() => {
        const molecule = this.unitsBeforeInit === SMILES ? this._smiles :
          (this.unitsBeforeInit === MOLV2000 || this.unitsBeforeInit === MOLV3000) ? this._molfile : this._smarts;
      if (molecule && molecule !== '' && !Sketcher.isEmptyMolfile(molecule))
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
      const lastSelecttedSketcher = window.localStorage.getItem(SKETCHER_LOCAL_STORAGE);
      this.sketcherFunctions = Func.find({tags: ['moleculeSketcher']});
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
        let canvas = this.createCanvas(width, height);
        canvas.style.height = '100%';
        canvas.style.width = '100%';
        ui.tooltip.bind(canvas, 'Click to edit');
        const clearButton = this.createClearSketcherButton(canvas);
        canvasMol(0, 0, width, height, canvas, this.getMolFile()!, null, {normalizeDepiction: true, straightenDepiction: true})
          .then((_) => {
            ui.empty(this.extSketcherDiv);
            this.extSketcherDiv.append(canvas);
            this.extSketcherDiv.append(clearButton);
          });
      }

      const sketchLinkStyle = {
        style: {
          width: `${width * 0.9}px`,
          height: `${height / 2}px`,
          textAlign: 'center',
          verticalAlign: 'middle',
          lineHeight: `${height / 2}px`,
          border: '1px solid #dbdcdf'
        }
      };
      let sketchLink = ui.divText('Click to edit', sketchLinkStyle);
      ui.tooltip.bind(sketchLink, 'Click to edit');
      sketchLink.style.paddingLeft = '0px';
      sketchLink.style.marginLeft = '0px';
      this._updateExtSketcherInnerHTML(sketchLink);
    };

    createClearSketcherButton(canvas: HTMLCanvasElement): HTMLButtonElement {
      const clearButton = ui.button('Clear', () => {
        this.setMolecule('');
        if (!this.sketcher) {
          this.onChanged.next(null);
        }
      });
      ui.tooltip.bind(clearButton, 'Clear sketcher');
      clearButton.style.position = 'absolute';
      clearButton.style.right = '0px';
      clearButton.style.fontSize = '10px';
      clearButton.style.visibility = 'hidden';
      clearButton.style.display = 'inline-flex';
      clearButton.onmouseover = () => {clearButton.style.visibility = 'visible';};
      canvas.onmouseenter = () => {clearButton.style.visibility = 'visible';};
      canvas.onmouseout = () => {clearButton.style.visibility = 'hidden';};
      return clearButton;
    }

    createExternalModeSketcher(): HTMLElement {
      this.extSketcherDiv = ui.div([], {style: {cursor: 'pointer'}});

      this.extSketcherDiv.onclick = () => {
        if (!this.sketcherDialogOpened) {
          this.sketcherDialogOpened = true;
          let savedMolFile = this.getMolFile();

          let dlg = ui.dialog();
          dlg.add(this.createInplaceModeSketcher(savedMolFile!))
            .onOK(() => {
              this.updateExtSketcherContent(this.extSketcherDiv);
              Sketcher.addRecent(savedMolFile!);
              this.sketcherDialogOpened = false;
            })
            .onCancel(() => {
              this.setMolFile(savedMolFile!);
              this.sketcherDialogOpened = false;
            })
            .show();
        }
      };

      ui.onSizeChanged(this.extSketcherDiv).subscribe((_) => {
        this.updateExtSketcherContent(this.extSketcherDiv);
      });

      this.updateExtSketcherContent(this.extSketcherDiv);
      return this.extSketcherDiv;
    }

    createInplaceModeSketcher(molStr?: string): HTMLElement {
      const molInputDiv = ui.div();
      $(this.molInput).attr('placeholder', 'SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');

      if (extractors == null) {
        grok.dapi.functions.filter('options.role="converter"').list()
          .then((res: Func[]) => {
            extractors = res.filter(it => it.outputs.filter(o => o.semType == SEMTYPE.MOLECULE).length);
          })
          .catch((_: any) => {
            extractors = [];
          });
      }

      const applyInput = async (e: any) => {
        const newSmilesValue: string = (e?.target as HTMLTextAreaElement).value;

        if (this.getSmiles() !== newSmilesValue)
          await this.setValue(newSmilesValue);

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
          .items(Sketcher.getRecent().map((m) => ui.tools.click(this.drawToCanvas(150, 60, m.molfile), () => this.setMolecule(m.molfile))), () => { })
          .endGroup()
          .group('Favorites')
          .item('Add to Favorites', () => Sketcher.addFavorite(this.getMolFile()))
          .separator()
          .items(Sketcher.getFavorites().map((m) => ui.tools.click(this.drawToCanvas(150, 60, m.molfile), () => this.setMolecule(m.molfile))), () => { })
          .endGroup()
          .separator()
          .items(this.sketcherFunctions.map((f) => f.friendlyName), (friendlyName: string) => {
              this.selectedSketcher = this.sketcherFunctions.filter(f => f.friendlyName === friendlyName)[0];
              grok.dapi.userDataStorage.postValue(STORAGE_NAME, KEY, this.selectedSketcher!.name, true);
              window.localStorage.setItem(SKETCHER_LOCAL_STORAGE, this.selectedSketcher!.name);
              this.setSketcher(this.getMolFile());
            },
            {isChecked: (item) => item === this.selectedSketcher?.friendlyName, toString: item => item})
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

    static getFavorites(): IStoredMolecule[] {
      return JSON.parse(localStorage.getItem(Sketcher.FAVORITES_KEY) ?? '[]');
    }

    static addFavorite(molecule: string) {
      const favorites = Sketcher.getRecent();
      Sketcher.checkDuplicatesAndAddToStorage(favorites, molecule, Sketcher.FAVORITES_KEY);
    }

    static getRecent(): IStoredMolecule[] {
      return JSON.parse(localStorage.getItem(Sketcher.RECENT_KEY) ?? '[]');
    }

    static addRecent(molecule: string) {
      const recent = Sketcher.getRecent();
      Sketcher.checkDuplicatesAndAddToStorage(recent, molecule, Sketcher.RECENT_KEY);
    }

    static checkDuplicatesAndAddToStorage(storage: IStoredMolecule[], molecule: string, localStorageKey: string) {
      let mol: any;
      grok.functions.call('Chem:getRdKitModule').then((rdKit: any) => {
        mol = rdKit.get_mol(molecule);
        const smiles = mol.get_smiles();         
        if (!storage.filter(mol => mol.smiles === smiles).length && !Sketcher.isEmptyMolfile(molecule)) {
          let s = JSON.stringify([...storage.slice(-9), {molfile: molecule, smiles: smiles}]);
          localStorage.setItem(localStorageKey, s);
        }
      }).finally(() => {
        mol.delete();
      });
    }

    static isEmptyMolfile(molFile: string): boolean {
      const rowWithAtomsAndNotation = molFile && molFile.split("\n").length >= 4 ? molFile.split("\n")[3] : '';
      return (molFile == null || molFile == '' ||
       (rowWithAtomsAndNotation.trimStart()[0] === '0' && rowWithAtomsAndNotation.trimEnd().endsWith('V2000')) ||
       (rowWithAtomsAndNotation.trimEnd().endsWith('V3000') && molFile.includes('COUNTS 0')));
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
          if (!Sketcher.isEmptyMolfile(molFile))
            grok.shell.o = SemanticValue.fromValueType(molFile, SEMTYPE.MOLECULE, UNITS.Molecule.MOLBLOCK);
        }
      });
    }

    drawToCanvas(w: number, h: number, molecule: string): HTMLElement{
      const imageHost = this.createCanvas(w, h);
      canvasMol(0, 0, w, h, imageHost, molecule, null, {normalizeDepiction: true, straightenDepiction: false});
      return imageHost;
    }

    createCanvas(w: number, h: number): HTMLCanvasElement {
      const imageHost = ui.canvas(w, h);
      const r = window.devicePixelRatio;
      imageHost.width = w * r;
      imageHost.height = h * r;
      imageHost.style.width = (w).toString() + 'px';
      imageHost.style.height = (h).toString() + 'px';
      return imageHost;
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
  export async function findSimilar(column: Column, molecule: string = '', settings = {
    limit: Number.MAX_VALUE,
    cutoff: 0.0
  }): Promise<DataFrame | null> {
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
      column, table, core, prefix: 'R'
    });
  }

  /**
   * Finds Most Common Substructure in the specified column.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
   * @async
   * @param {Column} column - Column with SMILES to analyze.
   * @returns {Promise<string>}
   * */
  export async function mcs(table: DataFrame, column: string, returnSmarts: boolean = false): Promise<string> {
    return await grok.functions.call('Chem:FindMCS', {
      'molecules': column,
      'df': table,
      'returnSmarts': returnSmarts
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
    options?: { [key: string]: boolean | number | string }
  ): HTMLDivElement {
    let root = document.createElement('div');
    // @ts-ignore
    import('openchemlib/full.js').then((OCL) => {
      let m = smiles.includes('M  END') ? OCL.Molecule.fromMolfile(smiles) : OCL.Molecule.fromSmiles(smiles);
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
    canvas: Object, molString: string, scaffoldMolString: string | null = null,
    options = {normalizeDepiction: true, straightenDepiction: true}
  ): Promise<void> {
    await grok.functions.call('Chem:canvasMol', {
      'x': x, 'y': y, 'w': w, 'h': h, 'canvas': canvas,
      'molString': molString, 'scaffoldMolString': scaffoldMolString ?? '',
      'options': options
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
}
