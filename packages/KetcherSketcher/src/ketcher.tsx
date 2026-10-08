/* eslint-disable max-len */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {Editor} from 'ketcher-react';
import {StandaloneStructServiceProvider} from 'ketcher-standalone';
import {Ketcher, MolSerializer, Pile, SettingsService, SupportedFormat} from 'ketcher-core';
import {asKetcherQuery, withQueries} from './query-molfile';
import 'ketcher-react/dist/index.css';
import '../css/editor.css';

type NotationKey = 'smiles' | 'molblock' | 'molblockV3000' | 'smarts';

// Ketcher's <RulerArea> reads SVGLength.value on a width="100%" canvas before the SVG can
// resolve relative units (e.g. while the macromolecules editor is still display:none),
// throwing NotSupportedError ("Could not resolve relative length") mid-render and blanking
// the editor. Return 0 in only that case — Ketcher already handles a 0-width canvas — while
// letting any unrelated error propagate. Workaround for upstream epam/ketcher#7568 / #3515;
// remove once Ketcher stops measuring the hidden canvas during initialization.
const _svgLenDesc = Object.getOwnPropertyDescriptor(SVGLength.prototype, 'value')!;
const _svgLenGet = _svgLenDesc.get!;
Object.defineProperty(SVGLength.prototype, 'value', {
  ..._svgLenDesc,
  get() {
    try {
      return _svgLenGet.call(this);
    } catch (e) {
      if (e instanceof DOMException)
        return 0;
      throw e;
    }
  },
});

/** How long an Indigo request is waited for before it is given up (`KetcherSketcher._inTurn`, `pageStructServiceProvider`). */
const INDIGO_ANSWER_MS = 15000;

/** `request`, once the requests before it have ended and `ready` has settled, given up after INDIGO_ANSWER_MS. */
function inTurnAfter<T>(before: Promise<unknown>, request: () => Promise<T>,
  ready: Promise<unknown> = Promise.resolve()): Promise<T> {
  return before.then(() => ready).then(() => new Promise<T>((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error('Indigo did not answer')), INDIGO_ANSWER_MS);
    request().then((value) => {
      clearTimeout(timer);
      resolve(value);
    }, (e) => {
      clearTimeout(timer);
      reject(e);
    });
  }));
}

/** The struct service's methods that send the worker nothing. */
const LOCAL_METHODS = new Set(['addKetcherId', 'getStandardServerOptions', 'callIndigoLoadedCallback',
  'callIndigoNoRenderLoadedCallback']);
let pageStructService: any = null;
let pageIndigoTurn: Promise<unknown> = Promise.resolve();

/** One struct service for every Ketcher of the page. ketcher-standalone runs Indigo in one worker per page and pairs a
 * reply with its request by the request's kind alone: a reply reaches only the IndigoService created last (each new one
 * takes over the worker's `onmessage`), so a request of a Ketcher mounted before was never answered; and of two requests
 * of a kind in flight, the reply for one drops the other unanswered. This service sends its requests one at a time, each
 * given up after INDIGO_ANSWER_MS once Indigo has first answered (its worker loads Indigo's WebAssembly first, which
 * takes what it takes); no Ketcher's unmount ends the page's worker (crux-sketch spike query-roundtrip, K3). */
const pageStructServiceProvider = {
  mode: 'standalone',
  createStructService(options: any): any {
    if (pageStructService === null) {
      const service: any = new StandaloneStructServiceProvider().createStructService(options);
      const loaded: Promise<unknown> = service.info().catch(() => {});
      pageStructService = new Proxy(service, {
        get(target, key) {
          const value = target[key];
          if (typeof value !== 'function')
            return value;
          if (key === 'destroy')
            return () => {};
          if (LOCAL_METHODS.has(key as string))
            return value.bind(target);
          return (...args: any[]) => {
            const result = inTurnAfter(pageIndigoTurn, () => value.apply(target, args), loaded);
            pageIndigoTurn = result.catch(() => {});
            return result;
          };
        },
      });
    }
    return pageStructService;
  },
};

/** ketcher-core's Ketcher subscribes to the page's one SettingsService (a singleton) in its constructor and never ends
 * that subscription, so every Ketcher ever mounted, its editor and its drawing stayed reachable from it: a page that had
 * mounted Ketcher some 50 to 150 times (dialogs opened, sketchers switched) filled its heap until it stopped answering
 * (crux-sketch spike query-roundtrip, K4; ketcher-core's own warning: "MaxListenersExceededWarning: Possible EventEmitter
 * memory leak detected. 11 settings:changed listeners"). Each subscription made while a sketcher mounts is its own, and
 * ends when the sketcher is suspended or detached (`_endSettingsSubscriptions`). */
const settingsSubscribe = SettingsService.prototype.subscribe;
SettingsService.prototype.subscribe = function(this: SettingsService, listener: any): () => void {
  const unsubscribe = settingsSubscribe.call(this, listener);
  KetcherSketcher._mounting?._settingsSubscriptions.push(unsubscribe);
  return unsubscribe;
};

export class KetcherSketcher extends grok.chem.SketcherBase {
  // ketcher-core is built on module-level singletons (CoreEditor, indigoWorker, ketcherProvider),
  // so only one live editor per page is possible (upstream limitation): mounting a new one
  // suspends all others behind a "Reload" placeholder instead of letting them silently break.
  private static _instances = new Set<KetcherSketcher>();
  private static _indigoTurn: Promise<unknown> = Promise.resolve();
  /** The sketcher mounting its editor, whose Ketcher's settings subscriptions are its own (K4). */
  static _mounting: KetcherSketcher | null = null;
  _settingsSubscriptions: (() => void)[] = [];
  _smiles: string | null = null;
  _molV2000: string | null = null;
  _molV3000: string | null = null;
  _smarts: string | null = null;
  _sketcher: Ketcher | null = null;
  ketcherHost: HTMLDivElement;
  reactRoot: ReactDOM.Root | null = null;
  updatingMolecule = false;
  private importedMoleculesCounter = 0;
  private _detached = false;
  private _suspended = false;
  private _exportId = 0;
  /** Settles once Ketcher's editor reports ready (its first onInit), or the sketcher is suspended or detached first:
   * `init` waits for it, so the host announces this sketcher ready (`sketcherReady`) only when it takes input. */
  private readonly _ready: Promise<void>;
  private _resolveReady: (() => void) | null = null;
  /** The last molecule set into Ketcher, settled once Ketcher has loaded it. */
  private _loading: Promise<unknown> = Promise.resolve();

  constructor() {
    super();
    this._ready = new Promise<void>((resolve) => this._resolveReady = resolve);
    this.ketcherHost = ui.div([], 'ketcher-host');
    this.root.appendChild(this.ketcherHost);
    KetcherSketcher._instances.add(this);
    this._mountEditor();
  }

  private _settleReady(): void {
    this._resolveReady?.();
    this._resolveReady = null;
  }

  private _mountEditor(): void {
    for (const other of KetcherSketcher._instances) {
      if (other !== this)
        other._suspend();
    }
    this._suspended = false;
    ui.empty(this.ketcherHost);
    KetcherSketcher._mounting = this;

    const props = {
      staticResourcesUrl: !_package.webRoot ?
        '' :
        _package.webRoot.substring(0, _package.webRoot.length - 1),
      structServiceProvider: pageStructServiceProvider as any,
      // its toggle is hidden (editor.css), and onInit would wait for its lazy chunk while the canvas already takes strokes
      disableMacromoleculesEditor: true,
      errorHandler: (message: string) => {
        console.log('Sketcher error', message);
      },
      onInit: (ketcher: Ketcher) => {
        // the Editor calls it again for the same Ketcher once the macromolecules editor it still mounts has loaded
        if (ketcher === this._sketcher)
          return;
        this._sketcher = ketcher;
        // workaround for sketcher not to be truncated when showed in a popup menu
        // in the end of the screen (on last dataframe column)
        if (this.host && this.host.isInPopupContainer()) {
          const ketcherRoot = this.ketcherHost.querySelector('.Ketcher-root');
          if (ketcherRoot)
            (ketcherRoot as HTMLElement).style.minWidth = '0px';
          this.ketcherHost.style.width = '100%';
        }
        // grok.dapi.userDataStorage.getValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, true).then((opts: string) => {
        //   if (opts) {
        //     this._sketcher?.editor.setOptions(opts);
        //   }
        // });
        (ketcher.editor as any).subscribe('change', () => {
          if (this._detached || this._suspended)
            return;
          this.updatingMolecule = false;
          // a molecule loaded into Ketcher answers with a change event of its own, which is not the user's edit
          if (this.importedMoleculesCounter > 0)
            this.importedMoleculesCounter --;
          else
            this.explicitMol = null;
          this._exportChange(this._sketcher!);
        });
        // Ketcher takes strokes before it reports ready: a drawing by then is the user's, and restoring the host's
        // molecule would load over it
        if (!KetcherSketcher._hasDrawing(ketcher))
          this._restoreMolecule();
        else {
          this.explicitMol = null;
          this._exportChange(ketcher);
        }
        this._settleReady();
      },
    };

    this.reactRoot = ReactDOM.createRoot(this.ketcherHost);
    this.reactRoot.render(React.createElement(Editor, props, null));
  }

  /** The template tool's floating preview sits in the structure as atoms too, and is not a drawing. */
  private static _hasDrawing(ketcher: Ketcher): boolean {
    const struct = ketcher.editor.struct();
    return [...struct.atoms.values()].some((a) => !a.isPreview) || struct.rxnArrows.size > 0 || struct.texts.size > 0;
  }

  /** Runs a conversion of this adapter's once the ones before it have ended: the standalone struct service hands a
   * worker reply to every pending conversion of the same input, whatever its format, and fails them all when the reply
   * is an error (SMILES of a query). One left unanswered is given up after INDIGO_ANSWER_MS, so that it holds up no
   * later conversion of the page's Ketchers: the queue stalled for good before, every SMARTS and V3000 export after it
   * never coming (crux-sketch spike query-roundtrip, K2; why one went unanswered: `pageStructServiceProvider`, K3). */
  private static _inTurn<T>(conversion: () => Promise<T>): Promise<T> {
    const result = inTurnAfter(KetcherSketcher._indigoTurn, conversion);
    KetcherSketcher._indigoTurn = result.catch(() => {});
    return result;
  }

  /** Exports the drawing as it is at the change. While the template tool hovers the canvas, Ketcher keeps the
   * template under the cursor in the structure itself (with no change event), so an export read later took it
   * in as a second copy of the template. The V2000 molblock is written at once and announced, so a dialog's
   * OK right after a stroke has it; the Indigo notations follow in turn and are announced again only when
   * they change the molblock (enhanced stereo). A conversion that starts after a detach finds no Ketcher
   * instance and fails, and the getters then convert the molblock instead. */
  private _exportChange(ketcher: Ketcher): void {
    if ((window as any).isPolymerEditorTurnedOn || ketcher.containsReaction())
      return;
    const struct = ketcher.editor.struct();
    const drawn = struct.clone(
      new Pile<number>(struct.atoms.keys()).filter((id) => !struct.atoms.get(id)!.isPreview),
      new Pile<number>(struct.bonds.keys()).filter((id) => !struct.bonds.get(id)!.isPreview));
    let molV2000: string;
    try {
      molV2000 = new MolSerializer().serialize(drawn);
    }
    catch {
      return;
    }
    // what ketcher-core's V2000 leaves out of a query, its SMARTS-only atom properties and its custom bond queries, put
    // back where RDKit reads them: V3000 with SMARTSQ groups where an atom has one (query-molfile.ts)
    try {
      molV2000 = withQueries(drawn, molV2000);
    } catch (e) {
      console.error(e);
    }
    const exportId = ++this._exportId;
    this._molV2000 = molV2000;
    this._smiles = this._molV3000 = this._smarts = null;
    this.onChanged.next(null);
    const molFile = this.molFile;
    const convert = (format: SupportedFormat) => KetcherSketcher._inTurn(() => exportId !== this._exportId ?
      Promise.reject(new Error('a later change is exported')) :
      ketcher.formatterFactory.create(format, ketcher.editor.serverSettings as any).getStringFromStructureAsync(drawn));
    Promise.all([convert(SupportedFormat.molV3000), convert(SupportedFormat.smarts), convert(SupportedFormat.smiles).catch(() => null)])
      .then(([molV3000, smarts, smiles]) => {
        if (exportId !== this._exportId)
          return;
        this._molV3000 = molV3000;
        this._smarts = smarts;
        this._smiles = smiles;
        if (this.molFile !== molFile)
          this.onChanged.next(null);
      }, () => {});
  }

  private _suspend(): void {
    if (this._suspended || this._detached || this.reactRoot === null)
      return;
    this._suspended = true;
    this._settleReady();
    try {
      this.reactRoot.unmount();
    } catch (e) {
      console.error(e);
    }
    this.reactRoot = null;
    this._sketcher = null;
    this._endSettingsSubscriptions();
    if (this.updatingMolecule) {
      this.updatingMolecule = false;
      this.onChanged.next(null);
    }
    ui.empty(this.ketcherHost);
    this.ketcherHost.appendChild(ui.divV([
      ui.divText('This sketcher was paused because another Ketcher sketcher was opened. ' +
        'Ketcher supports only one active editor per page.'),
      ui.button('Reload', () => this._mountEditor()),
    ], 'ketcher-suspended'));
  }

  /** Ends the settings subscriptions this sketcher's Ketchers made (K4): an unmounted Ketcher is then unreachable. */
  private _endSettingsSubscriptions(): void {
    for (const unsubscribe of this._settingsSubscriptions.splice(0))
      unsubscribe();
    if (KetcherSketcher._mounting === this)
      KetcherSketcher._mounting = null;
  }

  private _restoreMolecule(): void {
    if (this._smiles === null && this._smarts !== null) {
      this.smarts = this._smarts;
      return;
    }
    const mol = this.molFile;
    if (mol)
      this.molFile = mol;
    else
      this.setMoleculeFromHost();
  }

  async init(host: grok.chem.Sketcher) {
    this.host = host;
    if (this.host.isResizing)
      this.ketcherHost.classList.add('ketcher-resizing');
    await this._ready;
  }

  get supportedExportFormats() {
    return ['smiles', 'mol', 'smarts'];
  }

  get smiles() {
    if (this.explicitMol?.notation === 'smiles')
      return this.explicitMol.value;
    if (this._smiles !== null)
      return this._smiles;
    if (this._molV3000 !== null)
      return DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
    if (this._molV2000 !== null)
      return DG.chem.convert(this._molV2000, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
    if (this._smarts !== null)
      return DG.chem.smilesFromSmartsWarning();
    return '';
  }

  set smiles(smiles: string) {
    //in case we opened sketcher in filter, draw something and clicked Cancel -> ketcher will be detached
    // and we will not get into onChange event, so update inner structures to prevent loosing the information
    this._smiles = smiles;
    this._molV2000 = null;
    this._molV3000 = null;
    this._smarts = null;
    this._setNotation('smiles', smiles);
  }

  get molFile() {
    if (this.explicitMol?.notation === 'molblock')
      return this.explicitMol.value;
    if (this._molV2000 !== null) {
      // a query written V3000 (query-molfile.ts) is the molblock as it is
      if (this._molV3000 !== null && this._molV3000.includes('MDLV30/STE') && !this._molV2000.includes('V3000'))
        return DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.MolBlock);
      return this._molV2000;
    }
    if (this._molV3000 !== null)
      return DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.MolBlock);
    if (this._smiles !== null)
      return DG.chem.convert(this._smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    if (this._smarts !== null)
      return DG.chem.convert(this._smarts, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
    return '';
  }

  set molFile(molfile: string) {
    this._molV2000 = molfile;
    this._smiles = null;
    this._molV3000 = null;
    this._smarts = null;
    this._setNotation('molblock', molfile);
  }

  get molV3000() {
    if (this.explicitMol?.notation === 'molblockV3000')
      return this.explicitMol.value;
    // a query written V3000 with its SMARTSQ groups (query-molfile.ts), which Indigo's V3000 leaves out
    if (this._molV2000?.includes('V3000'))
      return this._molV2000;
    if (this._molV3000 !== null)
      return this._molV3000;
    if (this._molV2000 !== null)
      return DG.chem.convert(this._molV2000, DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock);
    if (this._smiles !== null)
      return DG.chem.convert(this._smiles, DG.chem.Notation.Smiles, DG.chem.Notation.V3KMolBlock);
    if (this._smarts !== null)
      return DG.chem.convert(this._smarts, DG.chem.Notation.Smarts, DG.chem.Notation.V3KMolBlock);
    return '';
  }

  set molV3000(molfile: string) {
    this._molV3000 = molfile;
    this._molV2000 = null;
    this._smiles = null;
    this._smarts = null;
    this._setNotation('molblockV3000', molfile);
  }

  async getSmarts(): Promise<string> {
    // the caller's SMARTS until the user's first edit, as the other getters give the caller's string
    if (this.explicitMol?.notation === 'smarts')
      return this.explicitMol.value;
    const ketcher = this._sketcher;
    if (ketcher) {
      // a molecule set loads asynchronously: export what it drew, not the canvas before it
      await this._loading;
      return !this._detached && this._sketcher === ketcher ?
        await KetcherSketcher._inTurn(() => ketcher.getSmarts()) : this._smarts ?? '';
    }
    return this._smarts ?? '';
  }

  set smarts(smarts: string) {
    this._smarts = smarts;
    this._molV3000 = null;
    this._molV2000 = null;
    this._smiles = null;
    this._setNotation('smarts', smarts);
  }

  get isInitialized() {
    return this._sketcher !== null;
  }

  resize() {
    this.ketcherHost.classList.add('ketcher-resizing');
  }

  setKetcherMolecule(molecule: string) {
    try {
      // a query molfile with SMARTSQ groups is shown with them as Ketcher's own query properties (query-molfile.ts)
      let shown = molecule;
      try {
        shown = asKetcherQuery(molecule);
      } catch (e) {
        console.error(e);
      }
      const loading = this._sketcher?.setMolecule(shown);
      if (loading)
        this._loading = loading.catch((e) => console.error(e));
    } catch (e) {
      console.error(e);
    }
  }

  setMoleculeFromHost(): void {
    const host = this.host;
    if (!host) return;
    if (host._molfile !== null) {
      if (host.molFileUnits === DG.chem.Notation.MolBlock)
        this.molFile = host._molfile;
      if (host.molFileUnits === DG.chem.Notation.V3KMolBlock)
        this.molV3000 = host._molfile;
      return;
    }
    if (host._smiles !== null) {
      this.smiles = host._smiles;
      return;
    }
    if (host._smarts !== null) {
      this.smarts = host._smarts;
      return;
    }
  }

  private _setNotation(notation: NotationKey, value: string): void {
    this._exportId++;
    //@ts-ignore
    this.explicitMol = {notation, value};
    const ketcher = this._sketcher;
    // An empty molecule changes nothing on an empty canvas, but its load is asynchronous: a stroke made meanwhile
    // was taken for the load's change event, and the load then wiped it.
    if (ketcher !== null && (!value?.trim() || grok.chem.Sketcher.isEmptyMolfile(value)) && !KetcherSketcher._hasDrawing(ketcher))
      return;
    this.updatingMolecule = true;
    if (ketcher === null)
      return;
    this.importedMoleculesCounter++;
    this.setKetcherMolecule(value);
  }

  detach() {
    this._detached = true;
    this._settleReady();
    KetcherSketcher._instances.delete(this);
    // grok.dapi.userDataStorage.postValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, JSON.stringify(this._sketcher?.editor.options()), true);
    this.reactRoot?.unmount();
    this.reactRoot = null;
    this._endSettingsSubscriptions();
    super.detach();
    //if detach occured while setting molecule into ketcher, send onChange, since we will not enter ketcher's onChange handler
    if (this.updatingMolecule)
      this.onChanged.next(null);
  }
}
