/* eslint-disable max-len */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {Editor} from 'ketcher-react';
import {StandaloneStructServiceProvider} from 'ketcher-standalone';
import {Ketcher} from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import '../css/editor.css';
import {KETCHER_MOLV2000, KETCHER_MOLV3000} from './constants';

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

export class KetcherSketcher extends grok.chem.SketcherBase {
  // ketcher-core is built on module-level singletons (CoreEditor, indigoWorker, ketcherProvider),
  // so only one live editor per page is possible (upstream limitation): mounting a new one
  // suspends all others behind a "Reload" placeholder instead of letting them silently break.
  private static _instances = new Set<KetcherSketcher>();
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

  constructor() {
    super();
    this.ketcherHost = ui.div([], 'ketcher-host');
    this.root.appendChild(this.ketcherHost);
    KetcherSketcher._instances.add(this);
    this._mountEditor();
  }

  private _mountEditor(): void {
    for (const other of KetcherSketcher._instances) {
      if (other !== this)
        other._suspend();
    }
    this._suspended = false;
    ui.empty(this.ketcherHost);

    const structServiceProvider = new StandaloneStructServiceProvider();

    const props = {
      staticResourcesUrl: !_package.webRoot ?
        '' :
        _package.webRoot.substring(0, _package.webRoot.length - 1),
      structServiceProvider: structServiceProvider,
      errorHandler: (message: string) => {
        console.log('Sketcher error', message);
      },
      onInit: (ketcher: Ketcher) => {
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
        this._restoreMolecule();
        (this._sketcher.editor as any).subscribe('change', async () => {
          if (this._detached || this._suspended)
            return;
          this.updatingMolecule = false;
          // we do not reset explicit mol in case this is the first change event called after ketcher was created
          // since change event is fired not only when user changes the molecule but also when the molecule is
          // initially set into ketcher
          if (this.importedMoleculesCounter > 0)
            this.importedMoleculesCounter --;
          else
            this.explicitMol = null;
          try {
            this._smiles = await this._sketcher!.getSmiles();
          } catch { //in case we are working with smarts - getSmiles() will fail with exception
            this._smiles = null;
          }
          // detach() (e.g. clicking OK on the cell editor dialog) unmounts the React <Editor>
          // while the awaits above are still pending; ketcher-core then drops its singleton
          // instance and any getMolfile() still in flight throws "couldnt find ketcher instance N".
          try {
            if (this._detached || this._suspended)
              return;
            this._molV2000 = await this._sketcher!.getMolfile(KETCHER_MOLV2000);
            this._molV3000 = await this._sketcher!.getMolfile(KETCHER_MOLV3000);
            this._smarts = await this._sketcher!.getSmarts();
          } catch {
            return;
          }
          this.onChanged.next(null);
        });
      },
    };

    this.reactRoot = ReactDOM.createRoot(this.ketcherHost);
    this.reactRoot.render(React.createElement(Editor, props, null));
  }

  private _suspend(): void {
    if (this._suspended || this._detached || this.reactRoot === null)
      return;
    this._suspended = true;
    try {
      this.reactRoot.unmount();
    } catch (e) {
      console.error(e);
    }
    this.reactRoot = null;
    this._sketcher = null;
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
    this.importedMoleculesCounter++;
    this._setNotation('smiles', smiles);
  }

  get molFile() {
    if (this.explicitMol?.notation === 'molblock')
      return this.explicitMol.value;
    if (this._molV2000 !== null) {
      if (this._molV3000 !== null && this._molV3000.includes('MDLV30/STE'))
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
    this.importedMoleculesCounter++;
    this._setNotation('molblock', molfile);
  }

  get molV3000() {
    if (this.explicitMol?.notation === 'molblockV3000')
      return this.explicitMol.value;
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
    this.importedMoleculesCounter++;
    this._setNotation('molblockV3000', molfile);
  }

  async getSmarts(): Promise<string> {
    if (this._sketcher)
      return !this._detached ? await this._sketcher.getSmarts() : this._smarts ?? '';
    return this._smarts ?? '';
  }

  set smarts(smarts: string) {
    this._smarts = smarts;
    this._molV3000 = null;
    this._molV2000 = null;
    this._smiles = null;
    this.importedMoleculesCounter++;
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
      this._sketcher?.setMolecule(molecule);
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
    this.updatingMolecule = true;
    this.setKetcherMolecule(value);
    //@ts-ignore
    this.explicitMol = {notation, value};
  }

  detach() {
    this._detached = true;
    KetcherSketcher._instances.delete(this);
    // grok.dapi.userDataStorage.postValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, JSON.stringify(this._sketcher?.editor.options()), true);
    this.reactRoot?.unmount();
    this.reactRoot = null;
    super.detach();
    //if detach occured while setting molecule into ketcher, send onChange, since we will not enter ketcher's onChange handler
    if (this.updatingMolecule)
      this.onChanged.next(null);
  }
}
