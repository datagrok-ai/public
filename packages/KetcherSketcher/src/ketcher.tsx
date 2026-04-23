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

export class KetcherSketcher extends grok.chem.SketcherBase {
  _smiles: string | null = null;
  _molV2000: string | null = null;
  _molV3000: string | null = null;
  _smarts: string | null = null;
  _sketcher: Ketcher | null = null;
  ketcherHost: HTMLDivElement;
  reactRoot: ReactDOM.Root | null = null;

  constructor() {
    super();
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
        if (this.host && this.host.isInPopupContainer()) { //workaround for sketcher not to be truncated when showed in a popup menu in the end of the screen (on last dataframe column)
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
        this.setMoleculeFromHost();
        (this._sketcher.editor as any).subscribe('change', async () => {
          // we do not reset explicit mol in case this is the first change event called after ketcher was created
          // since change event is fired not only when user changes the molecule but also when the molecule is
          // initially set into ketcher
          if (this._smiles !== null || this._molV2000 !== null || this._molV3000 !== null || this._smarts !== null)
            this.explicitMol = null;
          try {
            this._smiles = await this._sketcher!.getSmiles();
          } catch { //in case we are working with smarts - getSmiles() will fail with exception
            this._smiles = null;
          }
          this._molV2000 = await this._sketcher!.getMolfile(KETCHER_MOLV2000);
          this._molV3000 = await this._sketcher!.getMolfile(KETCHER_MOLV3000);
          this.onChanged.next(null);
        });
      },
    };

    this.ketcherHost = ui.div([], 'ketcher-host');

    const component = React.createElement(Editor, props, null);
    this.reactRoot = ReactDOM.createRoot(this.ketcherHost);
    this.reactRoot.render(component);

    this.root.appendChild(this.ketcherHost);
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
    if (this._molV2000 !== null)
      return DG.chem.convert(this._molV2000, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
    if (this._molV3000 !== null)
      return DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
    if (this._smarts !== null)
      return DG.chem.smilesFromSmartsWarning();
    return '';
  }

  set smiles(smiles: string) {
    this._setNotation('smiles', smiles);
  }

  get molFile() {
    if (this.explicitMol?.notation === 'molblock')
      return this.explicitMol.value;
    if (this._molV2000 !== null)
      return this._molV2000;
    if (this._molV3000 !== null)
      return DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.MolBlock);
    if (this._smiles !== null)
      return DG.chem.convert(this._smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    if (this._smarts !== null)
      return DG.chem.convert(this._smarts, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
    return '';
  }

  set molFile(molfile: string) {
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
    this._setNotation('molblockV3000', molfile);
  }

  async getSmarts(): Promise<string> {
    return this._sketcher ? await this._sketcher.getSmarts() : this._smarts ?? '';
  }

  set smarts(smarts: string) {
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
    this.setKetcherMolecule(value);
    if (notation !== 'smarts')
      this.explicitMol = {notation, value};
  }

  detach() {
    // grok.dapi.userDataStorage.postValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, JSON.stringify(this._sketcher?.editor.options()), true);
    this.reactRoot?.unmount();
    this.reactRoot = null;
    super.detach();
  }
}
