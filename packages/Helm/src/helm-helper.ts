import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

/* eslint-disable max-len */
import {
  App, Point, HelmType, IHelmBio, HelmAtom, HelmBond, HelmMol, GetMonomerFunc, GetMonomerResType,
  MonomerExplorer, TabDescType, MonomersFuncs, MonomerSetType, ISeqMonomer, PolymerType, MonomerType,
  DojoType, JSDraw2ModuleType, IHelmEditorOptions, IHelmDrawOptions,
  Editor,
  IBio,
} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HelmConvertRes, HelmInputBase, HelmNotSupportedError, IHelmHelper, IHelmInputInitOptions
} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLibBase, IMonomerLinkData, IMonomerSetPlaceholder} from '@datagrok-libraries/bio/src/types/monomer-library';
import {GAP_SYMBOL, GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {OrgHelmModule, ScilModule} from './types';
import {HelmInput} from './widgets/helm-input';
import {getHoveredMonomerFromEditorMol} from './utils/get-hovered';
/* eslint-enable max-len */

import {_package} from './package';
import {IEditorOptions} from '@datagrok-libraries/js-draw-lite/src/types/jsdraw2';

// Phase 3 (hwe migration): the interactive editor dialog is backed by the
// standalone `@datagrok-libraries/hwe` library through its Datagrok adapter.
// Import from the bare entry (the adapter is re-exported there in addition to
// the `/datagrok` sub-entry) so the Datagrok toolchain — ts-loader with classic
// `node` moduleResolution — resolves it without `exports`-map subpath support.
import {HelmHelperAdapter, bridgeMonomerLib, HelmService} from '@datagrok-libraries/hwe';
import type {IMonomerLibBaseLike} from '@datagrok-libraries/hwe';

declare const dojo: DojoType;
declare const JSDraw2: JSDraw2ModuleType;
declare const scil: ScilModule;
declare const org: OrgHelmModule;

const HELM_GAP_SYMBOL: string = GapOriginals[NOTATION.HELM];

export class HelmHelper implements IHelmHelper {
  private static instanceCount: number = 0;

  constructor(
    public readonly seqHelper: ISeqHelper,
    private readonly logger: ILogger
  ) {
    // Watchdog for singleton
    if ((++HelmHelper.instanceCount) > 1)
      throw new Error(`HelmHelper must be a single.`);
  }

  createHelmInput(name?: string, options?: IHelmInputInitOptions): HelmInputBase {
    try {
      const monomerLib: IMonomerLibBase = _package._libHelper!.getMonomerLib();
      return HelmInput.create(this, monomerLib, name, options);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(`Helm: HelmHelper.createHelmInput(), Error: ${errMsg}`, undefined, errStack);
      throw err;
    }
  }

  // Phase 6 (hwe migration): the inline view-only editor (used by the
  // `HelmInput` widget, the Bio HELM substructure-filter, and PolyTool) is now
  // the standalone `@datagrok-libraries/hwe` editor wrapped in a
  // `LegacyEditorWrapper` (presents the legacy `editor.div/m/helm/setHelm/
  // setMol/getHelm/redraw/setSize/getDrawOptions` shape over `HelmService`).
  createHelmWebEditor(host?: HTMLDivElement, options?: Partial<IHelmEditorOptions>): IHelmWebEditor {
    return this.editorAdapter.createHelmWebEditor(host, {
      ...(options?.width !== undefined ? {width: options.width} : {}),
      ...(options?.height !== undefined ? {height: options.height} : {}),
    }) as unknown as IHelmWebEditor;
  }

  // -- HWE editor adapter (Phase 3 migration) --
  // The full editor dialog is now the standalone `@datagrok-libraries/hwe`
  // editor. Its `HelmService` reads Bio's monomer library via `bridgeMonomerLib`
  // so the editor resolves exactly the same monomers as the rest of the
  // platform. Built lazily on first `createWebEditorApp` (after `completeInit`
  // has populated `_package._libHelper`).
  private _editorAdapter: HelmHelperAdapter | null = null;
  private get editorAdapter(): HelmHelperAdapter {
    if (this._editorAdapter === null) {
      const monomerLib = _package._libHelper!.getMonomerLib() as unknown as IMonomerLibBaseLike;
      const service = new HelmService({monomerLib: bridgeMonomerLib(monomerLib)});
      this._editorAdapter = new HelmHelperAdapter({service, seqHelper: this.seqHelper});
    }
    return this._editorAdapter;
  }

  /**
   * Phase 3 (hwe migration): returns the standalone `@datagrok-libraries/hwe`
   * editor (palette + toolbar + canvas + Sequence/HELM/Properties/Atomic tabs)
   * wrapped in a `LegacyAppWrapper` that presents the legacy Pistoia `App`
   * shape — `canvas.getHelm(true)`, `canvas.helm.jsd.m.atoms[i].selected`,
   * `sequence` / `notation` / `properties.parent`, `mex`, `toolbarheight`,
   * `calculateSizes()` — so the existing dialog consumers swap in unchanged.
   * The hwe app self-lays-out (no manual size dance), and the legacy
   * Datagrok-only "Placeholders" monomer-explorer tab is dropped (hwe ships
   * its own palette).
   */
  createWebEditorApp(host: HTMLDivElement, helm?: string): App {
    return this.editorAdapter.createWebEditorApp(host, helm) as unknown as App;
  }

  // -- MonomerExplorer.onShowTab --

  public onShowTabPlaceholders(mex: MonomerExplorer, div: HTMLDivElement): void {
    const d = scil.Utils.createElement(div, 'div', null, {paddingTop: '5px'});

    /* eslint-disable max-len*/
    // if (this.options.canvastoolbar == false) {
    //   const b = scil.Utils.createElement(d, 'div', '<img src=\'' + scil.Utils.imgSrc('helm/arrow.png') + '\' style=\'vertical-align:middle\'>Mouse Pointer', {cursor: 'pointer', padding: '2px', border: 'solid 1px gray', margin: '5px'});
    //   scil.connect(b, 'onclick', function() {
    //     mex.plugin.jsd.doCmd('lasso');
    //   });
    // }
    /* eslint-enable max-len*/

    const phTabs: Partial<TabDescType>[] = [];
    phTabs.push({caption: 'Meta', tabkey: 'placeholders-meta'});
    phTabs.push({caption: 'Sets', tabkey: 'placeholders-sets'});
    const placeholdersTabs = new scil.Tabs(d, {
      onShowTab(td: HTMLTableCellElement) {
        mex.onShowTab(td);
      },
      tabpadding: '5px 2px 1px 2px',
      tabs: phTabs,
      marginBottom: `0`,
    });
  }

  public onShowTabPlaceholdersMeta(mex: MonomerExplorer, div: HTMLDivElement): void {

  }

  public onShowTabPlaceholdersSets(mex: MonomerExplorer, div: HTMLDivElement): void {
    const monomerSets = _package._libHelper!.getMonomerSets();

    const phSet: { [polymerType: string]: { [helmType: string]: IMonomerSetPlaceholder[] } } = {};
    for (const ph of monomerSets.placeholders) {
      let ofPolymerTypeSet = phSet[ph.polymerType];
      if (!ofPolymerTypeSet) ofPolymerTypeSet = phSet[ph.polymerType] = {};

      const phHelmType = this.toHelmType(ph.polymerType, ph.monomerType);
      let phList = ofPolymerTypeSet[phHelmType];
      if (!phList) phList = ofPolymerTypeSet[phHelmType] = [];
      phList.push(ph);
    }

    const fillPlaceholders = (
      parentDiv: HTMLDivElement, phList: IMonomerSetPlaceholder[], helmType: HelmType
    ): void => {
      for (const ph of phList) {
        const phDiv = mex.createMonomerDiv(parentDiv, ph.symbol, helmType);
        ui.tooltip.bind(phDiv, () => {
          const phTooltipDiv = ui.div([
            ui.table(ph.monomerLinks, (ml: IMonomerLinkData, idx: number) => {
              return [ml.symbol, ml.source];
            }, ['Symbol', 'Library source'])
          ]);
          return phTooltipDiv;
        });
      }
    };

    const polymerTypes: [string, PolymerType][] = [
      ['Chem', PolymerTypes.CHEM],
      ['Peptide', PolymerTypes.PEPTIDE],
      ['RNA', PolymerTypes.RNA],
      ['Blob', PolymerTypes.BLOB],
      ['G', PolymerTypes.G],
    ];
    const polymerTypeAcc = ui.accordion();
    for (const [name, polymerType] of polymerTypes) {
      if (polymerType in phSet) {
        polymerTypeAcc.addPane(name, () => {
          const polymerTypeDiv = ui.div();
          const ofPolymerTypeSet = phSet[polymerType];
          if (Object.keys(ofPolymerTypeSet).length == 1) {
            const [helmType, phList] = Object.entries(ofPolymerTypeSet)[0];
            fillPlaceholders(polymerTypeDiv, phList, helmType as HelmType);
          } else {
            const helmTypeAcc = ui.accordion();
            for (const [helmTypeName, helmType] of Object.entries(HelmTypes)) {
              if (helmType in ofPolymerTypeSet) {
                helmTypeAcc.addPane(helmTypeName, () => {
                  const helmTypeDiv = ui.div();

                  const phList = ofPolymerTypeSet[helmType];
                  fillPlaceholders(helmTypeDiv, phList, helmType as HelmType);

                  return helmTypeDiv;
                });
              }
            }
            polymerTypeDiv.appendChild(helmTypeAcc.root);
          }
          return polymerTypeDiv;
        }, true, undefined, false);
      }
    }
    div.appendChild(polymerTypeAcc.root);
  }

  // -- GetMonomer --

  public originalMonomersFuncs: MonomersFuncs | null = null;

  public revertOriginalMonomersFuncs(): MonomersFuncs {
    const logPrefix = `Helm: revertOriginalGetMonomer()`;
    if (this.originalMonomersFuncs === null)
      throw new Error('Unable to revert original getMonomer');

    const monomers = org.helm.webeditor.Monomers;
    const overriddenMonomersFuncs: MonomersFuncs = {
      getMonomer: monomers.getMonomer,
      getMonomerSet: monomers.getMonomerSet,
    };

    monomers.getMonomer = this.originalMonomersFuncs.getMonomer;
    monomers.getMonomerSet = this.originalMonomersFuncs.getMonomerSet;
    this.originalMonomersFuncs = null;

    return overriddenMonomersFuncs;
  }

  public overrideMonomersFuncs(monomersFuncs: MonomersFuncs): MonomersFuncs {
    const logPrefix = `Helm: HelmHelper.overrideGetMonomer()`;

    if (this.originalMonomersFuncs != null)
      throw new Error(`${logPrefix}, originalGetMonomer is overridden already`);

    const monomers = org.helm.webeditor.Monomers;
    this.originalMonomersFuncs = {
      getMonomer: monomers.getMonomer.bind(monomers),
      getMonomerSet: monomers.getMonomerSet.bind(monomers),
    };
    monomers.getMonomer = monomersFuncs.getMonomer;
    monomers.getMonomerSet = monomersFuncs.getMonomerSet;
    return this.originalMonomersFuncs!;
  }

  public buildMonomersFuncsFromLib(monomerLib: IMonomerLib): MonomersFuncs {
    const logPrefix = `Helm: HelmHelper.buildGetMonomerFromLib()`;

    return {
      getMonomer: (a: HelmAtom | HelmType, name?: string): GetMonomerResType => {
        const logPrefixInt = `${logPrefix}, org.helm.webeditor.Monomers.getMonomer()`;
        try {
          // logger.debug(`${logPrefixInt}, a: ${JSON.stringify(a, helmJsonReplacer)}, name: '${name}'`);
          if (name?.startsWith('[') && name.endsWith(']'))
            name = name.substring(1, name.length - 1);
          // Creates monomers in lib
          const dgWem = monomerLib.getWebEditorMonomer(a, name);

          return dgWem; //dgWem;
        } catch (err) {
          const [errMsg, errStack] = errInfo(err);
          this.logger.error(`${logPrefixInt}, Error: ${errMsg}`, undefined, errStack);
          throw err;
        }
      },
      getMonomerSet: (a: HelmAtom | HelmType | null): MonomerSetType | null => {
        // if (a == null)
        //   return null;
        // const atom = a as HelmAtom;
        // if (atom.T === 'ATOM')
        //   a = atom.biotype();
        // return monomerLib.getMonomerSet(a as HelmType);
        return this.originalMonomersFuncs!.getMonomerSet(a);
      }
    };
  }

  // helm2type(seqM: ISeqMonomer): HelmType {
  //   if (seqM.polymerType == "PEPTIDE")
  //     return org.helm.webeditor.HELM.AA;
  //   else if (seqM.polymerType == "CHEM")
  //     return org.helm.webeditor.HELM.CHEM;
  //   else if (seqM.polymerType == "RNA") {
  //     if (seqM.mt == "Branch")
  //       return org.helm.webeditor.HELM.BASE;
  //     if (seqM.mt == "Backbone") {
  //       if (seqM.na == "P" || seqM.na == "p")
  //         return org.helm.webeditor.HELM.LINKER;
  //       else
  //         return org.helm.webeditor.HELM.SUGAR;
  //     }
  //   }
  //   throw new Error(`Not supported seq monomer ${JSON.stringify(seqM)}.`);
  // }
  private toHelmType(polymerType: PolymerType, monomerType: MonomerType, naturalAnalog?: string) {
    if (polymerType == PolymerTypes.PEPTIDE)
      return org.helm.webeditor.HELM.AA;
    else if (polymerType == PolymerTypes.CHEM)
      return org.helm.webeditor.HELM.CHEM;
    else if (polymerType == PolymerTypes.RNA) {
      if (monomerType == MonomerTypes.BRANCH)
        return org.helm.webeditor.HELM.BASE;
      if (monomerType == MonomerTypes.BACKBONE) {
        if (naturalAnalog == 'P' || naturalAnalog == 'p')
          return org.helm.webeditor.HELM.LINKER;
        else
          return org.helm.webeditor.HELM.SUGAR;
      }
    }
    throw new Error(`Not supported monomer polymerType: '${polymerType}', monomerType: '${monomerType}'.`);
  }

  public getHoveredAtom(x: number, y: number, mol: HelmMol, height: number): HelmAtom | null {
    return getHoveredMonomerFromEditorMol(x, y, mol, height);
  }

  //if molfiles are accessed too often, it is better to cache it
  private _molfilesEditor: Editor<unknown, IBio<unknown>, IEditorOptions> | null = null;
  private _molfilesEditorHost: HTMLDivElement | null = null;
  private get molfilesEditor(): Editor<unknown, IBio<unknown>, IEditorOptions> {
    if (this._molfilesEditorRemovingTimer)
      clearTimeout(this._molfilesEditorRemovingTimer);
    if (this._molfilesEditor == null) {
      this._molfilesEditorHost = ui.div([], {style: {width: '0px', height: '0px', position: 'fixed'}});
      document.documentElement.appendChild(this._molfilesEditorHost);
      this._molfilesEditor = new JSDraw2.Editor(this._molfilesEditorHost, {viewonly: true});
    }
    this._molfilesEditorRemovingTimer = setTimeout(() => {
      this.removeMolfilesEditor();
    }, 1000);
    return this._molfilesEditor;
  }

  // after batch conversion, remove the editor to avoid sliding pages
  private _molfilesEditorRemovingTimer: any = null;
  private removeMolfilesEditor(): void {
    if (this._molfilesEditorHost) {
      $(this._molfilesEditorHost).empty();
      this._molfilesEditorHost.remove();
      this._molfilesEditorHost = null;
      this._molfilesEditor = null;
    }
  }

  public getMolfiles(helmStrList: string[]): string[] {
    // Phase 5 (hwe migration): pseudo-molfiles (monomers as atoms) now come
    // from the standalone hwe library via the adapter, not the legacy JSDraw2
    // editor.
    return this.editorAdapter.getMolfiles(helmStrList);
  }

  public parse(helm: string, origin?: Point): HelmMol {
    const molHandler = new JSDraw2.MolHandler<HelmType, IHelmBio, IHelmEditorOptions>();
    const plugin = new org.helm.webeditor.Plugin(molHandler);
    org.helm.webeditor.IO.parseHelm(plugin, helm, origin ?? new JSDraw2.Point(0, 0));
    return plugin.jsd.m;
  }

  removeGaps(srcHelm: string): HelmConvertRes {
    // Phase 5 (hwe migration): gap removal + monomer-index remap now runs in the
    // standalone hwe library via the adapter (immutable-mol based), replacing
    // the legacy mutable-HelmMol splice/relink + org.helm.webeditor.IO.getHelm.
    return this.editorAdapter.removeGaps(srcHelm);
  }
}
