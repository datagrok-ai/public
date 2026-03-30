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
import {HelmTabKeys, HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HelmConvertRes, HelmInputBase, IHelmHelper, IHelmInputInitOptions
} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLibBase, IMonomerLinkData, IMonomerSetPlaceholder} from '@datagrok-libraries/bio/src/types/monomer-library';
import {GAP_SYMBOL, GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {HelmWebEditor} from './helm-web-editor';
import {OrgHelmModule, ScilModule} from './types';
import {HelmInput} from './widgets/helm-input';
import {getHoveredMonomerFromEditorMol} from './utils/get-hovered';
/* eslint-enable max-len */

import {_package} from './package';
import {IEditorOptions} from '@datagrok-libraries/js-draw-lite/src/types/jsdraw2';

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

  createHelmWebEditor(host?: HTMLDivElement, options?: Partial<IHelmEditorOptions>): IHelmWebEditor {
    return new HelmWebEditor(host, options);
  }

  createWebEditorApp(host: HTMLDivElement, helm: string): App {
    org.helm.webeditor.MolViewer.molscale = 0.8;
    const webEditorApp: App = new org.helm.webeditor.App(host, {
      showabout: false,
      mexfontsize: '90%',
      mexrnapinontab: true,
      topmargin: 20,
      mexmonomerstab: true,
      sequenceviewonly: false,
      mexfavoritefirst: true,
      mexfilter: true,
      currentTabKey: HelmTabKeys.Helm,
      overrideTabs: (tabs: Partial<TabDescType>[]): Partial<TabDescType>[] => {
        const res: Partial<TabDescType>[] = [...tabs,
          {caption: 'Placeholders', tabkey: 'placeholders'},
        ];
        return res;
      },
      onShowTab: (mex: MonomerExplorer, div: HTMLDivElement, key: string): void => {
        switch (key) {
        case 'placeholders': {
          this.onShowTabPlaceholders(mex, div);
          break;
        }
        case 'placeholders-meta': {
          this.onShowTabPlaceholdersMeta(mex, div);
          break;
        }
        case 'placeholders-sets': {
          this.onShowTabPlaceholdersSets(mex, div);
          break;
        }
        }
        dojo.connect(div, 'onmousedown', function(e: MouseEvent) {
          mex.select(e);
        });
        dojo.connect(div, 'ondblclick', function(e: MouseEvent) {
          mex.dblclick(e);
        });
      }
    });
    const sizes = webEditorApp.calculateSizes();
    webEditorApp.canvas!.resize(sizes.rightwidth - 100, sizes.topheight - 210);
    let s: Partial<CSSStyleDeclaration> = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
    scil.apply(webEditorApp.sequence.style, s);
    scil.apply(webEditorApp.notation!.style, s);
    s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + webEditorApp.toolbarheight) + 'px'};
    scil.apply(webEditorApp.properties!.parent.style, s);
    webEditorApp.structureview!.resize(sizes.rightwidth, sizes.bottomheight + webEditorApp.toolbarheight);
    webEditorApp.mex!.resize(sizes.topheight - 80);
    webEditorApp.canvas?.setHelm(helm);
    return webEditorApp;
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
    try {
      const editor = this.molfilesEditor;
      const resList = helmStrList.map((helmStr, i) => {
        //no-draw instruction to avoid time spent on rendering
        editor.helm && (editor.helm.noDraw = true);
        editor.setHelm(helmStr);
        const mol = editor.getMolfile();
        return mol;
      });
      return resList;
    } finally {
      // $(host).empty();
      // host.remove();
    }
  }

  public parse(helm: string, origin?: Point): HelmMol {
    const molHandler = new JSDraw2.MolHandler<HelmType, IHelmBio, IHelmEditorOptions>();
    const plugin = new org.helm.webeditor.Plugin(molHandler);
    org.helm.webeditor.IO.parseHelm(plugin, helm, origin ?? new JSDraw2.Point(0, 0));
    return plugin.jsd.m;
  }

  removeGaps(srcHelm: string): HelmConvertRes {
    const monomerMap = new Map<number, number>();
    const mol: HelmMol = this.parse(srcHelm);
    for (let aI: number = mol.atoms.length - 1; aI >= 0; --aI) {
      const a = mol.atoms[aI];
      if (a.elem === HELM_GAP_SYMBOL /* '*' - original */) {
        // Collect all bonds connected to this gap atom, tracking which side the gap is on
        const connectedBonds: { bI: number, gapSide: 1 | 2 }[] = [];
        for (let bI: number = mol.bonds.length - 1; bI >= 0; --bI) {
          const b = mol.bonds[bI];
          if (b.a1 === a) connectedBonds.push({bI, gapSide: 1}); // gap is a1
          else if (b.a2 === a) connectedBonds.push({bI, gapSide: 2}); // gap is a2
        }

        // Separate into backbone bonds (R1/R2) and other bonds (R3+ cross-polymer)
        // A backbone bond connects to the gap's R1 or R2 attachment point
        let backboneR1Bond: { bI: number, gapSide: 1 | 2 } | null = null; // bond using gap's R1
        let backboneR2Bond: { bI: number, gapSide: 1 | 2 } | null = null; // bond using gap's R2
        const otherBonds: number[] = []; // bond indices to just remove

        for (const cb of connectedBonds) {
          const bond = mol.bonds[cb.bI];
          // The R-group on the gap's side tells us which attachment point is used
          const gapR = cb.gapSide === 1 ? bond.r1 : bond.r2;
          const gapRNum = typeof gapR === 'string' ? parseInt(gapR.replace(/\D/g, '')) : gapR;
          if (gapRNum === 1)
            backboneR1Bond = cb;
          else if (gapRNum === 2)
            backboneR2Bond = cb;
          else
            otherBonds.push(cb.bI);
        }

        // Remove other (non-backbone) bonds first (iterate in reverse to keep indices valid)
        otherBonds.sort((a, b) => b - a);
        for (const bI of otherBonds)
          mol.bonds.splice(bI, 1);

        // Recalculate bond indices after splice (otherBonds were removed)
        // We need to adjust backboneR1Bond/R2Bond indices
        const adjustIdx = (origBi: number): number => {
          let adj = origBi;
          for (const removedBi of otherBonds)
            if (removedBi < origBi) adj--;
          return adj;
        };

        const r1 = backboneR1Bond ? {...backboneR1Bond, bI: adjustIdx(backboneR1Bond.bI)} : null;
        const r2 = backboneR2Bond ? {...backboneR2Bond, bI: adjustIdx(backboneR2Bond.bI)} : null;

        if (r1 && r2) {
          // Re-link: the neighbor on the R1 side connects directly to the neighbor on the R2 side
          const r1Bond = mol.bonds[r1.bI];
          const r2Bond = mol.bonds[r2.bI];
          // The "other" atom (not the gap) on each bond
          const r1Neighbor = r1.gapSide === 1 ? r1Bond.a2 : r1Bond.a1;
          const r2Neighbor = r2.gapSide === 1 ? r2Bond.a2 : r2Bond.a1;
          const r1NeighborR = r1.gapSide === 1 ? r1Bond.r2 : r1Bond.r1;
          const r1NeighborApo = r1.gapSide === 1 ? r1Bond.apo2 : r1Bond.apo1;
          const r2NeighborR = r2.gapSide === 1 ? r2Bond.r2 : r2Bond.r1;
          const r2NeighborApo = r2.gapSide === 1 ? r2Bond.apo2 : r2Bond.apo1;

          // Rewrite the R1 bond to connect r1Neighbor ↔ r2Neighbor directly
          r1Bond.a1 = r1Neighbor;
          r1Bond.r1 = r1NeighborR;
          r1Bond.apo1 = r1NeighborApo;
          r1Bond.a2 = r2Neighbor;
          r1Bond.r2 = r2NeighborR;
          r1Bond.apo2 = r2NeighborApo;

          // Remove the R2 bond (no longer needed)
          mol.bonds.splice(r2.bI, 1);
        } else if (r1 && !r2) {
          // Gap at the end (only R1 bond) — just remove
          mol.bonds.splice(r1.bI, 1);
        } else if (!r1 && r2) {
          // Gap at the start (only R2 bond) — just remove
          mol.bonds.splice(r2.bI, 1);
        }
        // else: isolated gap, no backbone bonds — just remove atom

        mol.atoms.splice(aI, 1);
      }
    }

    for (let aI: number = 0; aI < mol.atoms.length; ++aI) {
      const a: HelmAtom = mol.atoms[aI];
      monomerMap.set(a.bio!.continuousId - 1, aI);
    }

    const resHelm = org.helm.webeditor.IO.getHelm(mol)!;
    return {srcHelm, resHelm, monomerMap};
  }
}
