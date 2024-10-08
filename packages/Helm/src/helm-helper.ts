import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

/* eslint-disable max-len */
import {
  App, Point, HelmAtom, HelmBond, HelmMol, HelmType, GetMonomerFunc, GetMonomerResType,
  MonomerExplorer, TabDescType, MonomersFuncs, MonomerSetType, ISeqMonomer, PolymerType, MonomerType,
  DojoType, JSDraw2ModuleType,
} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HelmConvertRes, HelmInputBase, HelmNotSupportedError, IHelmHelper, IHelmInputInitOptions
} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLibBase, IMonomerLinkData, IMonomerSetPlaceholder} from '@datagrok-libraries/bio/src/types/index';
import {HelmTabKeys, IHelmDrawOptions} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';
import {GAP_SYMBOL, GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';

import {HelmWebEditor} from './helm-web-editor';
import {OrgHelmModule, ScilModule} from './types';
import {HelmInput} from './widgets/helm-input';
import {getHoveredMonomerFromEditorMol} from './utils/get-hovered';
/* eslint-enable max-len */

import {_package} from './package';

declare const dojo: DojoType;
declare const JSDraw2: JSDraw2ModuleType;
declare const scil: ScilModule;
declare const org: OrgHelmModule;

const HELM_GAP_SYMBOL: string = GapOriginals[NOTATION.HELM];

export class HelmHelper implements IHelmHelper {
  private static instanceCount: number = 0;

  constructor(
    private readonly logger: ILogger
  ) {
    // Watchdog for singleton
    if ((++HelmHelper.instanceCount) > 1)
      throw new Error(`HelmHelper must be a single.`);
  }

  createHelmInput(name?: string, options?: IHelmInputInitOptions): HelmInputBase {
    try {
      const monomerLib: IMonomerLibBase = _package.libHelper!.getMonomerLib();
      return HelmInput.create(this, monomerLib, name, options);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(`Helm: HelmHelper.createHelmInput(), Error: ${errMsg}`, undefined, errStack);
      throw err;
    }
  }

  createHelmWebEditor(host?: HTMLDivElement, drawOptions?: Partial<IHelmDrawOptions>): IHelmWebEditor {
    return new HelmWebEditor(host, drawOptions);
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
    const monomerSets = _package.libHelper!.getMonomerSets();

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

  public getMolfiles(helmStrList: string[]): string[] {
    const host = ui.div([], {style: {width: '0', height: '0'}});
    document.documentElement.appendChild(host);
    try {
      const editor = new JSDraw2.Editor(host, {viewonly: true});
      const resList = helmStrList.map((helmStr, i) => {
        editor.setHelm(helmStr);
        const mol = editor.getMolfile();
        return mol;
      });
      return resList;
    } finally {
      $(host).empty();
      host.remove();
    }
  }

  public parse(helm: string, origin?: Point): HelmMol {
    const molHandler = new JSDraw2.MolHandler<HelmType>();
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
        const leftBondList: { aI: number, bI: number }[] = [];
        const rightBondList: { aI: number, bI: number }[] = [];
        for (let bI: number = mol.bonds.length - 1; bI >= 0; --bI) {
          const b = mol.bonds[bI];
          if (b.a1 !== a && b.a2 === a) leftBondList.push({aI, bI});
          if (b.a1 === a && b.a2 !== a) rightBondList.push({aI, bI});
        }
        if (leftBondList.length > 1 || rightBondList.length > 1)
          throw new HelmNotSupportedError(`Removing a gap monomer #${aI} with more than two bonds is unsupported.`);

        const lb = leftBondList[0];
        const rb = rightBondList[0];
        if (lb && !rb) {
          mol.bonds.splice(lb.bI, 1);
          mol.atoms.splice(lb.aI, 1);
        } else if (!lb && rb) {
          const rb = rightBondList[0];
          mol.bonds.splice(rb.bI, 1);
          mol.atoms.splice(rb.aI, 1);
        } else {
          if (lb.aI !== rb.aI)
            throw new Error('Something is really wrong here.');

          // is not enough, breaks the simple polymer
          //mol.delAtom(a, true);

          const leftBond = mol.bonds[lb.bI];
          const rightBond = mol.bonds[rb.bI];
          const a2 = leftBond.a2 = rightBond.a2; // right atom
          leftBond.r2 = rightBond.r2;
          leftBond.apo2 = rightBond.apo2;

          if (a2.bonds)
            throw new Error('Bond list of the atom is not corrected.');
          // a2.bonds!.splice(a2.bonds!.indexOf(rightBond), 1); // remove the right bond from links of the right atom
          // a2.bonds!.push(leftBond); // put the left bond to links of the right atom

          mol.bonds.splice(rb.bI, 1); // remove right bond
          mol.atoms.splice(lb.aI, 1);
        }
      }
    }

    for (let aI: number = 0; aI < mol.atoms.length; ++aI) {
      const a: HelmAtom = mol.atoms[aI];
      monomerMap.set(parseInt(a.bio!.continuousId as string) - 1, aI);
    }

    const resHelm = org.helm.webeditor.IO.getHelm(mol)!;
    return {srcHelm, resHelm, monomerMap};
  }
}
