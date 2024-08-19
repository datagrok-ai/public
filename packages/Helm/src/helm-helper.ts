import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  App, Atom, HelmType, GetMonomerFunc, GetMonomerResType, HelmString, HelmMol,
  MonomerExplorer, TabDescType, MonomersFuncs, MonomerSetType, HelmAtom, ISeqMonomer, PolymerType, MonomerType, type DojoType
} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes, MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {HelmInputBase, IHelmHelper, IHelmInputInitOptions} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLinkData, IMonomerSetPlaceholder} from '@datagrok-libraries/bio/src/types/index';
import {HelmTabKeys, IHelmDrawOptions} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';

import {HelmWebEditor} from './helm-web-editor';
import {OrgHelmModule, ScilModule} from './types';
import {getWebEditorMonomer} from './utils/get-monomer';
import {HelmInput} from './widgets/helm-input';
import {getHoveredMonomerFromEditorMol} from './utils/get-hovered';

import {_package} from './package';

declare const dojo: DojoType;
declare const scil: ScilModule;
declare const org: OrgHelmModule;

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
      const monomerLib: IMonomerLib = _package.libHelper!.getMonomerLib();
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
      // currentTabKey: HelmTabKeys.Helm,
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
          const dgWem = getWebEditorMonomer(monomerLib, a, name);

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
}
