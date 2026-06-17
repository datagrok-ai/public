/* eslint-disable max-len */
import {
  App, Point, HelmAtom, HelmMol, MonomersFuncs,
  IHelmEditorOptions,
} from '@datagrok-libraries/bio/src/helm/types';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HelmConvertRes, HelmInputBase, IHelmHelper, IHelmInputInitOptions
} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {HelmInput} from './widgets/helm-input';
import {getHoveredMonomerFromEditorMol} from './utils/get-hovered';
/* eslint-enable max-len */

import {_package} from './package';

// hwe migration (Phases 3-7): every editor / parse / calc path now runs through
// the standalone `@datagrok-libraries/hwe` library via its Datagrok adapter
// (`HelmHelperAdapter` over `HelmService`). The legacy forked Pistoia HELM Web
// Editor (`org.helm.webeditor` / JSDraw2 / Dojo) is gone — this class no longer
// touches any of those globals.
import {HelmHelperAdapter, bridgeMonomerLib, HelmService} from '@datagrok-libraries/hwe';
import type {IMonomerLibBaseLike, MonomersFuncsLike} from '@datagrok-libraries/hwe';

export class HelmHelper implements IHelmHelper {
  private static instanceCount: number = 0;

  constructor(
    public readonly seqHelper: ISeqHelper,
    private readonly logger: ILogger,
    // Datagrok's RDKit module (from `getRdKitModule()`), forwarded to hwe so the
    // editor can render monomer structures (atom / palette / RNA-builder tooltips
    // and the bond-props dialog). hwe consumes only `get_mol(...).get_svg(...)`.
    private readonly rdKitModule?: RDModule,
  ) {
    // Watchdog for singleton
    if ((++HelmHelper.instanceCount) > 1)
      throw new Error(`HelmHelper must be a single.`);
  }

  // -- hwe editor adapter --
  // The adapter owns one `HelmService` reading Bio's monomer library through
  // `bridgeMonomerLib`, so the editor / parser / calculators resolve exactly
  // the same monomers as the rest of the platform. Built lazily on first use
  // (after `completeInit` has populated `_package._libHelper`).
  private _editorAdapter: HelmHelperAdapter | null = null;
  private get editorAdapter(): HelmHelperAdapter {
    if (this._editorAdapter === null) {
      const monomerLib = _package._libHelper!.getMonomerLib() as unknown as IMonomerLibBaseLike;
      const service = new HelmService({monomerLib: bridgeMonomerLib(monomerLib), rdkitModule: this.rdKitModule});
      this._editorAdapter = new HelmHelperAdapter({service, seqHelper: this.seqHelper});
    }
    return this._editorAdapter;
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
    return this.editorAdapter.createHelmWebEditor(host, {
      ...(options?.width !== undefined ? {width: options.width} : {}),
      ...(options?.height !== undefined ? {height: options.height} : {}),
    }) as unknown as IHelmWebEditor;
  }

  /**
   * Returns the standalone hwe editor (palette + toolbar + canvas +
   * Sequence/HELM/Properties/Atomic tabs) wrapped in a `LegacyAppWrapper` that
   * presents the legacy Pistoia `App` shape — `canvas.getHelm(true)`,
   * `canvas.helm.jsd.m.atoms[i].selected`, `sequence` / `notation` /
   * `properties.parent`, etc. — so the existing dialog consumers swap in
   * unchanged.
   */
  createWebEditorApp(host: HTMLDivElement, helm?: string): App {
    return this.editorAdapter.createWebEditorApp(host, helm) as unknown as App;
  }

  // -- GetMonomer (monomer-funcs hooks) --
  // Delegated to the adapter. The legacy implementation patched the global
  // `org.helm.webeditor.Monomers` registry so the forked editor used Bio's
  // library; the hwe editor reads `bridgeMonomerLib(...)` directly, so these
  // are now thin pass-throughs kept for the `IHelmHelper` contract.

  public get originalMonomersFuncs(): MonomersFuncs | null {
    return this.editorAdapter.originalMonomersFuncs as unknown as MonomersFuncs | null;
  }

  public revertOriginalMonomersFuncs(): MonomersFuncs {
    return this.editorAdapter.revertOriginalMonomersFuncs() as unknown as MonomersFuncs;
  }

  public overrideMonomersFuncs(monomersFuncs: MonomersFuncs): MonomersFuncs {
    return this.editorAdapter.overrideMonomersFuncs(
      monomersFuncs as unknown as MonomersFuncsLike) as unknown as MonomersFuncs;
  }

  public buildMonomersFuncsFromLib(monomerLib: IMonomerLib): MonomersFuncs {
    return this.editorAdapter.buildMonomersFuncsFromLib(
      monomerLib as unknown as IMonomerLibBaseLike) as unknown as MonomersFuncs;
  }

  public getHoveredAtom(x: number, y: number, mol: HelmMol, height: number): HelmAtom | null {
    return getHoveredMonomerFromEditorMol(x, y, mol, height);
  }

  public getMolfiles(helmStrList: string[]): string[] {
    return this.editorAdapter.getMolfiles(helmStrList);
  }

  public parse(helm: string, origin?: Point): HelmMol {
    return this.editorAdapter.parse(helm, origin) as unknown as HelmMol;
  }

  removeGaps(srcHelm: string): HelmConvertRes {
    return this.editorAdapter.removeGaps(srcHelm);
  }
}
