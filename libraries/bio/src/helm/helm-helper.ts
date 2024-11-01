import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable} from 'rxjs';

import type {Point, App, HelmType, IHelmEditorOptions} from './types';
import {
  GetMonomerFunc, MonomersFuncs, HelmMol, HelmString, IHelmWebEditor, HelmAtom, IHelmDrawOptions
} from './types';
import {IMonomerLibBase} from '../types/index';
import {ISeqHelper} from '../utils/seq-helper';
import {SeqValueBase} from '../utils/macromolecule/seq-handler';

export type IHelmInputInitOptions = ui.input.IInputInitOptions<SeqValueBase> & {
  editorOptions: Partial<IHelmEditorOptions>;
  editable: boolean;
};

export abstract class HelmInputBase extends DG.JsInputBase<SeqValueBase> {
  abstract get molValue(): HelmMol;
  abstract set molValue(value: HelmMol);

  abstract get onMouseMove(): Observable<MouseEvent>;

  abstract get onClick(): Observable<MouseEvent>;

  abstract redraw(): void;

  abstract showTooltip(content: HTMLElement | string, a: HelmAtom): void;
}

/**
 * @property {Map<number, number>} monomerMap srcPosIdx -> resPosIdx
 */
export type HelmConvertRes = {
  srcHelm: string;
  resHelm: string;
  monomerMap: Map<number, number> | null;
}

export const HelmNotSupportedErrorType = 'HelmNotSupportedError';

export class HelmNotSupportedError extends Error {
  public readonly type = HelmNotSupportedErrorType;

  constructor(message?: string) {
    super(message);
  }
}

export interface IHelmHelper {
  get seqHelper(): ISeqHelper;

  createHelmInput(name: string, options?: IHelmInputInitOptions): HelmInputBase;

  createHelmWebEditor(host?: HTMLElement, options?: Partial<IHelmEditorOptions>): IHelmWebEditor;

  createWebEditorApp(host: HTMLDivElement, helm?: string): App;

  get originalMonomersFuncs(): MonomersFuncs | null;

  buildMonomersFuncsFromLib(monomerLib: IMonomerLibBase): MonomersFuncs;

  overrideMonomersFuncs(monomersFuncs: MonomersFuncs): MonomersFuncs;
  revertOriginalMonomersFuncs(): MonomersFuncs;

  getHoveredAtom(x: number, y: number, mol: HelmMol, height: number): HelmAtom | null;

  /** Gets pseudo molfiles with monomers as atoms */
  getMolfiles(helmStrList: string[]): string[];

  parse(helm: string, origin?: Point): HelmMol;
  removeGaps(helm: string): HelmConvertRes;
}

export async function getHelmHelper(): Promise<IHelmHelper> {
  const packageName = 'Helm';
  const funcList = DG.Func.find({package: packageName, name: `getHelmHelper`});
  if (funcList.length === 0)
    throw new Error(`Package '${packageName}' must be installed for HelmHelper.`);
  const res = (await funcList[0].prepare().call()).getOutputParamValue() as IHelmHelper;
  return res;
}

declare module 'datagrok-api/ui' {
  export namespace input {
    /** To create HelmInput synchronously, get the {@link IHelmHelper} object
     * via {@link getHelmHelper} in advance, and call {@link IHelmHelper.createHelmInput}.
     */
    export function helmAsync(name: string, options?: IHelmInputInitOptions): Promise<HelmInputBase>;
  }
}

ui.input.helmAsync = async function(
  name: string, options?: IHelmInputInitOptions): Promise<HelmInputBase> {
  return (await getHelmHelper()).createHelmInput(name, options);
};

/** Inputs logic */
export function getMonomerHandleArgs(
  a: HelmAtom | HelmType, name?: string
): [/** biotype */ HelmType, /** elem */ string] {
  if (!a)
    throw new Error(`Argument 'a' of type Atom or HelmType is mandatory.`);
  let biotype: HelmType;
  let elem: string;
  const aa = a as HelmAtom;
  if (aa.T === 'ATOM') {
    biotype = aa.biotype()!;
    elem = aa.elem;
  } else {
    biotype = a as HelmType;
    elem = name!;
  }
  return [biotype, elem];
}
