import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import type {App} from '@datagrok-libraries/helm-web-editor/helm/App';
import {GetMonomerFunc, MonomersFuncs, HelmMol, HelmString, IHelmWebEditor} from './types';
import {IMonomerLib} from '../types/index';

export type IHelmInputInitOptions = ui.input.IInputInitOptions<HelmString | HelmMol>;

export abstract class HelmInputBase extends DG.JsInputBase<HelmString> {
  abstract get molValue(): HelmMol;
  abstract set molValue(value: HelmMol);
}

export interface IHelmHelper {
  createHelmInput(name: string, options?: IHelmInputInitOptions): HelmInputBase;

  createHelmWebEditor(host?: HTMLElement): IHelmWebEditor;

  createWebEditorApp(host: HTMLDivElement, helm?: string): App;

  get originalMonomersFuncs(): MonomersFuncs | null;

  buildMonomersFuncsFromLib(monomerLib: IMonomerLib): MonomersFuncs;

  overrideMonomersFuncs(monomersFuncs: MonomersFuncs): MonomersFuncs;
  revertOriginalMonomersFuncs(): MonomersFuncs;
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
