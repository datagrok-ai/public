import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import type {App} from '@datagrok-libraries/helm-web-editor/helm/App';
import {HelmMol, IHelmWebEditor} from './types';
import {IMonomerLib} from '../types/index';
import {GetMonomerFunc} from '../helm/types';

// TODO: Remove duplicate while the one ui is exported
export interface IInputInitOptions<T = any> {
  value?: T;
  property?: DG.Property;
  nullable?: boolean;
  elementOptions?: DG.ElementOptions;
  tooltipText?: string;
  onCreated?: (input: DG.InputBase<T>) => void;
  onValueChanged?: (input: DG.InputBase<T>) => void;
}

export interface IHelmHelper {
  createHelmInput(name: string, options?: IInputInitOptions<HelmMol>): DG.InputBase<HelmMol>;

  createHelmWebEditor(host?: HTMLElement): IHelmWebEditor;

  createWebEditorApp(host: HTMLDivElement, helm?: string): App;

  get originalGetMonomer(): GetMonomerFunc | null;

  buildGetMonomerFromLib(monomerLib: IMonomerLib): GetMonomerFunc;
  overrideGetMonomer(getMonomerFunc: GetMonomerFunc): GetMonomerFunc;
  revertOriginalGetMonomer(): GetMonomerFunc;
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
    export function helmAsync(name: string, options?: IInputInitOptions<HelmMol>): Promise<DG.InputBase<HelmMol>>;
  }
}

ui.input.helmAsync = async function(
  name: string, options?: IInputInitOptions<HelmMol>): Promise<DG.InputBase<HelmMol>> {
  return (await getHelmHelper()).createHelmInput(name, options);
};
