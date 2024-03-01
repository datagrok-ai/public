import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {IHelmWebEditor, IWebEditorApp} from './types';


export interface IHelmHelper {
  createHelmWebEditor(): IHelmWebEditor;

  createWebEditorApp(host: HTMLDivElement, helm?: string): IWebEditorApp;
}

export async function getHelmHelper(): Promise<IHelmHelper> {
  const packageName = 'Helm';
  const funcList = DG.Func.find({package: packageName, name: `getHelmHelper`});
  if (funcList.length === 0)
    throw new Error(`Package '${packageName}' must be installed for HelmHelper.`);
  const res = (await funcList[0].prepare().call()).getOutputParamValue() as IHelmHelper;
  return res;
}
