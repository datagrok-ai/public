import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PropsBase, RenderServiceBase} from '../utils/cell-renderer-async-base';

export class HelmProps extends PropsBase {
  public constructor(
    public readonly gridCell: DG.GridCell,
    backColor: number, width: number, height: number
  ) {
    super(backColor, width, height);
  }
}

export abstract class HelmServiceBase extends RenderServiceBase<HelmProps> {}

export async function getHelmService(): Promise<HelmServiceBase> {
  const funcList = DG.Func.find({package: 'Helm', name: 'getHelmService'});
  if (funcList.length === 0)
    throw new Error('Package "Helm" must be installed for Helm services.');

  const svc: HelmServiceBase = (await funcList[0].prepare().call()).getOutputParamValue() as HelmServiceBase;
  return svc;
}
