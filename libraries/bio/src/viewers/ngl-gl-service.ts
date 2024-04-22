import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PropsBase, RenderServiceBase} from '../utils/cell-renderer-async-base';

export class NglGlProps extends PropsBase {
  public constructor(
    public readonly pdb: string,
    backColor: number, width: number, height: number,
  ) {
    super(backColor, width, height);
  }
}

export abstract class NglGlServiceBase extends RenderServiceBase<NglGlProps> {}


export async function getNglGlService(): Promise<NglGlServiceBase> {
  const funcList = DG.Func.find({package: 'BiostructureViewer', name: 'getNglGlService'});
  if (funcList.length === 0)
    throw new Error('Package "BiostructureViewer" must be installed for NglGl services.');

  const svc: NglGlServiceBase = (await funcList[0].prepare().call()).getOutputParamValue() as NglGlServiceBase;
  return svc;
}
