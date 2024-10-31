import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {HelmMol, HelmType, IHelmBio, Mol} from '../helm/types';
import {PropsBase, RenderServiceBase} from '../utils/cell-renderer-async-base';
import {IMonomerLibBase} from '../types/index';

export class HelmProps extends PropsBase {
  public constructor(
    public readonly helm: string,
    public readonly monomerLib: IMonomerLibBase,
    backColor: number, width: number, height: number
  ) {
    super(backColor, width, height);
  }
}

export type HelmAux = {
  /** The molecule made of atoms (monomers) and bonds */ mol: HelmMol,
  /** SVG bounding box */ bBox: DG.Rect,
  /** Draw bounding box */ dBox: DG.Rect,
  /** Cell box {0, 0, w, h} */ cBox: DG.Rect,
}

export abstract class HelmServiceBase extends RenderServiceBase<HelmProps, HelmAux> {}

export async function getHelmService(): Promise<HelmServiceBase> {
  const funcList = DG.Func.find({package: 'Helm', name: 'getHelmService'});
  if (funcList.length === 0)
    throw new Error('Package "Helm" must be installed for Helm services.');

  const svc: HelmServiceBase = (await funcList[0].prepare().call()).getOutputParamValue() as HelmServiceBase;
  return svc;
}
