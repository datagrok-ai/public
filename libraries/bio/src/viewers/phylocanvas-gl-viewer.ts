import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';
import {Observable} from 'rxjs';
import {MjolnirPointerEvent} from 'mjolnir.js';
import {PickingInfo} from '@deck.gl/core/typed';


export type NodeStyleType = { [propName: string]: any };
export type StylesType = { [nodeName: string]: NodeStyleType };

export enum TreeTypesNames {
  Radial = 'Radial',
  /** Rectangular edges, leaves listed __vertically__ */
  Rectangular = 'Rectangular',
  Polar = 'Polar',
  Diagonal = 'Diagonal',
  /** Rectangular edges, leaves listed __horizontally__ */
  Orthogonal = 'Orthogonal',
}

export interface IPhylocanvasGlViewer extends IViewer {
  get nwkDf(): DG.DataFrame;

  set nwkDf(value: DG.DataFrame);

  setProps(updater: { [propName: string]: any }): void;

  get onAfterRender(): Observable<{ gl: WebGLRenderingContext }>;

  get onHover(): Observable<{ info: PickingInfo, event: MjolnirPointerEvent }>;
}

// export interface IPhylocanvasGlRenderer {
//   get onAfterRender(): Observable<HTMLCanvasElement>;
// }

export type PhylocanvasGlTask = {
  name: string,
  backColor: number, props: { [propName: string]: any },
  onAfterRender: CanvasCallback
};

export type CanvasCallback = (canvas: HTMLCanvasElement) => void;

export abstract class PhylocanvasGlServiceBase {
  public static noneSource: { type: string, data: any } = {
    type: 'biojs',
    data: {name: 'none', branch_length: 1, children: []}
  };

  /** Queues PhylocanvasGL render task
   * @param key  Specify to skip previously queued tasks with the same key
   */
  abstract render(args: PhylocanvasGlTask, key?: keyof any): void;

  /** Default implementation of rendering tree on grid cell
   * @param gCtx    Context to draw on grid
   * @param bd      Bound rect to clip drawing on task moment
   * @param gCell   Grid cell to draw
   * @param canvas  Image of the tree rendered
   */
  abstract renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource): void;
}


export async function getPhylocanvasGlService(): Promise<PhylocanvasGlServiceBase> {
  const funcList = DG.Func.find({package: 'PhyloTreeViewer', name: 'getPhylocanvasGlService'});
  if (funcList.length === 0)
    throw new Error('Package "PhyloTreeViewer"" must be installed for PhylocanvasGL services.');

  const svc: PhylocanvasGlServiceBase = (await funcList[0].prepare().call())
    .getOutputParamValue() as PhylocanvasGlServiceBase;
  return svc;
}
