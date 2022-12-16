import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';
import {Observable} from 'rxjs';
import {MjolnirPointerEvent} from 'mjolnir.js';
import {PickingInfo} from '@deck.gl/core/typed';


export type NodeStyleType = { [propName: string]: any };
export type StylesType = { [nodeName: string]: NodeStyleType };

export interface INglGlViewer extends IViewer {
  get nwkDf(): DG.DataFrame;

  set nwkDf(value: DG.DataFrame);

  setProps(updater: { [propName: string]: any }): void;

  get onAfterRender(): Observable<{ gl: WebGLRenderingContext }>;

  get onHover(): Observable<{ info: PickingInfo, event: MjolnirPointerEvent }>;
}


export type NglGlTask = {
  name: string,
  backColor: number,
  props: { [propName: string]: any },
  onAfterRender: (canvas: HTMLCanvasElement) => Promise<void>
};

export abstract class NglGlServiceBase {
  public static noneSource: { type: string, data: any } = {
    type: 'biojs',
    data: {name: 'none', branch_length: 1, children: []}
  };

  /** Queues NglGl render task
   * @param key  Specify to skip previously queued tasks with the same key
   */
  abstract render(args: NglGlTask, key?: keyof any): void;

  /** Default implementation of rendering tree on grid cell
   * @param gCtx    Context to draw on grid
   * @param bd      Bound rect to clip drawing on task moment
   * @param gCell   Grid cell to draw
   * @param canvas  Image of the tree rendered
   */
  abstract renderOnGridCell(
    gCtx: CanvasRenderingContext2D, bd: DG.Rect, gCell: DG.GridCell, canvas: CanvasImageSource): void;
}


export async function getNglGlService(): Promise<NglGlServiceBase> {
  const funcList = DG.Func.find({package: 'BiostructureViewer', name: 'getNglGlService'});
  if (funcList.length === 0)
    throw new Error('Package "BiostructureViewer" must be installed for NglGl services.');

  const svc: NglGlServiceBase = (await funcList[0].prepare().call()).getOutputParamValue() as NglGlServiceBase;
  return svc;
}
