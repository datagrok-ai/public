import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type NglGlTask = {
  name: string,
  backColor: number,
  props: { width: number, height: number } & { [propName: string]: any },
  onAfterRender: (canvas: HTMLCanvasElement) => void
};

export abstract class NglGlServiceBase {
  abstract get errorCount(): number;

  /** Disposes MolstarViewer to free WebGL context */
  abstract reset(): Promise<void>;

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
