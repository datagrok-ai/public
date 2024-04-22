import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';

import {HelmProps, HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';

import {IEditor} from '../helm-monomer-placer';

import {_package} from '../package';

export class HelmService extends HelmServiceBase {
  private readonly hostDiv: HTMLDivElement;

  private editor: IEditor | null = null;

  constructor() {
    super(_package.logger);

    this.hostDiv = ui.box();
    // this.hostDiv.style.display = 'none'; // Disables drawing at all
    this.hostDiv.style.position = 'absolute';
    this.hostDiv.style.left = '0px';
    this.hostDiv.style.right = '0px';
    this.hostDiv.style.width = '0px';
    this.hostDiv.style.height = '0px';
    this.hostDiv.style.visibility = 'hidden';
    document.body.appendChild(this.hostDiv);
  }

  protected override async requestRender(
    key: keyof any | undefined, task: RenderTask<HelmProps>, renderHandler: () => void
  ): Promise<[Unsubscribable, number, () => void]> {
    const logPrefix = `${this.toLog()}.requestRender()`;
    this.logger.debug(`${logPrefix}, ` + `key: ${key?.toString()}`);
    const emptyCanvasHash: number = 0;

    if (!this.editor) {
      this.editor = new JSDraw2.Editor(this.hostDiv, {width: 300, height: 300, skin: 'w8', viewonly: true}) as IEditor;
    }

    const gridCell = task.props.gridCell;

    // let tableCol: DG.Column | null = null;
    // try { tableCol = gridCell.tableColumn; } catch { }
    // if (!tableCol) return [undefined as unknown as SignalBinding, -1, () => {}];
    //
    // const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
    // const missedColor = 'red';
    // const monomerColor: string = '#404040';
    // const frameColor: string = '#C0C0C0';
    //
    // const seq: string = !gridCell.cell.value ? '' : removeGapsFromHelm(gridCell.cell.value);
    // const monomerList = parseHelm(seq);
    // const monomers: Set<string> = new Set<string>(monomerList);
    // const missedMonomers: Set<string> = findMonomers(monomerList);
    // const helmPlacer = HelmMonomerPlacer.getOrCreate(tableCol);
    //
    // let w: number = task.props.width;
    // let h: number = task.props.height;
    //
    // if (missedMonomers.size == 0) {
    //   // Recreate host to avoid hanging in window.dojox.gfx.svg.Text.prototype.getTextWidth
    //   const host = gridCell.element = ui.div([],
    //     {style: {width: `${w - 2}px`, height: `${h - 2}px`, margin: `${1}px`, backgroundColor: '#FFE0E0'}});
    //   host.setAttribute('dataformat', 'helm');
    //   host.setAttribute('data', seq /* gaps skipped */);
    //   // if grid has neighbour to the left, then shift host to the left
    //   if (host.parentElement && (gridCell.grid?.canvas?.offsetLeft ?? 0) > 0) {
    //     host.parentElement.style.left =
    //       `${(gridCell.grid?.canvas?.offsetLeft ?? 0) + host.parentElement.offsetLeft}px`;
    //   }
    //
    //   // Recreate editor to avoid hanging in window.dojox.gfx.svg.Text.prototype.getTextWidth
    //   const editor = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true}) as IEditor;
    //   helmPlacer.setEditor(gridCell.tableRowIndex!, editor);
    //
    //   helmPlacer.skipCell(gridCell.tableRowIndex!);
    //   return [undefined as unknown as SignalBinding, ];
    // }
    //
    // if (missedMonomers.size > 0) {
    //   if (!grid) {
    //     const r = window.devicePixelRatio;
    //     h = 28;
    //     g.canvas.height = h * r;
    //     g.canvas.style.height = `${h}px`;
    //   }
    //
    //   w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
    //   //g.save();
    //   g.beginPath();
    //   g.rect(x, y, w, h);
    //   g.clip();
    //   g.transform(1, 0, 0, 1, x, y);
    //   g.font = '12px monospace';
    //   g.textBaseline = 'top';
    //   const [allParts, _lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);
    //
    //   for (let i = 0; i < allParts.length; ++i) {
    //     const part: string = allParts[i];
    //     const color: string =
    //       part === '.' || part.endsWith('{') || part.startsWith('}') ? frameColor :
    //         missedMonomers.has(part) ? missedColor :
    //           monomers.has(part) ? monomerColor :
    //             frameColor;
    //     g.fillStyle = color;
    //     printLeftOrCentered(sumLengths[i], 0, w, h, g, allParts[i], color, 0, true, 1.0,
    //       undefined, undefined, undefined, undefined, undefined,
    //       undefined, undefined, undefined, helmPlacer.monomerTextSizeMap);
    //   }
    // }

    const trigger = () => {};
    return [undefined as unknown as Unsubscribable, emptyCanvasHash, trigger];
  }

  protected onRendered(
    key: keyof any | undefined, task: RenderTask<HelmProps>, emptyCanvasHash: number
  ): boolean {
    return true;
  }

  async reset(): Promise<void> { }
}
