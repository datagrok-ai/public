import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {
  CellRendererBackAsyncBase, RenderServiceBase
} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {HelmProps} from '@datagrok-libraries/bio/src/viewers/helm-service';

import {HelmMonomerPlacer, IEditor, ISeqMonomer} from '../helm-monomer-placer';
import {findMonomers, parseHelm, removeGapsFromHelm} from './index';
import {getHoveredMonomerFallback, getHoveredMonomerFromEditor} from './get-hovered';
import {_getHelmService} from '../package-utils';

import {_package} from '../package';

export const enum Temps {
  renderer = '.renderer.helm'
}

class HelmGridCellRendererBack extends CellRendererBackAsyncBase<HelmProps> {
  constructor(gridCol: DG.GridColumn) {
    super(gridCol, _package.logger);
  }

  protected override getRenderService(): RenderServiceBase<HelmProps> {
    return _getHelmService();
  }

  protected override getRenderTaskProps(gridCell: DG.GridCell, dpr: number): HelmProps {
    return new HelmProps(
      gridCell,
      gridCell.grid.props.backColor,
      gridCell.gridColumn.width * dpr,
      gridCell.grid.props.rowHeight * dpr);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    /* Can not do anything without tableColumn containing temp */
    let tableCol: DG.Column | null = null;
    try { tableCol = gridCell.tableColumn; } catch { }
    if (!tableCol) return;

    /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
    const argsX = e.offsetX - gcb.x;
    const argsY = e.offsetY - gcb.y;

    //@ts-ignore
    const helmPlacer: HelmMonomerPlacer = HelmMonomerPlacer.getOrCreate(tableCol);
    const editor: IEditor | null = helmPlacer.getEditor(gridCell.tableRowIndex!);
    let seqMonomer: ISeqMonomer | null;
    let missedMonomers: Set<string> = new Set<string>(); // of .size = 0
    if (editor)
      seqMonomer = getHoveredMonomerFromEditor(argsX, argsY, gridCell, editor);
    else {
      const seq: string = !gridCell.cell.value ? '' : removeGapsFromHelm(gridCell.cell.value as string);
      const monomerList = parseHelm(seq);
      missedMonomers = findMonomers(monomerList);
      const parsedMonomers = new Set<string>(monomerList);
      seqMonomer = getHoveredMonomerFallback(argsX, argsY, gridCell, helmPlacer);
      if (seqMonomer && !parsedMonomers.has(seqMonomer.symbol)) seqMonomer = null;
      if (seqMonomer) {
        const textSize = helmPlacer.monomerTextSizeMap[seqMonomer.symbol];
        if (textSize) {
          const textBaseLine = gcb.height / 2 -
            (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2 + 1;
          const textTop = textBaseLine - textSize.fontBoundingBoxAscent;
          const textBottom = textBaseLine + textSize.fontBoundingBoxDescent;
          if (argsY < textTop || textBottom < argsY) seqMonomer = null;
        }
      }
    }

    if (seqMonomer) {
      if (!missedMonomers.has(seqMonomer.symbol)) {
        const monomer = helmPlacer.getMonomer(seqMonomer);
        if (monomer) {
          const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
          const monomerSvg = grok.chem.svgMol(monomer.smiles, undefined, undefined, options);
          ui.tooltip.show(ui.divV([
            ui.divText(seqMonomer.symbol),
            monomerSvg,
          ]), e.x + 16, e.y + 16);
        }
      } else {
        ui.tooltip.show(ui.divV([
          ui.divText(`Monomer '${seqMonomer.symbol}' not found.`),
          ui.divText('Open the Context Panel, then expand Manage Libraries'),
        ]), e.x + 16, e.y + 16);
      }
    } else if (missedMonomers.size == 0) {
      ui.tooltip.hide();
      return;
    } else { // seqMonomer == null && missedMonomers.size > 0
      const mmStrList = wu(missedMonomers.keys()).toArray().sort()
        .filter((_, i) => i < 3);
      const missedMonomersStr = mmStrList.join(', ') + (missedMonomers.size > 3 ? ', ...' : '');
      ui.tooltip.show(ui.divV([ui.divText('Monomers missed in monomer libraries:'),
        ui.divText(missedMonomersStr)
      ]), e.x + 16, e.y + 16);
    }
  }

  static getOrCreate(gridCol: DG.GridColumn): HelmGridCellRendererBack {
    let res: HelmGridCellRendererBack = gridCol.temp[Temps.renderer];
    if (!res) res = gridCol.temp[Temps.renderer] = new HelmGridCellRendererBack(gridCol);
    return res;
  }
}

export class HelmGridCellRenderer extends DG.GridCellRenderer {
  get name() { return 'helm'; }

  get cellType() { return 'helm'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
    const logPrefix = `Helm: HelmGridCellRenderer.onMouseMove()`;
    try {
      HelmGridCellRendererBack.getOrCreate(gridCell.gridColumn).onMouseMove(gridCell, e);
    } catch (err: any) {
      const errMsg: string = errorToConsole(err);
      console.error(`${logPrefix} error:\n` + errMsg);
    } finally {
      e.preventDefault();
      e.stopPropagation();
    }
  }

  override render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ): void {
    HelmGridCellRendererBack.getOrCreate(gridCell.gridColumn).render(g, x, y, w, h, gridCell, cellStyle);
  }
}
