import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {printLeftOrCentered} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {findMonomers, parseHelm} from './utils';
import {HelmMonomerPlacer} from './helm-monomer-placer';

const enum tempTAGS {
  helmSumMaxLengthWords = 'helm-sum-maxLengthWords',
  helmMaxLengthWords = 'helm-maxLengthWords',

  helmPlacer = 'bio-helmPlacer',
}

// Global flag is for replaceAll
const helmGapStartRe = /\{(\*\.)+/g;
const helmGapIntRe = /\.(\*\.)+/g;
const helmGapEndRe = /(\.\*)+\}/g;

type TempType = { [tagName: string]: any };

/** Helm cell renderer in case of no missed monomer draws with JSDraw2.Editor (webeditor),
 * in case of missed monomers presented, draws linear sequences aligned in width per monomer.
 */
export class HelmCellRenderer extends DG.GridCellRenderer {
  get name() { return 'helm'; }

  get cellType() { return 'helm'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    try {
      /* Can not do anything without tableColumn containing temp */
      let tableCol: DG.Column | null = null;
      try { tableCol = gridCell.tableColumn; } catch { }
      if (!tableCol) return;

      const helmPlacer = HelmMonomerPlacer.getOrCreate(tableCol);
      const [allParts, lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);

      const maxIndex = Object.values(lengths).length - 1;
      const argsX = e.offsetX - gridCell.bounds.x;
      let left = 0;
      let right = maxIndex;
      let found = false;
      let iterCount: number = 0;

      let mid = 0;
      if (argsX > sumLengths[0]) {
        while (!found && iterCount < sumLengths.length) {
          mid = Math.floor((right + left) / 2);
          if (argsX >= sumLengths[mid] && argsX <= sumLengths[mid + 1]) {
            left = mid;
            found = true;
          } else if (argsX < sumLengths[mid])
            right = mid - 1;
          else if (argsX > sumLengths[mid + 1])
            left = mid + 1;

          if (left == right)
            found = true;

          iterCount++;
        }
      }
      left = (argsX >= sumLengths[left]) ? left : left - 1; // correct left to between sumLengths

      const seq: string = !gridCell.cell.value ? '' : gridCell.cell.value
        .replaceAll(helmGapStartRe, '{').replaceAll(helmGapIntRe, '.').replaceAll(helmGapEndRe, '}')
        .replace('{*}', '{}');
      const monomerList = parseHelm(seq);
      const monomers = new Set<string>(monomerList);
      const missedMonomers = findMonomers(monomerList);

      const tooltipMessage: HTMLElement[] = [];
      for (const [part, partI] of wu.enumerate(allParts)) {
        if (missedMonomers.has(part)) {
          tooltipMessage[partI] = ui.divV([
            ui.divText(`Monomer ${allParts[partI]} not found.`),
            ui.divText('Open the Context Panel, then expand Manage Libraries')
          ]);
        } else if (monomers.has(part)) {
          const elList = [ui.div(part)];
          const monomer = helmPlacer.getMonomer(part);
          if (monomer) {
            const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
            const monomerSvg = grok.chem.svgMol(monomer.smiles, undefined, undefined, options);
            elList.push(monomerSvg);
          }
          tooltipMessage[partI] = ui.divV(elList);
        }
      }

      (((tooltipMessage[left]?.childNodes.length ?? 0) > 0)) ?
        ui.tooltip.show(tooltipMessage[left], e.x + 16, e.y + 16) :
        ui.tooltip.hide();
    } catch (err: any) {
      const errMsg: string = errorToConsole(err);
      console.error('Helm: HelmCellRenderer.onMouseMove() error:\n' + errMsg);
    }
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    g.save();
    try {
      /* Can not do anything without tableColumn containing temp */
      let tableCol: DG.Column | null = null;
      try { tableCol = gridCell.tableColumn; } catch { }
      if (!tableCol) return;

      const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
      const missedColor = 'red';
      const monomerColor: string = '#404040';
      const frameColor: string = '#C0C0C0';

      const seq: string = !gridCell.cell.value ? '' : gridCell.cell.value
        .replaceAll(helmGapStartRe, '{').replaceAll(helmGapIntRe, '.').replaceAll(helmGapEndRe, '}')
        .replace('{*}', '{}');
      const monomerList = parseHelm(seq);
      const monomers: Set<string> = new Set<string>(monomerList);
      const missedMonomers: Set<string> = findMonomers(monomerList);
      const helmPlacer = HelmMonomerPlacer.getOrCreate(tableCol);

      if (missedMonomers.size == 0) {
        helmPlacer.skipCell(gridCell.tableRowIndex!);
        const host = ui.div([], {style: {width: `${w}px`, height: `${h}px`}});
        host.setAttribute('dataformat', 'helm');
        host.setAttribute('data', seq /* gaps skipped */);
        gridCell.element = host;
        //@ts-ignore
        const canvas = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true});
        return;
      }

      if (missedMonomers.size > 0) {
        if (!grid) {
          const r = window.devicePixelRatio;
          h = 28;
          g.canvas.height = h * r;
          g.canvas.style.height = `${h}px`;
        }

        w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
        g.save();
        g.beginPath();
        g.rect(x, y, w, h);
        g.clip();
        g.transform(1, 0, 0, 1, x, y);
        g.font = '12px monospace';
        g.textBaseline = 'top';
        const [allParts, lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);

        for (let i = 0; i < allParts.length; ++i) {
          const part: string = allParts[i];
          const color: string =
            part === '.' || part.endsWith('{') || part.startsWith('}') ? frameColor :
              missedMonomers.has(part) ? missedColor :
                monomers.has(part) ? monomerColor :
                  frameColor;
          g.fillStyle = color;
          printLeftOrCentered(sumLengths[i], 0, w, h, g, allParts[i], color, 0, true, 1.0);
        }
      }
    } finally {
      g.restore();
    }
  }
}
