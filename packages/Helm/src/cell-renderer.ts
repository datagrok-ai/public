import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import {findMonomers, parseHelm, getParts} from './utils';
import {printLeftOrCentered} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

const enum tempTAGS {
  helmSumMaxLengthWords = 'helm-sum-maxLengthWords',
  helmMaxLengthWords = 'helm-maxLengthWords',
}

type TempType = { [tagName: string]: any };

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

      const colTemp: TempType[] = tableCol.temp ?? new Array<TempType>(tableCol.length);
      // Exit if no missed monomers (tags are not presented in colTemp)
      if (!colTemp || Object.keys(colTemp).length == 0) return;

      const maxLengthWordsSum: { [pos: number]: number } = colTemp[tempTAGS.helmSumMaxLengthWords];
      const maxLengthWords: { [pos: number]: number } = colTemp[tempTAGS.helmMaxLengthWords];

      const maxIndex = Object.values(maxLengthWords).length - 1;
      const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCell.bounds.x);
      let left = 0;
      let right = maxIndex;
      let found = false;
      maxLengthWordsSum[maxIndex + 1] = argsX + 1;
      let mid = 0;
      if (argsX > maxLengthWordsSum[0]) {
        while (!found) {
          mid = Math.floor((right + left) / 2);
          if (argsX >= maxLengthWordsSum[mid] && argsX <= maxLengthWordsSum[mid + 1]) {
            left = mid;
            found = true;
          } else if (argsX < maxLengthWordsSum[mid]) {
            right = mid - 1;
          } else if (argsX > maxLengthWordsSum[mid + 1]) {
            left = mid + 1;
          }
          if (left == right) {
            found = true;
          }
        }
      }
      left = (argsX >= maxLengthWordsSum[left]) ? left + 1 : left;
      const monomers = findMonomers(gridCell.cell.value);
      const s: string = gridCell.cell.value ?? '';
      const subParts: string[] = parseHelm(s);
      const allParts: string[] = getParts(subParts, s);
      const tooltipMessage: HTMLElement[] = [];
      for (let partI = 0; partI < allParts.length; ++partI) {
        if (monomers.has(allParts[partI]))
          tooltipMessage[partI] = ui.divV([
            ui.divText(`Monomer ${allParts[partI]} not found.`),
            ui.divText('Open the Context Panel, then expand Manage Libraries')
          ]);
      }
      (((tooltipMessage[left]?.childNodes.length ?? 0) > 0))
        ? ui.tooltip.show(ui.div(tooltipMessage[left]), e.x + 16, e.y + 16)
        : ui.tooltip.hide();
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
      const undefinedColor = 'rgb(100,100,100)';
      const grayColor = '#808080';

      const missedMonomers = findMonomers(gridCell.cell.value);
      let s: string = gridCell.cell.value ?? '';
      let subParts: string[] = parseHelm(s);

      if (missedMonomers.size == 0) {
        const host = ui.div([], {style: {width: `${w}px`, height: `${h}px`}});
        host.setAttribute('dataformat', 'helm');
        host.setAttribute('data', gridCell.cell.value);
        gridCell.element = host;
        //@ts-ignore
        const canvas = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true});
        return;
      }

      if (missedMonomers.size > 0) {
        const maxLengthWords: number[] = tableCol.temp[tempTAGS.helmMaxLengthWords] ?? [];
        if (subParts.length > maxLengthWords.length)
          maxLengthWords.push(...(new Array<number>(subParts.length - maxLengthWords.length).fill(-1)));

        w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
        g.save();
        g.beginPath();
        g.rect(x, y, w, h);
        g.clip();
        g.font = '12px monospace';
        g.textBaseline = 'top';
        let x1 = x;
        const allParts: string[] = getParts(subParts, s);
        for (let i = 0; i < allParts.length; ++i) {
          maxLengthWords[i] = Math.max(maxLengthWords[i], allParts[i].length * 7); /* What is 7, width of char ? */
          const color = missedMonomers.has(allParts[i]) ? 'red' : grayColor;
          g.fillStyle = undefinedColor;
          x1 = printLeftOrCentered(x1, y, w, h, g, allParts[i], color, 0, true, 1.0);
        }

        const maxLengthWordSum: number[] = new Array<number>(maxLengthWords.length);
        maxLengthWordSum[0] = maxLengthWords[0];
        for (let partI = 1; partI < allParts.length; partI++) {
          maxLengthWordSum[partI] = maxLengthWordSum[partI - 1] + maxLengthWords[partI];
        }
        tableCol.temp = {
          [tempTAGS.helmSumMaxLengthWords]: maxLengthWordSum,
          [tempTAGS.helmMaxLengthWords]: maxLengthWords
        };
        return;
      }
    } finally {
      g.restore();
    }
  }
}
