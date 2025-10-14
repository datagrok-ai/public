import * as DG from 'datagrok-api/dg';
import {MMPA} from './mmp-analysis/mmpa';
import {MMP_NAMES} from './mmp-viewer/mmp-constants';

export class MmmpaActivityRenderer extends DG.GridCellRenderer {
  mmpa: MMPA;
  activityIndex: number;
  color: string;

  constructor(mmpa: MMPA, activityIndex: number, color: string) {
    super();
    this.mmpa = mmpa;
    this.activityIndex = activityIndex;
    this.color = color;
  }

  get name(): string {return 'MMPA activity cell renderer';}
  get cellType(): string {return 'MmpaActivity';}
  get defaultWidth() {return 150;}
  get defaultHeight() {return 150;}

  render(g: any, x: number, y: number, w: number, h: number, cell: DG.GridCell): void {
    const ruleIdx = cell.tableRowIndex!;
    const value = cell.cell.value;
    const ctx = (g.canvas as HTMLCanvasElement).getContext('2d');

    if (!ctx) return;

    const rule = this.mmpa.rules.rules[ruleIdx];
    const pairs = rule.pairs;

    if (!pairs || pairs.length === 0) return;

    // Get all activity values
    const allValues: number[] = [];
    for (const pair of pairs) {
      const fromValue = this.mmpa.allCasesBased.pairedActivities[this.activityIndex][0][pair.id!];
      const toValue = this.mmpa.allCasesBased.pairedActivities[this.activityIndex][1][pair.id!];
      if (fromValue !== DG.FLOAT_NULL && toValue !== DG.FLOAT_NULL)
        allValues.push(fromValue, toValue);
    }

    if (allValues.length === 0) return;

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);

    const drawValue = () => {
      if (value !== DG.FLOAT_NULL) {
        ctx.fillStyle = '#000';
        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        ctx.fillText(value.toFixed(2), x + w/2, y + h/2);
      }
    };

    // If canvas width is too small (20px or less), just show the value
    if (w <= 35) {
      drawValue();
      return;
    }

    // Calculate boundaries with padding
    const range = max - min;

    // Position constants
    const leftX = x + 5; // Right-shifted left anchor
    const rightX = x + w - 5; // Left-shifted right anchor
    const topY = y + 10;
    const bottomY = y + h - 10;
    const visHeight = bottomY - topY;

    const valueToY = (value: number) => {
      const normalized = 1 - ((value - min) / range);
      return topY + (normalized * visHeight);
    };

    // Draw connection lines
    const lineOpacity = 0.5;
    const colorWithOpacity = `${this.color}${Math.floor(lineOpacity * 255).toString(16).padStart(2, '0')}`;
    ctx.strokeStyle = colorWithOpacity;
    ctx.lineWidth = 1;

    for (const pair of pairs) {
      const fromValue = this.mmpa.allCasesBased.pairedActivities[this.activityIndex][0][pair.id!];
      const toValue = this.mmpa.allCasesBased.pairedActivities[this.activityIndex][1][pair.id!];

      if (fromValue === DG.FLOAT_NULL || toValue === DG.FLOAT_NULL) continue;

      const fromY = valueToY(fromValue);
      const toY = valueToY(toValue);

      ctx.beginPath();
      ctx.moveTo(leftX, fromY);
      ctx.lineTo(rightX, toY);
      ctx.stroke();
    }

    // Draw mean value
    drawValue();
  }
}
