import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Subject} from 'rxjs';

const COMBINED_MEASUREMENTS_CELL_TYPE = 'combined-measurements';
export const MAX_MARKER_RADIUS = 3;
export const MIN_MARKER_RADIUS = 1.5;
export const PADDING_X = 4;
/** Vertical padding inside cells (exported for overlay drawing). */
export const PADDING_Y = 10;
const HIT_RADIUS = 5;

export function getMarkerRadius(armWidth: number): number {
  return Math.max(MIN_MARKER_RADIUS, Math.min(MAX_MARKER_RADIUS, armWidth * 0.25));
}

export function getEnabledArms(df: DG.DataFrame): string[] {
  const tag = df.getTag('armOrder');
  const globalArms: string[] = tag ? JSON.parse(tag) : [];
  const disabledTag = df.getTag('disabledArms');
  const disabled: string[] = disabledTag ? JSON.parse(disabledTag) : [];
  return globalArms.filter((a) => !disabled.includes(a));
}

export const subjectClicked$ = new Subject<string>();

type ArmEntry = {values: number[], subjects: string[]};
type ArmValues = {[arm: string]: ArmEntry};

function parseArmData(value: any): ArmValues | null {
  if (!value)
    return null;
  try {
    const data = typeof value === 'string' ? JSON.parse(value) : value;
    if (!data || typeof data !== 'object')
      return null;
    return data as ArmValues;
  }
  catch {
    return null;
  }
}

function getGlobalArms(gridCell: DG.GridCell): string[] {
  const tag = gridCell.grid.dataFrame.getTag('armOrder');
  return tag ? JSON.parse(tag) : [];
}

function getSelectedSubject(gridCell: DG.GridCell): string {
  return gridCell.grid.dataFrame.getTag('selectedSubject') ?? '';
}

export function getArmColorStr(armIndex: number, alpha: number = 1): string {
  const color = DG.Color.getCategoricalColor(armIndex);
  const r = (color >> 16) & 0xFF;
  const g = (color >> 8) & 0xFF;
  const b = color & 0xFF;
  return `rgba(${r},${g},${b},${alpha})`;
}

function getRowMinMax(gridCell: DG.GridCell): {min: number, max: number} {
  const df = gridCell.grid.dataFrame;
  const tableRowIdx = gridCell.cell.rowIndex;
  const minCol = df.col('~rowMin');
  const maxCol = df.col('~rowMax');
  if (minCol && maxCol)
    return {min: minCol.get(tableRowIdx), max: maxCol.get(tableRowIdx)};
  return {min: 0, max: 0};
}

type PointInfo = {arm: string, value: number, subject: string, px: number, py: number, armIndex: number};

function computePointPositions(
  armData: ArmValues, globalArms: string[], enabledArms: string[],
  x: number, y: number, w: number, h: number,
  min: number, max: number,
): PointInfo[] {
  const arms = Object.keys(armData);
  if (arms.length === 0 || enabledArms.length === 0)
    return [];

  const range = max - min;
  const armWidth = (w - 2 * PADDING_X) / enabledArms.length;
  const points: PointInfo[] = [];

  for (const arm of arms) {
    const enabledIdx = enabledArms.indexOf(arm);
    if (enabledIdx < 0)
      continue;

    const entry = armData[arm];
    if (!entry || entry.values.length === 0)
      continue;

    const globalIdx = globalArms.indexOf(arm);
    const armCenterX = x + PADDING_X + enabledIdx * armWidth + armWidth / 2;

    for (let i = 0; i < entry.values.length; i++) {
      const v = entry.values[i];
      const subject = entry.subjects[i] ?? '';
      let py: number;
      if (range === 0)
        py = y + h / 2;
      else
        py = y + h - PADDING_Y - ((v - min) / range) * (h - 2 * PADDING_Y);
      points.push({arm, value: v, subject, px: armCenterX, py, armIndex: globalIdx >= 0 ? globalIdx : 0});
    }
  }

  return points;
}

@grok.decorators.cellRenderer({
  name: 'combinedMeasurementsRenderer',
  cellType: 'combined-measurements',
})
export class CombinedMeasurementsCellRenderer extends DG.GridCellRenderer {
  get name(): string {
    return 'combinedMeasurementsRenderer';
  }

  get cellType(): string {
    return COMBINED_MEASUREMENTS_CELL_TYPE;
  }

  get defaultWidth(): number {
    return 120;
  }

  get defaultHeight(): number {
    return 60;
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    const armData = parseArmData(gridCell.cell.value);
    if (!armData)
      return;

    const globalArms = getGlobalArms(gridCell);
    const enabledArms = getEnabledArms(gridCell.grid.dataFrame);
    if (enabledArms.length === 0)
      return;
    const selectedSubject = getSelectedSubject(gridCell);
    const {min, max} = getRowMinMax(gridCell);
    const points = computePointPositions(armData, globalArms, enabledArms, x, y, w, h, min, max);
    const armWidth = (w - 2 * PADDING_X) / enabledArms.length;
    const radius = getMarkerRadius(armWidth);

    g.lineWidth = 1.5;
    for (const pt of points) {
      const isSelected = selectedSubject && pt.subject === selectedSubject;
      if (isSelected) {
        g.strokeStyle = 'rgba(0,0,0,1)';
        g.fillStyle = 'rgba(0,0,0,1)';
        g.beginPath();
        g.arc(pt.px, pt.py, radius, 0, 2 * Math.PI);
        g.fill();
        g.stroke();
      }
      else {
        g.strokeStyle = getArmColorStr(pt.armIndex, 0.3);
        g.beginPath();
        g.arc(pt.px, pt.py, radius, 0, 2 * Math.PI);
        g.stroke();
      }
    }
  }

  private findClosestPoint(gridCell: DG.GridCell, e: MouseEvent): PointInfo | null {
    const armData = parseArmData(gridCell.cell.value);
    if (!armData)
      return null;

    const bounds = gridCell.bounds;
    const globalArms = getGlobalArms(gridCell);
    const enabledArms = getEnabledArms(gridCell.grid.dataFrame);
    const {min, max} = getRowMinMax(gridCell);
    const points = computePointPositions(
      armData, globalArms, enabledArms, bounds.x, bounds.y, bounds.width, bounds.height, min, max);

    const mx = e.offsetX;
    const my = e.offsetY;

    let closestPt: PointInfo | null = null;
    let closestDist = HIT_RADIUS;

    for (const pt of points) {
      const dx = mx - pt.px;
      const dy = my - pt.py;
      const dist = Math.sqrt(dx * dx + dy * dy);
      if (dist < closestDist) {
        closestDist = dist;
        closestPt = pt;
      }
    }

    return closestPt;
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const closestPt = this.findClosestPoint(gridCell, e);
    if (closestPt)
      ui.tooltip.show(ui.divV([
        ui.divText(`Arm: ${closestPt.arm}`),
        ui.divText(`Value: ${closestPt.value}`),
        ui.divText(`Subject: ${closestPt.subject}`),
      ]), e.x + 16, e.y + 16);
    else
      ui.tooltip.hide();
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    const closestPt = this.findClosestPoint(gridCell, e);
    if (closestPt && closestPt.subject)
      subjectClicked$.next(closestPt.subject);
  }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    ui.tooltip.hide();
  }
}
