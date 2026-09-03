import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {desirabilityScore, isNumerical, PropertyDesirability} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {
  ColumnGroup,
  createBaseInputs,
  createTooltip, getRenderColor,
  getSettingsBase,
  Hit,
  isSummarySettingsBase,
  SparklineType,
  SummaryColumnColoringType,
  SummarySettingsBase,
  NormalizationType, getScaledNumber, scaleSettings, getSparklinesContextPanel
} from './shared';
import {VlaaiVisEditor} from '../vlaaivis/editor';
import {LABELS, STYLE_INFO} from '../vlaaivis/constants';

let minRadius: number;

enum PieChartStyle {
  Radius = 'Radius',
  Angle = 'Angle',
  Vlaaivis = 'VlaaiVis'
}

// PropertyDesirability is a union (numerical | categorical), and an interface cannot extend
// one — an intersection distributes over it and keeps every member's shape.
export type Subsector = PropertyDesirability & {
  name: string;
};

export interface Sector {
  name: string;
  sectorColor: string;
  subsectors: Subsector[];
}

export interface PieChartSettings extends SummarySettingsBase {
  radius: number;
  style: PieChartStyle.Radius | PieChartStyle.Angle | PieChartStyle.Vlaaivis;
  sectors?: {
      name?: string;
      lowerBound: number;
      upperBound: number;
      sectors: Sector[]; // Use the Sector interface here
      values: string | null;
  };
}

function getSettings(gc: DG.GridColumn): PieChartSettings {
  const sectors = gc.settings.sectors;
  const settings: PieChartSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
    gc.settings[SparklineType.PieChart] ??= getSettingsBase(gc, SparklineType.PieChart);
  settings.style ??= PieChartStyle.Radius;
  settings.sectors ??= sectors;
  settings.colorCode ??= SummaryColumnColoringType.Bins;
  settings.normalization ??= NormalizationType.Column;
  return settings;
}

function getColumnsSum(cols: DG.Column[], row: number) {
  let sum = 0;
  for (let i = 0; i < cols.length; i++) {
    if (cols[i].isNone(row))
      continue;
    sum += cols[i].getNumber(row);
  }
  return sum;
}

function normalizeValue(subsector: Subsector, col: DG.Column, row: number): number | null {
  if (col.isNone(row)) {
    const missing = subsector.missingValues;
    return missing?.strategy === 'default' ? missing.score : null;
  }
  if (!isNumerical(subsector))
    return subsector.categories.find((c) => c.name === col.get(row))?.desirability ?? null;
  return desirabilityScore(subsector, col.getNumber(row));
}

function renderSubsector(
  g: CanvasRenderingContext2D, box: DG.Rect, sectorColor: string,
  sectorAngle: number, currentAngle: number, subsector: Subsector,
  minRadius: number, cols: DG.Column[], row: number,
  sectorWeight: number
): number {
  const normalizedSubsectorWeight = subsector.weight / sectorWeight;
  const subsectorAngle = sectorAngle * normalizedSubsectorWeight;
  const radiusFactor = Math.min(box.width, box.height) / 2;
  const subsectorCol = cols.find((col) => col.name === subsector.name);
  const score = subsectorCol ? normalizeValue(subsector, subsectorCol, row) : null;
  const r = Math.max((score ?? 0) * radiusFactor, minRadius);

  g.beginPath();
  g.moveTo(box.midX, box.midY);
  g.arc(box.midX, box.midY, r, currentAngle, currentAngle + subsectorAngle);
  g.closePath();
  g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
  g.lineWidth = 0.6;
  g.stroke();
  if (score !== null) {
    g.fillStyle = hexToRgbA(sectorColor, 0.6);
    g.fill();
  }
  return currentAngle + subsectorAngle;
}

function hexToRgbA(hex: string, opacity: number): string {
  const bigint = parseInt(hex.substring(1), 16);
  const r = (bigint >> 16) & 255;
  const g = (bigint >> 8) & 255;
  const b = bigint & 255;
  return `rgba(${r},${g},${b},${opacity})`;
}

function calculateSectorWeight(sector: { sectorColor: string; subsectors: Subsector[]; }): number {
  return sector.subsectors.reduce((acc, subsector) => acc + subsector.weight, 0);
}

function getColumns(gridCell: DG.GridCell, settings: PieChartSettings): DG.Column[] {
  return gridCell.grid.dataFrame.columns.byNames(settings.columnNames).filter((c) => c != null);
}

function onHit(gridCell: DG.GridCell, e: MouseEvent, settings: PieChartSettings): Hit {
  const cols = getColumns(gridCell, settings);
  const vectorX = e.offsetX - gridCell.bounds.midX;
  const vectorY = e.offsetY - gridCell.bounds.midY;
  const distance = Math.sqrt(vectorX * vectorX + vectorY * vectorY);
  const atan2 = Math.atan2(vectorY, vectorX);
  const angle = atan2 < 0 ? atan2 + 2 * Math.PI : atan2;
  let activeColumn = -1;
  const row: number = gridCell.cell.row.idx;

  let r: number = (gridCell.bounds.width - 4) / 2;
  if (settings.style == PieChartStyle.Radius && !settings.sectors) {
    activeColumn = Math.floor((angle * cols.length) / (2 * Math.PI));
    if (cols[activeColumn] !== null) {
      const scaledNumber = getScaledNumber(cols, row, cols[activeColumn],
        scaleSettings(settings, cols[activeColumn])
      );
      r = scaledNumber * (gridCell.bounds.width - 4) / 2;
      r = Math.max(r, minRadius);
    }
  } else if (settings.sectors) {
    const {sectors} = settings.sectors;
    let currentAngle = 0;
    const totalSectorWeight = sectors.reduce((acc, sector) => acc + calculateSectorWeight(sector), 0);
    for (const sector of sectors) {
      const sectorWeight = calculateSectorWeight(sector);
      const normalizedSectorWeight = sectorWeight / totalSectorWeight;
      const sectorAngle = 2 * Math.PI * normalizedSectorWeight;
      const sectorStartAngle = currentAngle;
      const sectorEndAngle = currentAngle + sectorAngle;

      if (angle >= sectorStartAngle && angle < sectorEndAngle) {
        const subsectors = sector.subsectors;
        const totalSubsectorWeight = subsectors.reduce((acc, subsector) => acc + subsector.weight, 0);
        let subsectorStartAngle = sectorStartAngle;
        for (const subsector of subsectors) {
          const subsectorWeight = subsector.weight;
          const normalizedSubsectorWeight = subsectorWeight / totalSubsectorWeight;
          const subsectorAngle = sectorAngle * normalizedSubsectorWeight;
          const subsectorEndAngle = subsectorStartAngle + subsectorAngle;
          if (angle >= subsectorStartAngle && angle < subsectorEndAngle) {
            activeColumn = cols.findIndex((col) => col && col.name === subsector.name);
            break;
          }
          subsectorStartAngle = subsectorEndAngle;
        }
        break;
      }
      currentAngle += sectorAngle;
    }
  } else {
    const sum = getColumnsSum(cols, row);
    r = (gridCell.bounds.width - 4) / 2;

    let currentAngle = 0;
    for (let i = 0; i < cols.length; i++) {
      if (cols[i].isNone(gridCell.cell.row.idx))
        continue;
      const endAngle = currentAngle + 2 * Math.PI * cols[i].getNumber(row) / sum;
      if ((angle > currentAngle) && (angle < endAngle)) {
        activeColumn = i;
        break;
      }
      currentAngle = endAngle;
    }
  }

  return {
    isHit: (r >= distance),
    activeColumn: activeColumn,
    row: row,
    cols: cols,
  };
}

function columnGroups(settings: PieChartSettings, cols: DG.Column[]): ColumnGroup[] {
  if (!settings.sectors)
    return [{name: '', cols}];
  const sectors = settings.sectors.sectors;
  const assigned = new Set(sectors.flatMap((s) => s.subsectors.map((p) => p.name)));
  return [
    ...sectors.map((s) => ({name: s.name, color: s.sectorColor,
      cols: s.subsectors.flatMap((p) => cols.filter((c) => c.name === p.name))})),
    {name: LABELS.UNASSIGNED, cols: cols.filter((c) => !assigned.has(c.name))},
  ].filter((g) => g.cols.length > 0);
}

export class PieChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'pie ts'; }

  get cellType() { return SparklineType.PieChart; }

  // getPreferredCellSize(col: DG.GridColumn) {
  //   return new Size(80,80);
  // }

  get defaultWidth(): number | null { return 80; }

  get defaultHeight(): number | null { return 80; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const settings = getSettings(gridCell.gridColumn);
    const hitData = onHit(gridCell, e, settings);
    if (!hitData.isHit) {
      ui.tooltip.hide();
      return;
    }
    const groups = columnGroups(settings, hitData.cols);
    ui.tooltip.show(createTooltip(hitData.cols, hitData.activeColumn, hitData.row, groups), e.x + 16, e.y + 16);
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 5 || h < 5 || df === void 0) return;

    const settings = getSettings(gridCell.gridColumn);
    let row: number = gridCell.cell.row.idx;
    let cols = getColumns(gridCell, settings);
    const box = new DG.Rect(x, y, w, h).fitSquare().inflate(-2, -2);
    minRadius = Math.min(box.width, box.height) / 10;
    if (settings.style == PieChartStyle.Radius && !settings.sectors) {
      for (let i = 0; i < cols.length; i++) {
        if (cols[i] === null || row === -1 || cols[i].isNone(row))
          continue;

        const scaledNumber = getScaledNumber(cols, row, cols[i],
          scaleSettings(settings, cols[i])
        );
        let r = scaledNumber * box.width / 2;
        r = Math.max(r, minRadius);
        g.beginPath();
        g.moveTo(box.midX, box.midY);
        g.arc(box.midX, box.midY, r,
          2 * Math.PI * i / cols.length, 2 * Math.PI * (i + 1) / cols.length);
        g.closePath();

        g.fillStyle = DG.Color.toRgb(getRenderColor(settings, DG.Color.blue, {column: cols[i], colIdx: i, rowIdx: row}));
        g.fill();
        g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
        g.stroke();
      }
    } else if (settings.sectors) {
      const {lowerBound, upperBound, sectors, values} = settings.sectors;
      cols = values ? Array.from(DG.DataFrame.fromCsv(values).columns) : cols;
      row = values ? 0 : row;
      let currentAngle = 0;
      const totalSectorWeight = sectors.reduce((acc, sector) => acc + calculateSectorWeight(sector), 0);

      for (const sector of sectors) {
        const sectorWeight = calculateSectorWeight(sector);
        const normalizedSectorWeight = sectorWeight / totalSectorWeight;
        const sectorAngle = 2 * Math.PI * normalizedSectorWeight;
        const radiusFactor = Math.min(box.width, box.height) / 2;
        const arcEnd = currentAngle + sectorAngle;

        // Render inner circle representing the range
        g.beginPath();
        g.arc(box.midX, box.midY, lowerBound * radiusFactor, currentAngle, arcEnd);
        g.arc(box.midX, box.midY, upperBound * radiusFactor, arcEnd, currentAngle, true);
        g.fillStyle = hexToRgbA(sector.sectorColor, 0.2);
        g.fill();

        // Render subsectors
        let subsectorCurrentAngle = currentAngle;
        for (const subsector of sector.subsectors)
          subsectorCurrentAngle = renderSubsector(g, box, sector.sectorColor, sectorAngle, subsectorCurrentAngle, subsector, minRadius, cols, row, sectorWeight);

        currentAngle += sectorAngle;
      }
    } else {
      const sum = getColumnsSum(cols, row);
      let currentAngle = 0;
      for (let i = 0; i < cols.length; i++) {
        if (cols[i].isNone(row))
          continue;
        const r = box.width / 2;
        const endAngle = currentAngle + 2 * Math.PI * cols[i].getNumber(row) / sum;
        g.beginPath();
        g.moveTo(box.midX, box.midY);
        g.arc(box.midX, box.midY, r, currentAngle, endAngle);
        g.closePath();

        g.fillStyle = DG.Color.toRgb(getRenderColor(settings, DG.Color.blue, {column: cols[i], colIdx: i, rowIdx: row}));
        g.fill();
        g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
        g.stroke();
        currentAngle = endAngle;
      }
    }
  }

  renderSettings(gc: DG.GridColumn): Element {
    const settings: PieChartSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
      gc.settings[SparklineType.PieChart] ??= getSettings(gc);

    const elementsDiv = ui.div([]);
    const styleInfo = ui.icons.info(() => {});
    styleInfo.hidden = true;
    ui.tooltip.bind(styleInfo,
      () => ui.divV([ui.divText(STYLE_INFO.SUMMARY), ui.link(STYLE_INFO.LEARN_MORE, STYLE_INFO.HELP_URL)]));
    let editor: VlaaiVisEditor | null = null;
    let stashedSectors: PieChartSettings['sectors'];

    const showEditor = (style: PieChartStyle) => {
      editor?.detach();
      editor = null;
      ui.empty(elementsDiv);
      const isVlaaivis = style === PieChartStyle.Vlaaivis;
      styleInfo.hidden = !isVlaaivis;
      for (const input of scalingInputs)
        input.visible = !isVlaaivis;
      if (!isVlaaivis) {
        if (settings.sectors)
          stashedSectors = settings.sectors;
        delete settings.sectors;
        gc.grid.invalidate();
        return;
      }
      settings.sectors ??= stashedSectors;
      editor = new VlaaiVisEditor(settings, gc);
      inputs.append(editor.profileInput.root, ...editor.boundsInputs.map((input) => input.root));
      elementsDiv.appendChild(editor.root);
      gc.grid.invalidate();
    };

    const [columnsInput, ...scalingInputs] = createBaseInputs(gc, settings);
    columnsInput.onChanged.subscribe(() => {
      if (editor)
        editor.refresh();
    });

    const style = settings.style ?? PieChartStyle.Radius;
    const styleInput = ui.input.choice('Style', {value: style,
      items: [PieChartStyle.Angle, PieChartStyle.Radius, PieChartStyle.Vlaaivis],
      onValueChanged: (value) => {
        settings.style = value;
        showEditor(value);
      }});
    styleInput.addOptions(styleInfo);
    const inputs = ui.inputs([
      styleInput,
      columnsInput,
      ...scalingInputs,
    ]);
    if (style === PieChartStyle.Vlaaivis)
      showEditor(style);

    return ui.divV([inputs, elementsDiv], 'power-grid-pie-settings');
  }

  hasContextValue(gridCell: DG.GridCell): boolean { return true; }
  async getContextValue(gridCell: DG.GridCell): Promise<any> {
    const settings = getSettings(gridCell.gridColumn);
    const groups = columnGroups(settings, getColumns(gridCell, settings));
    return getSparklinesContextPanel(gridCell, settings.columnNames, groups);
  }
}
