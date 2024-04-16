import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {
  getSettingsBase,
  names,
  SparklineType,
  SummarySettingsBase,
  createTooltip,
  Hit,
  isSummarySettingsBase
} from './shared';

let minRadius: number;

enum PieChartStyle {
  Radius = 'Radius',
  Angle = 'Angle'
}

interface PieChartSettings extends SummarySettingsBase {
  radius: number;
  style: PieChartStyle.Radius | PieChartStyle.Angle;
  sectors: {
    sectorColor: string;
    subsectors: { name: string; radius: number }[];
  }[];
}

function getSettings(gc: DG.GridColumn): PieChartSettings {
  const sectors = gc.settings.sectors;
  const settings: PieChartSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
    gc.settings[SparklineType.PieChart] ??= getSettingsBase(gc, SparklineType.PieChart);
  settings.style ??= PieChartStyle.Radius;
  settings.sectors ??= sectors;
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

function onHit(gridCell: DG.GridCell, e: MouseEvent): Hit {
  const settings = getSettings(gridCell.gridColumn);
  const cols = gridCell.grid.dataFrame.columns.byNames(settings.columnNames);
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
    r = cols[activeColumn].scale(row) * (gridCell.bounds.width - 4) / 2;
    r = Math.max(r, minRadius);
  } else if (settings.sectors) {
    const sectors = settings.sectors;
    const totalSectors = sectors.length;
    const sectorAngle = (2 * Math.PI) / totalSectors;
    const normalizedAngle = (angle + Math.PI * 2) % (Math.PI * 2);
    for (let i = 0; i < totalSectors; ++i) {
      const sector = sectors[i];
      const sectorStartAngle = i * sectorAngle;
      const sectorEndAngle = (i + 1) * sectorAngle;
      if (normalizedAngle >= sectorStartAngle && normalizedAngle < sectorEndAngle) {
        const subsectors = sector.subsectors;
        const totalSubsectors = subsectors.length;
        const subsectorAngle = sectorAngle / totalSubsectors;
        let subsectorStartAngle = sectorStartAngle;
        for (let j = 0; j < totalSubsectors; ++j) {
          const subsector = subsectors[j];
          const subsectorEndAngle = subsectorStartAngle + subsectorAngle;
          if (normalizedAngle >= subsectorStartAngle && normalizedAngle < subsectorEndAngle) {
            activeColumn = cols.findIndex((col) => col && subsector && col.name === subsector.name);
            break;
          }
          subsectorStartAngle = subsectorEndAngle;
        }
        break;
      }
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

export class PieChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'pie ts'; }

  get cellType() { return SparklineType.PieChart; }

  // getPreferredCellSize(col: DG.GridColumn) {
  //   return new Size(80,80);
  // }

  get defaultWidth(): number | null { return 80; }

  get defaultHeight(): number | null { return 80; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const hitData = onHit(gridCell, e);
    if (hitData.isHit)
      ui.tooltip.show(ui.divV(createTooltip(hitData.cols, hitData.activeColumn, hitData.row)), e.x + 16, e.y + 16);
    else
      ui.tooltip.hide();
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 5 || h < 5 || df === void 0) return;

    const settings = getSettings(gridCell.gridColumn);
    const row: number = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);
    const box = new DG.Rect(x, y, w, h).fitSquare().inflate(-2, -2);
    minRadius = Math.min(box.width, box.height) / 10;
    if (settings.style == PieChartStyle.Radius && !settings.sectors) {
      for (let i = 0; i < cols.length; i++) {
        if (cols[i].isNone(row))
          continue;

        let r = cols[i].scale(row) * box.width / 2;
        r = Math.max(r, minRadius);
        g.beginPath();
        g.moveTo(box.midX, box.midY);
        g.arc(box.midX, box.midY, r,
          2 * Math.PI * i / cols.length, 2 * Math.PI * (i + 1) / cols.length);
        g.closePath();

        g.fillStyle = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
        g.fill();
        g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
        g.stroke();
      }
    } else if (settings.sectors) {
      let currentAngle = 0;
      const sectors = settings.sectors;
      for (let i = 0; i < sectors.length; ++i) {
        const sector = sectors[i];
        const sectorAngle = (2 * Math.PI) / sectors.length;
        g.beginPath();
        g.moveTo(box.midX, box.midY);
        g.arc(box.midX, box.midY, Math.min(box.width, box.height) / 2, currentAngle, currentAngle + sectorAngle);
        g.closePath();
        g.fillStyle = this.hexToRgbA(sector.sectorColor, 0.4);
        g.fill();
        
        const totalSubsectors = sector.subsectors.length;
        let subsectorCurrentAngle = currentAngle;
        sector.subsectors.forEach(subsector => {
          const subsectorAngle = sectorAngle / totalSubsectors;
          const gap = 0.05;
          let r = Math.max(Math.min(box.width, box.height) / 2, minRadius);
          const subsectorName = subsector.name;
          const subsectorCol = cols.find(col => col.name === subsectorName);
          if (subsectorCol) {
            r = subsectorCol.scale(row) * (Math.min(box.width, box.height) / 2);
            r = Math.max(r, minRadius);
          }
          g.beginPath();
          g.moveTo(box.midX, box.midY);
          g.arc(box.midX, box.midY, r, subsectorCurrentAngle + gap, subsectorCurrentAngle + subsectorAngle - gap);
          g.closePath();
          g.fillStyle = this.hexToRgbA(sector.sectorColor, 0.8);
          g.fill();
          subsectorCurrentAngle += subsectorAngle;
        });
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

        g.fillStyle = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
        g.fill();
        g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
        g.stroke();
        currentAngle = endAngle;
      }
    }
  }

  hexToRgbA(hex: string, opacity: number): string {
    const bigint = parseInt(hex.substring(1), 16);
    const r = (bigint >> 16) & 255;
    const g = (bigint >> 8) & 255;
    const b = bigint & 255;
    return `rgba(${r},${g},${b},${opacity})`;
  }
  
  renderSettings(gc: DG.GridColumn): Element {
    const settings: PieChartSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
      gc.settings[SparklineType.PieChart] ??= getSettings(gc);

    return ui.inputs([
      ui.columnsInput('Сolumns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
      ui.choiceInput('Style', PieChartStyle.Radius, [PieChartStyle.Angle, PieChartStyle.Radius],
        function(value: PieChartStyle) {
          settings.style = value;
          gc.grid.invalidate();
        }),
    ]);
  }
}
