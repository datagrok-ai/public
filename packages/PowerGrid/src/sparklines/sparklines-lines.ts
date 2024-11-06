import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {
  getSettingsBase,
  SparklineType,
  SummarySettingsBase,
  Hit,
  distance,
  createTooltip,
  isSummarySettingsBase, SummaryColumnColoringType, createBaseInputs, getRenderColor
} from './shared';

const minDistance = 5;

// interface for getPos function data
interface getPosConstants {
  b: DG.Rect;
  settings: SparklineSettings;
  cols: DG.Column[];
}


function getPos(col: number, row: number, constants: getPosConstants): DG.Point {
  const b = constants.b;
  const settings = constants.settings;
  const cols = constants.cols;
  const gmin = settings.globalScale ? Math.min(...cols.map((c: DG.Column) => c.min)) : 0;
  const gmax = settings.globalScale ? Math.max(...cols.map((c: DG.Column) => c.max)) : 0;
  const r: number = settings.globalScale ? (cols[col].getNumber(row) - gmin) / (gmax - gmin) : cols[col].scale(row);
  return new DG.Point(
    b.left + b.width * (cols.length == 1 ? 0 : col / (cols.length - 1)),
    (b.top + b.height) - b.height * r);
}

interface SparklineSettings extends SummarySettingsBase {
  globalScale: boolean;
}

function getSettings(gc: DG.GridColumn): SparklineSettings {
  const settings: SparklineSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
    gc.settings[SparklineType.Sparkline] ??= getSettingsBase(gc, SparklineType.Sparkline);
  settings.globalScale ??= false;
  settings.colorCode ??= SummaryColumnColoringType.Bins;
  return settings;
}

function onHit(gridCell: DG.GridCell, e: MouseEvent): Hit {
  const df = gridCell.grid.dataFrame;

  if (gridCell.bounds.width < 20 || gridCell.bounds.height < 10 || df === void 0) {
    return {
      isHit: false,
      row: -1,
      cols: [],
      activeColumn: -1
    };
  }

  const row = gridCell.cell.row.idx;
  const settings = getSettings(gridCell.gridColumn);
  const b = new DG.Rect(gridCell.bounds.x, gridCell.bounds.y, gridCell.bounds.width,
    gridCell.bounds.height).inflate(-3, -2);

  const cols = df.columns.byNames(settings.columnNames);
  const getPosConstants: getPosConstants = {
    b: b,
    settings: settings,
    cols: cols
  };

  const MousePoint = new DG.Point(e.offsetX, e.offsetY);
  const activeColumn = Math.floor((MousePoint.x - b.left + Math.sqrt(minDistance)) /
    b.width * (cols.length - 1 > 0 ? cols.length - 1 : 1));
  if (activeColumn >= cols.length || activeColumn === -1) {
    return {
      isHit: false,
      row: -1,
      cols: [],
      activeColumn: -1
    };
  }

  const activePoint = getPos(activeColumn, row, getPosConstants);
  return {
    isHit: distance(activePoint, MousePoint) < minDistance,
    activeColumn: activeColumn,
    row: row,
    cols: cols,
  };
}

export class SparklineCellRenderer extends DG.GridCellRenderer {
  get name() { return SparklineType.Sparkline; }

  get cellType() { return SparklineType.Sparkline; }

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

    if (w < 20 || h < 10 || df === void 0) return;

    const settings = getSettings(gridCell.gridColumn);
    const b = new DG.Rect(x, y, w, h).inflate(-3, -2);
    g.strokeStyle = 'lightgrey';
    g.lineWidth = 1;

    const row = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);

    const getPosConstants: getPosConstants = {
      b: b,
      settings: settings,
      cols: cols
    };

    g.beginPath();
    let started = false;
    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        const p = getPos(i, row, getPosConstants);

        if (!started) {
          g.moveTo(p.x, p.y);
          started = true;
        } else {
          g.lineTo(p.x, p.y);
        }
      }
    }
    g.stroke();

    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        const p = getPos(i, row, getPosConstants);
        DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, p.x, p.y, getRenderColor(settings, DG.Color.blue,
          {column: cols[i], colIdx: i, rowIdx: row}), 3);
      }
    }
  }

  renderSettings(gridColumn: DG.GridColumn): HTMLElement {
    const settings: SparklineSettings = isSummarySettingsBase(gridColumn.settings) ? gridColumn.settings :
      gridColumn.settings[SparklineType.Sparkline] ??= getSettings(gridColumn);

    const globalScaleProp = DG.Property.js('globalScale', DG.TYPE.BOOL, {
      description: 'Determines the way a value is mapped to the vertical scale.\n' +
        '- Global Scale OFF: bottom is column minimum, top is column maximum. Use when columns ' +
        'contain values in different units.\n' +
        '- Global Scale ON: uses the same scale. This lets you compare values ' +
        'across columns, if units are the same (for instance, use it for tracking change over time).'
    });

    const normalizeInput = DG.InputBase.forProperty(globalScaleProp, settings);
    normalizeInput.onChanged.subscribe(() => gridColumn.grid.invalidate());

    return ui.inputs([normalizeInput, ...createBaseInputs(gridColumn, settings)]);
  }
}
