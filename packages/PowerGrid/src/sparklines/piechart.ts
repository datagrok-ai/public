import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {
  getSettingsBase,
  names,
  SparklineType,
  SummarySettingsBase,
  createTooltip,
  Hit,
  isSummarySettingsBase
} from './shared';
import { generalProps, groupProps, subGroupProps } from './properties';

let minRadius: number;

enum PieChartStyle {
  Radius = 'Radius',
  Angle = 'Angle'
}

interface Subsector {
  name: string;
  lowThreshold: number;
  highThreshold: number;
  weight: number;
  applicability?: number;
  probabilities?: number[];
}

interface PieChartSettings extends SummarySettingsBase {
  radius: number;
  style: PieChartStyle.Radius | PieChartStyle.Angle;
  sectors?: {
    lowerBound: number;
    upperBound: number;
    sectors: {
      name: string;
      sectorColor: string;
      subsectors: Subsector[];
    }[];
    values: string | null;
  };
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

function normalizeValue(value: number, subsector: Subsector): number {
  const { lowThreshold, highThreshold } = subsector;
  const isMax = highThreshold > lowThreshold;
  if (isMax ? value < lowThreshold : value > lowThreshold)
    return 0;
  else if (isMax ? value > highThreshold : value < highThreshold)
    return 1;
  else
    return isMax
      ? (value - lowThreshold) / (highThreshold - lowThreshold)
      : (value - highThreshold) / (lowThreshold - highThreshold);
}

function renderDiagonalStripes(
  g: CanvasRenderingContext2D, box: DG.Rect, r: number,
  currentAngle: number, subsectorAngle: number
) {
  const patternSize = r;
  const patternCanvas = ui.canvas();
  patternCanvas.width = patternSize;
  patternCanvas.height = patternSize;
  const patternCtx = patternCanvas.getContext('2d')!;
  patternCtx.strokeStyle = '#535659';
  patternCtx.lineWidth = 0.5;
  const numLines = 15;
  const spacing = patternSize / (numLines + 1);
  for (let i = 1; i <= numLines; i++) {
    const y = i * spacing;
    patternCtx.beginPath();
    patternCtx.moveTo(0, y);
    patternCtx.lineTo(patternSize, y);
    patternCtx.stroke();
  }
  const pattern = g.createPattern(patternCanvas, 'repeat')!;
  g.beginPath();
  g.moveTo(box.midX, box.midY);
  g.arc(box.midX, box.midY, r, currentAngle, currentAngle + subsectorAngle);
  g.closePath();
  g.save();
  g.translate(box.midX, box.midY);
  g.rotate(Math.PI / 6);
  g.translate(-box.midX, -box.midY);
  g.fillStyle = pattern;
  g.fill();
  g.restore();
}

function renderSubsector(
  g: CanvasRenderingContext2D, box: DG.Rect, sectorColor: string,
  sectorAngle: number, currentAngle: number, subsector: Subsector,
  minRadius: number, cols: DG.Column[], row: number,
  sectorWeight: number
): number {
  const normalizedSubsectorWeight = subsector.weight / sectorWeight;
  const subsectorAngle = sectorAngle * normalizedSubsectorWeight;
  let r = Math.max(Math.min(box.width, box.height) / 2, minRadius);
  const subsectorName = subsector.name;
  const subsectorCol = cols.find((col) => col.name === subsectorName);
  let value;
  let erroneous = false;
  if (subsectorCol) {
    value = subsectorCol.get(row);
    //erroneous = subsector.probabilities[row] < subsector.applicability; 
    const normalizedValue = value /*&& !erroneous*/ ? normalizeValue(value, subsector) : 1;
    r = normalizedValue * (Math.min(box.width, box.height) / 2);
    r = Math.max(r, minRadius);
  }
  if (erroneous)
    renderDiagonalStripes(g, box, r, currentAngle, subsectorAngle);
  else {
    g.beginPath();
    g.moveTo(box.midX, box.midY);
    g.arc(box.midX, box.midY, r, currentAngle, currentAngle + subsectorAngle);
    g.closePath();
    g.strokeStyle = DG.Color.toRgb(DG.Color.lightGray);
    g.lineWidth = 0.6;
    g.stroke();
    g.fillStyle = value ? hexToRgbA(sectorColor, 0.6) : 'rgba(255, 255, 255, 1)';
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
    const { sectors } = settings.sectors;
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
    let row: number = gridCell.cell.row.idx;
    let cols = df.columns.byNames(settings.columnNames);
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
      const { lowerBound, upperBound, sectors, values } = settings.sectors;
      cols = values ? Array.from(DG.DataFrame.fromCsv(values).columns) : cols;
      row = values ? 0 : row;
      sectors.sort((a, b) => calculateSectorWeight(b) - calculateSectorWeight(a));

      let currentAngle = 0;
      const totalSectorWeight = sectors.reduce((acc, sector) => acc + calculateSectorWeight(sector), 0);

      for (const sector of sectors) {
        const sectorWeight = calculateSectorWeight(sector);
        const normalizedSectorWeight = sectorWeight / totalSectorWeight;
        const sectorAngle = 2 * Math.PI * normalizedSectorWeight;

        // Render sector
        g.beginPath();
        g.moveTo(box.midX, box.midY);
        g.arc(box.midX, box.midY, Math.min(box.width, box.height) / 2, currentAngle, currentAngle + sectorAngle);
        g.closePath();
        g.fillStyle = hexToRgbA(sector.sectorColor, 0.4);
        g.fill();

        // Render subsectors
        let subsectorCurrentAngle = currentAngle;
        for (const subsector of sector.subsectors)
          subsectorCurrentAngle = renderSubsector(g, box, sector.sectorColor, sectorAngle, subsectorCurrentAngle, subsector, minRadius, cols, row, sectorWeight);

        // Render inner circle representing the range
        g.beginPath();
        g.arc(box.midX, box.midY, lowerBound * (Math.min(box.width, box.height) / 2), currentAngle, currentAngle + sectorAngle);
        g.arc(box.midX, box.midY, upperBound * (Math.min(box.width, box.height) / 2), currentAngle + sectorAngle, currentAngle, true);
        g.closePath();
        g.fillStyle = hexToRgbA(sector.sectorColor, 0.4);
        g.fill();

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

  renderSettings(gc: DG.GridColumn): Element {
    const settings: PieChartSettings = isSummarySettingsBase(gc.settings) ? gc.settings :
      gc.settings[SparklineType.PieChart] ??= getSettings(gc);

    const elementsDiv = ui.div([]);

    const inputs = ui.inputs([
      ui.columnsInput('Сolumns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
      ui.choiceInput('Style', PieChartStyle.Radius, [PieChartStyle.Angle, PieChartStyle.Radius],
        function (value: PieChartStyle) {
          settings.style = value;
          gc.grid.invalidate();
        }),
      ui.choiceInput('Type', 'Standard', ['Standard', 'VlaaiVis'],
        function (value: string) {
          ui.empty(elementsDiv);
          const sectors = initializeSectors(settings.columnNames, grok.shell.tv.dataFrame);
          if (value === 'Vlaavis') {
            settings.sectors = sectors;
            gc.grid.invalidate();
            elementsDiv.appendChild(createTreeGroup(settings, gc));
          } else {
            delete settings.sectors;
            gc.grid.invalidate();
          }
        }
      )
    ]);

    return ui.divV([inputs, elementsDiv]);
  }
}

function generateNumber() {
  return Math.random();
}

function createInputForProperty(property: any, name: string, settings: PieChartSettings, gc: DG.GridColumn) {
  const prop = DG.Property.fromOptions(property.property);
  const input = DG.InputBase.forProperty(prop, {});

  // Find the initial value in settings.sectors
  let value: any;

  for (const sector of settings.sectors!.sectors) {
    if (sector.name === name) {
      value = sector[property.property.name as keyof typeof sector];
      break;
    }
    for (const subsector of sector.subsectors) {
      if (subsector.name === name) {
        value = subsector[property.property.name as keyof typeof subsector];
        break;
      }
    }

    if (name === '') {
      value = (settings.sectors as any)[property.property.name];
    }

    if (value !== undefined) break;
  }

  input.value = value !== undefined ? value : property.object[property.property.name];

  // Handle input changes
  input.onChanged(() => {
    for (const sector of settings.sectors!.sectors) {
      if (sector.name === name) {
        sector[property.property.name as keyof typeof sector] = input.value;
        break;
      }
      for (const subsector of sector.subsectors) {
        if (subsector.name === name) {
          subsector[property.property.name as keyof typeof subsector] = input.value as never;
          break;
        }
      }

      if (name === '') {
        (settings.sectors as any)[property.property.name] = input.value;
      }
    }
    settings.sectors = settings.sectors;
    gc.grid.invalidate();
  });

  return input.root;
}

/**
 5. think of better naming
 6. UI fixes
 */

// Function to create tree group with properties
function createTreeGroup(settings: PieChartSettings, gc: DG.GridColumn): Element {
  const df = grok.shell.tv.dataFrame;
  const tag = '.group-name';
  const tree = ui.tree();
  const inputs = ui.divV([]);

  // Map to store group names and corresponding columns
  const groupMap = new Map<string, DG.Column[]>();
  const propertiesMap = new Map<string, any>();
  const columns = df.columns.byNames(settings.columnNames);

  // Populate the group map
  //Should utilize columnNames from settings to be sync
  columns.forEach((col: DG.Column) => {
    const groupName = col.getTag(tag);
    if (groupName !== null) {
      if (!groupMap.has(groupName)) {
        groupMap.set(groupName, []);
        propertiesMap.set(groupName, {});
      }
      groupMap.get(groupName)!.push(col);
      propertiesMap.set(col.name, {});
    }
  });

  // Create tree nodes
  Array.from(groupMap.keys()).forEach((group) => {
    const groupNode = tree.group(group);
    groupNode.expanded = true;
    groupNode.enableCheckBox(false);

    const columns = groupMap.get(group);
    columns?.forEach((col) => {
      const colNode = groupNode.item(col.name);
      colNode.checked = false;
      colNode.enableCheckBox(false);
    });
  });

  tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
    ui.empty(inputs);
    if (settings.columnNames.includes(node.text)) {
      for (let prop of subGroupProps)
        inputs.appendChild(createInputForProperty(prop, node.text, settings, gc));
    } else {
      for (let prop of groupProps) {
        inputs.appendChild(createInputForProperty(prop, node.text, settings, gc));
        const propName = prop.property.name;
        propertiesMap.get(node.text)[propName] = (prop.object as any)[propName];
      }
    }
  });

  const generalInp = ui.divV([]);
  //routine to create general props
  for (let prop of generalProps)
    generalInp.appendChild(createInputForProperty(prop, '', settings, gc));

  //console.log('settings object');
  //console.log(JSON.stringify(createTreeSettings(groupMap, propertiesMap)));

  return ui.divV([generalInp, ui.divH([tree.root, inputs], 'ui-form')]);
}

const defaultGeneralProps = generalProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record<string, any>);

// Default properties based on groupProps and subGroupProps
const defaultGroupProps = groupProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record<string, any>);

function initializeSectors(columnNames: string[], dataFrame: DG.DataFrame): {
  lowerBound: number;
  upperBound: number;
  sectors: {
    name: string;
    sectorColor: string;
    subsectors: Subsector[];
  }[];
  values: string | null;
} {
  const groupMap = new Map<string, DG.Column[]>();
  const propertiesMap = new Map<string, any>();
  const columns = dataFrame.columns.byNames(columnNames);
  const tag = '.group-name';

  columns.forEach((col: DG.Column) => {
    const groupName = col.getTag(tag);
    if (groupName !== null) {
      if (!groupMap.has(groupName)) {
        groupMap.set(groupName, []);
        propertiesMap.set(groupName, {});
      }
      groupMap.get(groupName)!.push(col);
      propertiesMap.set(col.name, {});
    }
  });

  let sectors: { name: string; sectorColor: string; subsectors: Subsector[] }[] = [];

  groupMap.forEach((columns, groupName) => {
    let subsectors: Subsector[] = [];

    columns.forEach((col) => {
      subsectors.push({
        name: col.name,
        lowThreshold: parseFloat(col.getTag('.low-threshold')) ?? 0,
        highThreshold: parseFloat(col.getTag('.high-threshold')) ?? 1,
        weight: parseFloat(col.getTag('.weight')) ?? generateNumber(),
        //applicability: defaultSubGroupProps["applicability"] ?? DEFAULT_APPLICABILITY_VALUE,
        //probabilities: []
      });
    });

    sectors.push({
      name: groupName,
      sectorColor: defaultGroupProps["sectorColor"],
      subsectors
    });
  });

  return {
    lowerBound: defaultGeneralProps["lowerBound"] ?? DEFAULT_LOWER_VALUE,
    upperBound: defaultGeneralProps["upperBound"] ?? DEFAULT_UPPER_VALUE,
    sectors,
    values: null
  };
}

const DEFAULT_LOWER_VALUE = 0.7;
const DEFAULT_UPPER_VALUE = 0.9;