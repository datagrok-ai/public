import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import '../css/admetox.css';
import { TEMPLATES_FOLDER, Model, ModelColoring, Subgroup, DEFAULT_LOWER_VALUE, DEFAULT_UPPER_VALUE, TAGS } from './constants';
import { PieChartCellRenderer } from '@datagrok/power-grid/src/sparklines/piechart';
import { CellRenderViewer } from '@datagrok-libraries/utils/src/viewers/cell-render-viewer';

export let properties: any;

async function getAdmetoxContainer() {
  const admetoxContainer = await grok.dapi.docker.dockerContainers.filter('admetox').first();
  return admetoxContainer;
}

async function startAdmetoxContainer(containerId: string) {
  grok.shell.warning('ADMETox container has not started yet. Try again in a few seconds');
  grok.dapi.docker.dockerContainers.run(containerId);
}

async function sendRequestToContainer(containerId: string, path: string, params: RequestInit): Promise<string | null> {
  try {
    const response = await grok.dapi.docker.dockerContainers.request(containerId, path, params);
    return response;
  } catch (error) {
    //grok.log.error(error);
    return null;
  }
}

export async function runAdmetox(csvString: string, queryParams: string, addProbability: string): Promise<string | null> {
  /*const admetoxContainer = await getAdmetoxContainer();
  if (!admetoxContainer || (admetoxContainer.status !== 'started' && admetoxContainer.status !== 'checking')) {
    await startAdmetoxContainer(admetoxContainer?.id);
    return null;
  }*/

  const params: RequestInit = {
    method: 'POST',
    headers: {
      'Accept': 'text/csv',
      'Content-type': 'text/csv'
    },
    body: csvString
  };

  const path = `http://127.0.0.1:6678/df_upload?models=${queryParams}&probability=${addProbability}`;
  const response = await fetch(path, params);
  return await response.text();
  //return await sendRequestToContainer(admetoxContainer.id, path, params);
}

export async function setProperties(template?: string) {
  if (properties) return;
  const items = await grok.dapi.files.list(TEMPLATES_FOLDER);
  const fileName = template ? `${template}.json` : items[0].fileName;
  const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
  properties = JSON.parse(propertiesJson);
}

export async function performChemicalPropertyPredictions(molColumn: DG.Column, viewTable: DG.DataFrame, properties: string, template?: string, addPiechart?: boolean) {
  await setProperties(template);
  const progressIndicator = DG.TaskBarProgressIndicator.create('Running ADMETox...');
  const csvString = DG.DataFrame.fromColumns([molColumn]).toCsv();
  progressIndicator.update(10, 'Predicting...');

  try {
    const admetoxResults = await runAdmetox(csvString, properties, 'false');
    progressIndicator.update(80, 'Results are ready');
    const table = admetoxResults ? DG.DataFrame.fromCsv(admetoxResults) : null;
    table ? addResultColumns(table, viewTable, addPiechart) : grok.log.warning('');
  } catch (e) {
    grok.log.error(e);
  } finally {
    progressIndicator.close();
  }
}

function applyColorCoding(col: DG.GridColumn, model: Model): void {
  if (!model.coloring) return;
  col.isTextColorCoded = true;
  const { type, min, max, colors } = model.coloring;
  if (type === DG.COLOR_CODING_TYPE.LINEAR)
    col!.column!.meta.colors.setLinear(JSON.parse(colors!), {min: min, max: max});
  else if (type === 'Conditional')
    col!.column!.meta.colors.setConditional(createConditionalColoringRules(model.coloring));
}

export function addColorCoding(table: DG.DataFrame, columnNames: string[]): void {
  const tv = grok.shell.tableView(table.name);
  if (!tv) return;

  for (const columnName of columnNames) {
    const col = tv.grid.col(columnName);
    const model = properties.subgroup.flatMap((subg: Subgroup) => subg.models)
      .find((model: Model) => columnName.includes(model.name));
    if (model) applyColorCoding(col!, model);
  }
}

function createConditionalColoringRules(coloring: ModelColoring): { [index: string]: string | number } {
  const conditionalColors = Object.entries(coloring).slice(1);
  return conditionalColors.reduce((acc, [range, color]) => {
    acc[range] = color;
    return acc;
  }, {} as { [index: string]: string | number });
}

export async function getQueryParams(): Promise<string> {
  await setProperties();
  return properties.subgroup.flatMap((subg: Subgroup) => subg.models)
      .map((model: Model) => model.name).join(',');
}

function generateNumber(): number {
  return Math.random() * (1 - Number.MIN_VALUE) + Number.MIN_VALUE;
}

function createPieSettings(table: DG.DataFrame, columnNames: string[], properties: any): any {
  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'];
  let sectors: any[] = [];
  let sectorColorIndex = 0;

  for (const subgroup of properties.subgroup) {
    const sector: any = {
      name: subgroup.name,
      sectorColor: colors[sectorColorIndex],
      subsectors: []
    };
    
    for (const model of subgroup.models) {
      const modelName = columnNames.find((name: string) => name.includes(model.name));
      if (modelName) {
        let { min, max } = model;
        if (model.properties) {
          const directionProperty = model.properties.find((prop: any) => prop.property.name === 'direction');
          if (directionProperty && directionProperty.object.direction === 'Lower is better')
            [min, max] = [max, min];
        }

        const column = table.col(modelName);
        let weight = Number(generateNumber().toFixed(2));

        if (column) {
          min = column.getTag(TAGS.LOW) ? Number(column.getTag(TAGS.LOW)) : min;
          max = column.getTag(TAGS.HIGH) ? Number(column.getTag(TAGS.HIGH)) : max;
          weight = column.getTag(TAGS.WEIGHT) ? Number(column.getTag(TAGS.WEIGHT)) : weight;
          
          column.setTag(TAGS.GROUP_NAME, subgroup.name);
          column.setTag(TAGS.HIGH, max);
          column.setTag(TAGS.LOW, min);
          column.setTag(TAGS.WEIGHT, weight.toString());
          column.setTag(TAGS.SECTOR_COLOR, colors[sectorColorIndex]);
        }
          
        sector.subsectors.push({
          name: modelName,
          low: min,
          high: max,
          weight: weight,
        });
      }
    }

    if (sector.subsectors.length > 0)
      sectors.push(sector);

    sectorColorIndex = (sectorColorIndex + 1) % colors.length;  
  }
    
  return {
    sectors: {
      lowerBound: DEFAULT_LOWER_VALUE,
      upperBound: DEFAULT_UPPER_VALUE,
      values: '',
      sectors
    }
  };
}

export function addSparklines(table: DG.DataFrame, columnNames: string[]): void {
  const tv = grok.shell.tableView(table.name);
  if (!tv) return;

  const pie = tv.grid.columns.add({ cellType: 'piechart' });
  pie.settings = { columnNames: columnNames };
  pie.settings = createPieSettings(table, columnNames, properties);
}

function getTooltipContent(model: any, value: any): string {
  if (!model) return '';
  const rangesProp = model.properties.find((prop: any) => prop.property.name === 'ranges');
  let tooltipContent = '';
  if (model && rangesProp) {
    const ranges = rangesProp.object.ranges;
    for (const rangeKey in ranges) {
      const rangeParts = rangeKey.split(' ');
      const rangeType = rangeParts.includes('-') ? '-' : rangeParts[0];
      const rangeStart = rangeParts.includes('-') ? parseFloat(rangeParts[0]) : parseFloat(rangeParts[1]);
      const rangeEnd = parseFloat(rangeParts[2]);
      if ((rangeType === '-' && value >= rangeStart && value <= rangeEnd) ||
        (rangeType === '<' && value < rangeStart) ||
        (rangeType === '>' && value > rangeStart)) {
        tooltipContent += `${ranges[rangeKey as keyof typeof ranges]}\n`;
        break;
      }
    }
  } else {
    const direction = model.properties.find((prop: any) => prop.property.name === 'direction');
    const interpretation = direction ? direction.object.direction : '';
    tooltipContent += `${interpretation}\n`;
  }
  return tooltipContent;
}

export function addCustomTooltip(table: string): void {
  const view = grok.shell.tableView(table);
  view.grid.onCellTooltip((cell, x, y) => {
    if (cell.isTableCell) {
      const subgroup = cell.tableColumn!.name;
      const value = cell.cell.value;
      const model = properties.subgroup.flatMap((subg: Subgroup) => subg.models)
        .find((model: Model) => subgroup.includes(model.name));
      const tooltipContent = getTooltipContent(model, value);
      ui.tooltip.show(ui.divV([
        ui.divText(tooltipContent)
      ]), x, y);
      return true;
    }
  });
}

function updateColumnProperties(column: DG.Column, model: any, viewTable: DG.DataFrame): void {
  const newColumnName = viewTable.columns.getUnusedName(column.name);
  column.name = newColumnName;
  column.meta.format = '0.00';
  column.setTag(DG.TAGS.DESCRIPTION, model.properties.find((prop: any) => prop.property.name === 'description').object.description);
  column.meta.units = model.units;
}

export function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame, addPiechart: boolean = true): void {
  if (table.columns.length === 0) return;

  if (table.rowCount > viewTable.rowCount)
    table.rows.removeAt(table.rowCount - 1);

  const modelNames: string[] = table.columns.names();
  const updatedModelNames: string[] = [];
  const models = properties.subgroup.flatMap((subgroup: any) => subgroup.models.map((model: any) => model));

  for (let i = 0; i < modelNames.length; ++i) {
    let column: DG.Column = table.columns.byName(modelNames[i]);
    for (const model of models) {
      if (model.name === modelNames[i]) {
        updateColumnProperties(column, model, viewTable);
        break;
      }
    }
    updatedModelNames.push(column.name);
    viewTable.columns.add(column);
  }

  if (addPiechart)
    addSparklines(viewTable, updatedModelNames);

  addColorCoding(viewTable, updatedModelNames);
  addCustomTooltip(viewTable.name);
}

export async function getModelsSingle(smiles: string, semValue: DG.SemanticValue): Promise<DG.Accordion> {
  const acc = ui.accordion('ADME/Tox');
  await setProperties();

  const update = async (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties.subgroup.find((subg: any) => subg.name === modelName)
      ['models'].map((model: any) => model.name);

    if (smiles === 'MALFORMED_INPUT_VALUE') {
      result.appendChild(ui.divText('The molecule is possibly malformed'));
      return;
    }

    result.appendChild(ui.loader());
    try {
      const csvString = await runAdmetox(`smiles\n${smiles}`, queryParams.join(','), 'false');
      ui.empty(result);
      const table = DG.DataFrame.fromCsv(csvString!);
      const map: { [_: string]: any } = {};
      for (const model of queryParams) {
        map[model] = Number(table.col(model)?.get(0)).toFixed(2);
      }
      result.appendChild(ui.tableFromMap(map));
    } catch (e) {
      result.appendChild(ui.divText('Couldn\'t analyse properties'));
      //console.log(e);
    }
  };

  for (const subgroup of properties.subgroup) {
    const result = ui.div();
    acc.addPane(subgroup.name, () => {
      update(result, subgroup.name);
      return result;
    }, false);
  }

  acc.addPane('Summary', () => createSummaryPane(semValue), false);

  return acc;
}

function createSummaryPane(semValue: DG.SemanticValue): HTMLElement {
  const result = ui.div();
  result.append(ui.loader());

  createPieChartPane(semValue).then((canvas) => {
    ui.empty(result);
    result.appendChild(canvas);
  }).catch((error) => {
    ui.empty(result);
    result.appendChild(ui.divText('Error creating pie chart'));
    console.error(error);
  });

  return result;
}

async function createPieChartPane(semValue: DG.SemanticValue): Promise<HTMLElement> {
  const { cell } = semValue;
  const { dataFrame, column, rowIndex, value } = cell;

  const view = grok.shell.tableView(dataFrame.name);
  const gridCol = view.grid.col(column.name);
  const gridCell = view.grid.cell(column.name, rowIndex);

  // Use params when containers will work on dev
  const params = await getQueryParams();
  const query = `smiles\n${value}`;
  const result = await runAdmetox(query, params, 'false');

  const pieSettings = createPieSettings(dataFrame, params.split(','), properties);
  pieSettings.sectors.values = result!;
  gridCol!.settings = pieSettings;

  const pieChartRenderer = new PieChartCellRenderer();
  return CellRenderViewer.fromGridCell(gridCell, pieChartRenderer).root;
}
