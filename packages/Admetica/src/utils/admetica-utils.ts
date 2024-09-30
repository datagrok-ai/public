import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { TEMPLATES_FOLDER, Model, ModelColoring, Subgroup, DEFAULT_LOWER_VALUE, DEFAULT_UPPER_VALUE, TAGS, DEFAULT_TABLE_NAME } from './constants';
import { PieChartCellRenderer } from '@datagrok/power-grid/src/sparklines/piechart';
import { CellRenderViewer } from '@datagrok-libraries/utils/src/viewers/cell-render-viewer';
import { fetchWrapper } from '@datagrok-libraries/utils/src/fetch-utils';

import '../css/admetica.css';
import { FormStateGenerator } from './admetica-form';

export let properties: any;
export const tablePieChartIndexMap: Map<string, number> = new Map();
let piechartIndex = 0;

async function getAdmeticaContainer() {
  const admeticaContainer = await grok.dapi.docker.dockerContainers.filter('admetica').first();
  return admeticaContainer;
}

// can be removed, need to add to container info
async function startAdmeticaContainer(containerId: string) {
  grok.shell.warning('Admetica container has not started yet. Try again in a few seconds');
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

export async function runAdmetica(csvString: string, queryParams: string, addProbability: string): Promise<string | null> {
  const admeticaContainer = await getAdmeticaContainer();
  if (!admeticaContainer || (admeticaContainer.status !== 'started' && admeticaContainer.status !== 'checking')) {
    await startAdmeticaContainer(admeticaContainer?.id);
    return null;
  }

  const params: RequestInit = {
    method: 'POST',
    headers: {
      'Accept': 'text/csv',
      'Content-type': 'text/csv'
    },
    body: csvString
  };

  const path = `/df_upload?models=${queryParams}&probability=${addProbability}`;
  const response = await fetchWrapper(() => sendRequestToContainer(admeticaContainer.id, path, params));
  return await convertLD50(response!, DG.Column.fromStrings('smiles', csvString.split('\n').slice(1)));
}

export async function convertLD50(response: string, smilesCol: DG.Column): Promise<string> {
  const df = DG.DataFrame.fromCsv(response);
  if (!df.columns.names().includes('LD50')) return response;

  const ldCol = df.getCol('LD50');
  const rowCount = df.rowCount;

  const molWeights = await Promise.all(
    Array.from({ length: rowCount }, (_, i) => 
      grok.functions.call('Chem:getProperty', { molecule: smilesCol.get(i), prop: "MW" })
    )
  );

  ldCol.init((i) => {
    const molPerKg = Math.pow(10, -ldCol.get(i));
    const mgPerKg = molPerKg * molWeights[i] * 1000;
    return mgPerKg;
  });
  
  return df.toCsv();
}

export async function setProperties() {
  if (properties) return;
  const items = await grok.dapi.files.list(TEMPLATES_FOLDER);
  const fileName = items[0].fileName;
  const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
  properties = JSON.parse(propertiesJson);
}

export async function performChemicalPropertyPredictions(molColumn: DG.Column, viewTable: DG.DataFrame, models: string, template?: string, addPiechart?: boolean, addForm?: boolean) {
  properties = JSON.parse(template!);
  const progressIndicator = DG.TaskBarProgressIndicator.create('Running Admetica...');
  const csvString = DG.DataFrame.fromColumns([molColumn]).toCsv();
  progressIndicator.update(10, 'Predicting...');

  try {
    const admeticaResults = await runAdmetica(csvString, models, 'false');
    progressIndicator.update(80, 'Results are ready');
    const table = admeticaResults ? DG.DataFrame.fromCsv(admeticaResults) : null;
    const molColIdx = viewTable?.columns.names().findIndex((name) => name === molColumn.name);
    table ? addResultColumns(table, viewTable, addPiechart, addForm, molColIdx!) : grok.log.warning('');
  } catch (e) {
    grok.log.error(e);
  } finally {
    progressIndicator.close();
  }
}

function applyColumnColorCoding(column: DG.Column, model: Model): void {
  const { type, min, max, colors } = model.coloring;

  if (type === DG.COLOR_CODING_TYPE.LINEAR)
    column.meta.colors.setLinear(JSON.parse(colors!), { min, max });
  else if (type === 'Conditional')
    column.meta.colors.setConditional(createConditionalColoringRules(model.coloring));
}

export function addColorCoding(table: DG.DataFrame, columnNames: string[], showInPanel: boolean = false, props?: string): void {
  const tableView = grok.shell.tableView(table.name);
  if (!tableView && !showInPanel) return;

  for (const columnName of columnNames) {
    if (tableView) {
      const gridColumn = tableView.grid.col(columnName);
      if (gridColumn)
        gridColumn.isTextColorCoded = true;
    }

    const column = table.getCol(columnName);
    const matchingModel = (props ?? properties).subgroup
      .flatMap((subgroup: Subgroup) => subgroup.models)
      .find((model: Model) => columnName.includes(model.name));

    if (!matchingModel || !matchingModel.coloring) continue;

    applyColumnColorCoding(column, matchingModel);
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
  const colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'];
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

export function addSparklines(table: DG.DataFrame, columnNames: string[], molColIdx: number): void {
  const tv = grok.shell.tableView(table.name);
  if (!tv) return;

  const pie = tv.grid.columns.add({ cellType: 'piechart', index: molColIdx + 1 });
  let pieChartIdx: number;
  if (tablePieChartIndexMap.has(table.name)) {
    pieChartIdx = tablePieChartIndexMap.get(table.name)!;
    tablePieChartIndexMap.set(table.name, pieChartIdx + 1);
  } else {
    pieChartIdx = piechartIndex;
    tablePieChartIndexMap.set(table.name, pieChartIdx + 1);
  }

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
  column.meta.format = '0.000';
  column.setTag(DG.TAGS.DESCRIPTION, model.properties.find((prop: any) => prop.property.name === 'description').object.description);
  column.meta.units = model.units;
}

export function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame, addPiechart: boolean = true, addForm: boolean = true, molColIdx: number): void {
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
    addSparklines(viewTable, updatedModelNames, molColIdx);

  addColorCoding(viewTable, updatedModelNames);
  addCustomTooltip(viewTable.name);

  if (addForm)
    createDynamicForm(viewTable, updatedModelNames, viewTable.columns.names()[molColIdx], addPiechart);
}

export async function getModelsSingle(smiles: string, semValue: DG.SemanticValue): Promise<DG.Accordion> {
  const acc = ui.accordion('Admetica');
  await setProperties();

  const templates = await getTemplates();
  let props: string;
  
  const handleTemplateChange = async (value: string) => {
    props = JSON.parse(await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${value}.json`));
    
    for (const subgroup of properties.subgroup) {
      const pane = acc.getPane(subgroup.name);
      const container = pane.root.children.item(1) as HTMLDivElement;
      ui.empty(container);
      update(container, subgroup.name);
    }
  };
  
  const templatesInput = ui.input.choice('Template', {
    value: templates[0],
    items: templates,
    onValueChanged: handleTemplateChange,
  });
  acc.root.appendChild(templatesInput.root);

  const update = async (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties.subgroup.find((subg: any) => subg.name === modelName)
      ['models'].map((model: any) => model.name);

    if (smiles === 'MALFORMED_INPUT_VALUE') {
      result.appendChild(ui.divText('The molecule is possibly malformed'));
      return;
    }

    result.appendChild(ui.loader());
    try {
      const csvString = await runAdmetica(`smiles\n${smiles}`, queryParams.join(','), 'false');
      ui.empty(result);

      const table = DG.DataFrame.fromCsv(csvString!);
      table.name = DEFAULT_TABLE_NAME;
      addColorCoding(table, queryParams, true, props);

      const map: { [_: string]: any } = {};
      for (const model of queryParams) {
        const column = table.getCol(model);
        map[model] = ui.divText(column.convertTo(DG.TYPE.STRING, '0.000').get(0), {
          style: { color: DG.Color.toHtml(column.meta.colors.getColor(0)!) }
        });
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
  const { dataFrame, column, rowIndex, value } = cell ?? grok.shell.tv.dataFrame.currentCell;

  const view = grok.shell.tableView(dataFrame.name);
  const gridCol = view.grid.col(column.name);
  const gridCell = view.grid.cell(column.name, rowIndex);

  const params = await getQueryParams();
  const query = `smiles\n${value}`;
  const result = await runAdmetica(query, params, 'false');

  const pieSettings = createPieSettings(dataFrame, params.split(','), properties);
  pieSettings.sectors.values = result!;
  gridCol!.settings = pieSettings;

  const pieChartRenderer = new PieChartCellRenderer();
  return CellRenderViewer.fromGridCell(gridCell, pieChartRenderer).root;
}

function createDynamicForm(viewTable: DG.DataFrame, updatedModelNames: string[], molColName: string, addPiechart: boolean) {
  const form = DG.FormViewer.createDefault(viewTable, {columns: updatedModelNames});
  grok.shell.tv.dockManager.dock(form, DG.DOCK_TYPE.RIGHT, null, 'Form', 0.45);
  const mapping = FormStateGenerator.createCategoryModelMapping(properties, updatedModelNames);
  const generator = new FormStateGenerator(viewTable.name, mapping, molColName, addPiechart);
  const formState = generator.generateFormState();
  form.form.state = JSON.stringify(formState);
}

export async function getTemplates(): Promise<string[]> {
  const files = await grok.dapi.files.list(TEMPLATES_FOLDER);
  return files.map((file) => file.fileName.split('.')[0]);
}
