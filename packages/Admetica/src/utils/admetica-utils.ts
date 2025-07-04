import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  colorsDictionary,
  DEFAULT_LOWER_VALUE,
  DEFAULT_TABLE_NAME,
  DEFAULT_UPPER_VALUE,
  ERROR_MESSAGES,
  Model,
  ModelColoring,
  Subgroup,
  TEMPLATES_FOLDER
} from './constants';
import {FormStateGenerator} from './admetica-form';
import {CellRenderViewer} from '../viewers/cell-render-viewer';

import '../css/admetica.css';

export let properties: any;
export const tablePieChartIndexMap: Map<string, number> = new Map();
let piechartIndex = 0;


export async function runAdmeticaFunc(csv: string, models: string, raiseException: boolean): Promise<DG.DataFrame> {
  let df: DG.DataFrame = await grok.functions.call('Admetica:run_admetica', {'csv': csv, 'models': models, 'raiseException': raiseException});
  df = await convertLD50(df);
  return df;
}

export async function convertLD50(df: DG.DataFrame): Promise<DG.DataFrame> {
  const smilesCol = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (!df.columns.names().includes('LD50')) return df;

  const ldCol = df.getCol('LD50');
  const molWeights: DG.Column = await grok.functions.call('Chem: getMolProperty', { molecules: smilesCol, property: "MW" });

  ldCol.init((i) => {
    const molPerKg = Math.pow(10, -ldCol.get(i));
    return molPerKg * molWeights.get(i) * 1000;
  });

  return df;
}

export async function setProperties() {
  if (properties) return;
  const items = await grok.dapi.files.list(TEMPLATES_FOLDER);
  const fileName = items[0].fileName;
  const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
  properties = JSON.parse(propertiesJson);
}

export async function performChemicalPropertyPredictions(
  molColumn: DG.Column, viewTable: DG.DataFrame, models: string, template?: string,
  addPiechart: boolean = false, addForm: boolean = false, update: boolean = false, raiseException: boolean = false
) {
  if (template) {
    try {
      properties = JSON.parse(template);
    } catch (error) {
      grok.shell.error('Invalid JSON template.');
      throw new Error('Failed to parse template JSON.');
    }
  }
  const progressIndicator = DG.TaskBarProgressIndicator.create('Running Admetica...');
  try {
    const csvString = DG.DataFrame.fromColumns([molColumn]).toCsv();
    progressIndicator.update(10, 'Predicting...');
    const table = await runAdmeticaFunc(csvString, models, raiseException);
    progressIndicator.update(80, 'Results are ready');
    const molColIdx = viewTable.columns.names().findIndex(name => name === molColumn.name);
    if (table) {
      addResultColumns(table, viewTable, addPiechart, addForm, molColIdx ?? -1, update);
    }
  } catch (error) {
    if (raiseException) throw error;
    grok.log.error(error);
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
  const tableView = grok.shell.getTableView(table.name);
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

function createPieSettings(table: DG.DataFrame, columnNames: string[], properties: any): any {
  let sectors: any[] = [];

  for (const subgroup of properties.subgroup) {
    const subgroupColor = colorsDictionary[subgroup.name];
    const sector: any = {
      name: subgroup.name,
      sectorColor: subgroupColor,
      subsectors: []
    };

    for (const model of subgroup.models) {
      const modelName = columnNames.find((name: string) => name.includes(model.name));
      if (modelName) {
        const column = table.col(modelName);
        const { line, weight, min, max } = model;

        if (column) {
          const updatedMeta = {
            groupName: subgroup.name,
            weight: weight,
            line: line,
            sectorColor: subgroupColor,
            min: min,
            max: max
          };

          column.setTag('.vlaaivis-metadata', JSON.stringify(updatedMeta));
        }

        sector.subsectors.push({
          name: modelName,
          weight: weight,
          line: line,
        });
      }
    }

    if (sector.subsectors.length > 0)
      sectors.push(sector);
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

export function addSparklines(table: DG.DataFrame, columnNames: string[], index: number, name?: string): void {
  const tv = grok.shell.getTableView(table.name);;
  const { grid } = tv;
  if (!tv) return;

  let pieChartIdx: number;
  if (tablePieChartIndexMap.has(table.name)) {
    pieChartIdx = tablePieChartIndexMap.get(table.name)! + 1;
    tablePieChartIndexMap.set(table.name, pieChartIdx);
  } else {
    pieChartIdx = piechartIndex + 1;
    tablePieChartIndexMap.set(table.name, pieChartIdx);
  }

  const pieName = pieChartIdx === 0 ? "piechart" : `piechart (${pieChartIdx})`
  name ??= pieName;
  const pie = grid.columns.add({ gridColumnName: name, cellType: 'piechart' });

  pie.settings = { columnNames: columnNames };
  pie.settings = createPieSettings(table, columnNames, properties);

  grid.columns.byName(name)?.move(index);
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
      if (
        (rangeType === '-' && value >= rangeStart && value <= rangeEnd) ||
        (rangeType === '<' && value < rangeStart) ||
        (rangeType === '<=' && value <= rangeStart) ||
        (rangeType === '>' && value > rangeStart) ||
        (rangeType === '>=' && value >= rangeStart)
      ) {
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

export function addCustomTooltip(table: DG.DataFrame): void {
  const view = grok.shell.getTableView(table.name);
  view.grid.onCellTooltip((cell, x, y) => {
    if (cell.isTableCell && typeof cell.cell.value === 'number') {
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

function setAdmeGroups(table: DG.DataFrame, columnNames: string[]): void {
  const isColumnRelatedToSubgroup = (column: string, subgroup: Subgroup): boolean => {
    return subgroup.models.some(model => column.includes(model.name));
  };

  const createSubgroupDict = (properties: any, columnNames: string[]) => {
    return properties.subgroup.reduce((dict: { [key: string]: { color: string; columns: string[] } }, subgroup: Subgroup) => {
      dict[subgroup.name] = {
        color: colorsDictionary[subgroup.name],
        columns: columnNames.filter(column => isColumnRelatedToSubgroup(column, subgroup))
      };
      return dict;
    }, {});
  };

  const subgroupDict = createSubgroupDict(properties, columnNames);
  // Temporary workaround: Addresses the issue where setGroups does not trigger a refresh 
  // when adding to an existing grid. This fix will be removed once the issue is resolved.
  table.temp['.columnGroups'] = subgroupDict;
}

export function updateColumnProperties(gridCol: DG.GridColumn, model: any): void {
  if (!gridCol.column) return;

  const column = gridCol.column;
  column.meta.format = '0.000';

  const description = model.properties.find((prop: any) => prop.property.name === 'description')?.object.description;
  if (description)
    column.setTag(DG.TAGS.DESCRIPTION, description);

  column.meta.units = model.units;

  const subgroupName = properties.subgroup.find((subg: Subgroup) =>
    subg.models.some((m: Model) => model.name.includes(m.name))
  )?.name;
  if (subgroupName) {
    gridCol.headerCellStyle.textColor = DG.Color.fromHtml(colorsDictionary[subgroupName]);
    column.setTag('group', subgroupName);
  }
}

export function addResultColumns(
  table: DG.DataFrame,
  viewTable: DG.DataFrame,
  addPiechart: boolean = true,
  addForm: boolean = true,
  molColIdx: number,
  update: boolean = false
): void {
  if (table.columns.length === 0) return;

  if (table.rowCount > viewTable.rowCount) {
    table.rows.removeAt(table.rowCount - 1);
  }

  const modelNames: string[] = table.columns.names();
  const updatedModelNames: string[] = [];
  const models = properties.subgroup.flatMap((subgroup: Subgroup) => subgroup.models);

  for (const modelName of modelNames) {
    const column: DG.Column = table.columns.byName(modelName);
    const newColumnName = viewTable.columns.getUnusedName(column.name);
    column.name = newColumnName;

    const colToReplace = viewTable.col(modelName);
    const model = models.find((m: any) => m.name === modelName);

    if (model) {
      if (!update || !colToReplace)
        viewTable.columns.add(column);
      else
        viewTable.columns.replace(colToReplace, column);
      updatedModelNames.push(column.name);
    }
  }

  const tableView = grok.shell.getTableView(viewTable.name);
  const { grid } = tableView;
  models.forEach((model: Model) => {
    const columnName = updatedModelNames.find((name) => name.includes(model.name));
    if (!columnName) return;

    const column = grid.col(columnName);
    if (column) updateColumnProperties(column, model);
  });

  if (addPiechart)
    addSparklines(viewTable, updatedModelNames, molColIdx + 2);

  setAdmeGroups(viewTable, updatedModelNames);
  addColorCoding(viewTable, updatedModelNames);
  addCustomTooltip(viewTable);

  if (addForm) {
    const form = createDynamicForm(viewTable, updatedModelNames, viewTable.columns.names()[molColIdx], addPiechart);
    tableView.dockManager.dock(form, DG.DOCK_TYPE.RIGHT, null, 'Form', 0.45);
    grid.invalidate();
  }
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

    if (smiles === 'MALFORMED_INPUT_VALUE')
      return result.appendChild(ui.divText(ERROR_MESSAGES.MALFORMED));

    if (!smiles || DG.chem.Sketcher.isEmptyMolfile(smiles))
      return result.appendChild(ui.divText(ERROR_MESSAGES.EMPTY));

    result.appendChild(ui.loader());
    try {
      const table = await runAdmeticaFunc(`smiles\n${smiles}`, queryParams.join(','), false);
      ui.empty(result);
      table.name = DEFAULT_TABLE_NAME;
      addColorCoding(table, queryParams, true, props);

      const map: { [_: string]: any } = {};
      for (const param of queryParams) {
        const column = table.getCol(param);
        const firstValue = column.get(0);
        const color = column.meta.colors.getColor(0);

        const model: Model | undefined = properties.subgroup
          .flatMap((subg: Subgroup) => subg.models)
          .find((model: Model) => param.includes(model.name));

        if (model) {
          const units = model.units && model.units !== '-' ? model.units : '';
          const value = typeof firstValue === "string" ? firstValue : `${firstValue.toFixed(3)} ${units}`;
          map[param] = ui.divText(value, {
            style: { color: color ? DG.Color.toHtml(color) : 'black' }
          });
        }
      }
      result.appendChild(ui.tableFromMap(map));
    } catch (e) {
      ui.empty(result);
      result.appendChild(ui.divText('Couldn\'t analyze properties'));
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
  let { cell, units } = semValue;
  cell ??= grok.shell.o.cell;

  if (!cell)
    return ui.divText('Failed to identify the current object');

  const { dataFrame, column, rowIndex, value } = cell;

  const view = grok.shell.getTableView(dataFrame.name);
  const gridCol = view.grid.col(column.name);
  const gridCell = view.grid.cell(column.name, rowIndex);

  const parsedValue = units === DG.UNITS.Molecule.MOLBLOCK ? `"${value}"` : value;
  const params = await getQueryParams();
  const query = `smiles\n${parsedValue}`;
  const result: DG.DataFrame = await runAdmeticaFunc(query, params, false);
  const pieSettings = createPieSettings(dataFrame, params.split(','), properties);
  pieSettings.sectors.values = result.toCsv()!;
  gridCol!.settings = pieSettings;

  const pieChartRenderer = await grok.functions.call('PowerGrid:piechartCellRenderer');
  return CellRenderViewer.fromGridCell(gridCell, pieChartRenderer).root;
}

export function createDynamicForm(viewTable: DG.DataFrame, updatedModelNames: string[], molColName: string, addPiechart: boolean) {
  const form = DG.FormViewer.createDefault(viewTable, { columns: updatedModelNames });
  const mapping = FormStateGenerator.createCategoryModelMapping(properties, updatedModelNames);
  const generator = new FormStateGenerator(viewTable.name, mapping, molColName, addPiechart);
  const formState = generator.generateFormState();
  form.form.state = JSON.stringify(formState);

  ui.onSizeChanged(form.root).subscribe(() => {
    const containerWidth = form.root.clientWidth;
    if (!containerWidth) return;
    const updatedFormState = generator.generateFormState(containerWidth);
    form.form.state = JSON.stringify(updatedFormState);
  });
  return form;
}

export async function getTemplates(): Promise<string[]> {
  const files = await grok.dapi.files.list(TEMPLATES_FOLDER);
  return files.map((file) => file.fileName.split('.')[0]);
}
