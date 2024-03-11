import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';
import '../css/admetox.css';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
export let properties: any;

export async function runAdmetox(csvString: string, queryParams: string, addProbability: string): Promise<string | null> {
  const admetoxContainer = await grok.dapi.docker.dockerContainers.filter('admetox').first();
  if (admetoxContainer.status !== 'started' && admetoxContainer.status !== 'checking') {
    grok.shell.warning('ADMETox container has not started yet. Try again in a few seconds');
    grok.dapi.docker.dockerContainers.run(admetoxContainer.id);
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
  try {
    const response = await grok.dapi.docker.dockerContainers.request(admetoxContainer.id, path, params);
    return response;
  } catch (error) {
    //grok.log.error(error);
    return null;
  }
}

async function setProperties() {
  const items = (await grok.dapi.files.list('System:AppData/Admetox/templates')).map((file) => file.fileName.split('.')[0]);
  const propertiesJson = await grok.dapi.files.readAsText(`System:AppData/Admetox/templates/${items[0]}.json`);
  properties = properties ? properties : JSON.parse(propertiesJson);
}

export async function addCalculationsToTable(viewTable: DG.DataFrame) {
  openModelsDialog(await getSelected(), viewTable, async (selected: any, includeProbabilities: boolean, smilesCol: DG.Column) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    selected = await getSelected();
    await performChemicalPropertyPredictions(smilesCol, viewTable, selected.join(','), includeProbabilities);
  });
}

export async function performChemicalPropertyPredictions(molColumn: DG.Column, viewTable: DG.DataFrame, properties: string, includeProbabilities: boolean) {
  const progressIndicator = DG.TaskBarProgressIndicator.create('Running ADMETox...');
  const smilesColumn = await extractSmilesColumn(molColumn);
  const csvString = DG.DataFrame.fromColumns([smilesColumn]).toCsv();
  progressIndicator.update(10, 'Predicting...');

  try {
    const admetoxResults = await runAdmetox(csvString, properties, String(includeProbabilities));
    progressIndicator.update(80, 'Results are ready');
    const table = admetoxResults ? DG.DataFrame.fromCsv(admetoxResults) : null;
    table ? addResultColumns(table, viewTable) : grok.log.warning('');
  } catch (e) {
    //grok.log.error(e);
  } finally {
    progressIndicator.close();
  }
}

async function extractSmilesColumn(molColumn: DG.Column): Promise<DG.Column> {
  const isSmiles = molColumn?.getTag(DG.TAGS.UNITS) === DG.UNITS.Molecule.SMILES;
  const smilesList: string[] = new Array<string>(molColumn.length);
  for (let rowIndex = 0; rowIndex < molColumn.length; rowIndex++) {
    let el: string = molColumn?.get(rowIndex);
    if (!isSmiles) {
      try {
        el = await grok.functions.call('Chem:convertMolNotation', {
          molecule: el,
          sourceNotation: DG.chem.Notation.Unknown,
          targetNotation: DG.chem.Notation.Smiles
        });
      } catch {
        el = '';
      }
    }

    smilesList[rowIndex] = el;
  }
  const smilesColumn: DG.Column = DG.Column.fromStrings('smiles', smilesList);
  return smilesColumn;
}

export function addColorCoding(table: DG.DataFrame, columnNames: string[]) {
  const tv = grok.shell.tableView(table.name);
  if (!tv)
    return;

  for (const columnName of columnNames) {
    const col = tv.grid.col(columnName);
    const model = properties.subgroup.flatMap((subg: any) => subg.models)
      .find((model: any) => columnName.includes(model.name));

    if (!model!.coloring || model.coloring === "") continue;

    if (!col) continue;

    col.isTextColorCoded = true;
    if (model!.coloring!.type === 'Linear') {
      const min = model!.coloring!.min;
      const max = model!.coloring!.max;
      const colors = model!.coloring!.colors;
      col.column!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
      col.column!.tags[DG.TAGS.COLOR_CODING_LINEAR] = colors;
      col.column!.tags[DG.TAGS.COLOR_CODING_SCHEME_MIN] = min;
      col.column!.tags[DG.TAGS.COLOR_CODING_SCHEME_MAX] = max;
    } else if (model!.coloring!.type === 'Conditional') {
      const conditionalColors = Object.entries(model!.coloring!).slice(1);
      const conditionalColoring = `{${conditionalColors.map(([range, color]) => `"${range}":"${color}"`).join(",")}}`;
      col.column!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
      col.column!.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = conditionalColoring;
    }
  }
}

export async function addAllModelPredictions(molCol: DG.Column, viewTable: DG.DataFrame) {
  const queryParams = await getQueryParams();
  try {
    await performChemicalPropertyPredictions(molCol, viewTable, queryParams, false);
  } catch (e) {
    //grok.log.error(e);
  }
}

export async function getQueryParams(): Promise<string> {
  await setProperties();
  return properties.subgroup.flatMap((subg: any) => subg.models)
      .map((model: any) => model.name).join(',');
}

function normalizeValue(value: number, minValue: number, maxValue: number): number {
  return (value - minValue) / (maxValue - minValue);
}

function addNormalizedColumns(table: DG.DataFrame, columnNames: string[]) {
  const tv = grok.shell.tableView(table.name);
  for (let i = 0; i < columnNames.length; ++i) {
    const column = table.columns.byName(columnNames[i]);
    const name = table.columns.getUnusedName(columnNames[i]);
    const savedColumnName = '~' + name;
    columnNames[i] = savedColumnName;
    const newCol = table.columns.getOrCreate(savedColumnName, DG.TYPE.FLOAT, table.rowCount);
    newCol.init((i) => {
      const model = properties.subgroup.flatMap((subg: any) => subg.models)
      .find((model: any) => column.name.includes(model.name));
      const normalized = normalizeValue(column.get(i), column.stats.min, column.stats.max);
      if (model.preference === 'higher')
        return normalized;
      return 1 - normalized;
    });
    newCol.setTag(DG.Tags.IncludeInCsvExport, 'false');
    newCol.setTag(DG.Tags.IncludeInBinaryExport, 'false');
  }
  const bar = tv.grid.columns.add({cellType: 'barchart'});
  bar.settings = {columnNames: columnNames};
}

function addCustomTooltip(table: string) {
  const view = grok.shell.tableView(table);
  view.grid.onCellTooltip(function (cell, x, y) {
    if (cell.isTableCell) {
      const subgroup = cell.tableColumn!.name;
      const value = cell.cell.value;
      const model = properties.subgroup.flatMap((subg: any) => subg.models)
        .find((model: any) => subgroup.includes(model.name));
      if (!model) return;
      const rangesProp = model.properties.find((prop: any) => prop.property.name === 'ranges');

      let tooltipContent = '';
      if (model && rangesProp) {
        const ranges = rangesProp.object.ranges;
        for (let rangeKey in ranges) {
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

      ui.tooltip.show(ui.divV([
        ui.divText(tooltipContent)
      ]), x, y);
      return true;
    }
  });
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame) {
  if (table.columns.length === 0)
    return;

  if (table.rowCount > viewTable.rowCount)
    table.rows.removeAt(table.rowCount - 1);

  const modelNames: string[] = table.columns.names();
  const updatedModelNames: string[] = [];
  const models = properties.subgroup.flatMap((subgroup: any) => subgroup.models.map((model: any) => model));

  for (let i = 0; i < modelNames.length; ++i) {
    let column: DG.Column = table.columns.byName(modelNames[i]);
    const newColumnName = viewTable.columns.getUnusedName(modelNames[i]);
    column.name = newColumnName;
    column.setTag(DG.TAGS.FORMAT, '0.00');

    for (const model of models) {
      if (model.name === modelNames[i]) {
        column.setTag(DG.TAGS.DESCRIPTION, model.properties.find((prop: any) => prop.property.name === 'description').object.description);
        column.setTag(DG.TAGS.UNITS, model.properties.find((prop: any) => prop.property.name === 'units').object.units);
        break;
      }
    }

    updatedModelNames[i] = newColumnName;
    viewTable.columns.add(column);
  }

  addColorCoding(viewTable, updatedModelNames);
  addNormalizedColumns(viewTable, updatedModelNames);
  addCustomTooltip(viewTable.name);
}

export async function getModelsSingle(smiles: string): Promise<DG.Accordion> {
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

  return acc;
}

function createTabControl(template: any, selected: string[]) {
  let tabsV = ui.tabControl(null, true);
  template.subgroup.forEach((subgroup: any) => {
    subgroup.models.forEach((model: any) => {
      const inputs = ui.divV([]);
      inputs.classList.add('admetox-input-form');
      const properties = model.properties;
      properties.map((p: any) => {
        const object = p.property.inputType === DG.InputType.Map ? {} : p.object;
        const property = p.property;

        const prop = DG.Property.fromOptions(property);
        const input = DG.InputBase.forProperty(prop, object);
        input.enabled = property.enable;
        const key = property.name as keyof typeof p.object;;
        input.value = p.object[key];
        input.addCaption(property.name);
        input.onChanged(() => {
          p.object[key] = input.value;
          grok.shell.info(input.value);
        });
        inputs.appendChild(input.root);
      });
      tabsV.addPane(model.name, () => inputs);
    });
  });

  tabsV.panes.forEach((pane)=> {
    const functionCheck = ui.boolInput('', false, (v:boolean) => {
      if (v)
        selected.push(pane.name);
      else
        selected.splice(selected.indexOf(pane.name), 1);
    });
    pane.header.insertBefore(functionCheck.root, pane.header.firstChild);
    pane.header.classList.add('admetox-pane-header');
  });
  return {tabsV, selected};
}

async function openModelsDialog(selected: any, viewTable: DG.DataFrame, onOK: any): Promise<void> {
  let selectedItems: string[] = [];
  await setProperties();
  const items = (await grok.dapi.files.list('System:AppData/Admetox/templates')).map((file) => file.fileName.split('.')[0]);
  const result = createTabControl(properties, selectedItems);
  const tabsV = result.tabsV;
  selectedItems = result.selected;
  const templates = ui.choiceInput('Template', items[0], items);
  let templateCopy = properties;
  const loadButton = ui.button('Load', () => {
    DG.Utils.download('template.json', JSON.stringify(templateCopy));
  });
  const buttons = ui.divH([loadButton]);
  const host = ui.div([]);
  
  host.appendChild(tabsV.root);

  let smilesCol = viewTable.columns.bySemType(DG.SEMTYPE.MOLECULE);
  const molInput = ui.columnInput('Molecules', viewTable, smilesCol, async (col: DG.Column) => {
    smilesCol = col;
  }, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE} as ColumnInputOptions);

  const dlg = ui.dialog('ADME/Tox');
  dlg
    .add(molInput)
    .add(templates)
    .add(host)
    .add(buttons)
    .onOK(() => {
      properties = templateCopy;
      onOK(selectedItems, false, smilesCol);
    })
    .show()
}

async function getSelected() : Promise<any> {
  const str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  return selected;
}