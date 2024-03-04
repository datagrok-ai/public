import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { properties, models, template } from './const';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';
import $ from 'cash-dom';
import '../css/admetox.css';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
let propertiesNew: { [s: string]: unknown; } | ArrayLike<unknown>;

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
    /*const model = Object.values(propertiesNew)
      .flatMap((category: any) => category.models)
      .find((m: any) => columnName.includes(m.name));*/
    
    //@ts-ignore
    const model = template.subgroup.models.find(model => columnName.includes(model.name));

    //if (!model!.coloring || model.coloring === "") continue;

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
  const queryParams = getQueryParams();
  try {
    await performChemicalPropertyPredictions(molCol, viewTable, queryParams, false);
  } catch (e) {
    //grok.log.error(e);
  }
}

export function getQueryParams(): string {
  return Object.keys(properties)
    .flatMap(property => properties[property]['models'])
    .filter(obj => obj['skip'] !== true)
    .map(obj => obj['name'])
    .join(',');
}

function normalizeValue(value: number, minValue: number, maxValue: number): number {
  return (value - minValue) / (maxValue - minValue);
}

function minMaxNormalization(value: number, oldMin: number, oldMax: number, newMin: number = 0, newMax: number = 10): number {
  return ((value - oldMin) * (newMax - newMin)) / (oldMax - oldMin) + newMin;
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
      const model = Object.values(propertiesNew).flatMap((category: any) => category.models).find((m: any) => column.name.includes(m.name));
      const normalized = normalizeValue(column.get(i), column.stats.min, column.stats.max);
      if (model.preference === 'higher')
        return normalized;
      return 1 - normalized;
    });
    newCol.setTag(DG.Tags.IncludeInCsvExport, 'false');
    newCol.setTag(DG.Tags.IncludeInBinaryExport, 'false');
  }
  const sparkline = tv.grid.columns.add({cellType: 'sparkline'});
  sparkline.settings = {columnNames: columnNames};
  const radar = tv.grid.columns.add({cellType: 'radar'});
  radar.settings = {columnNames: columnNames};
  const bar = tv.grid.columns.add({cellType: 'barchart'});
  bar.settings = {columnNames: columnNames};
  const pie = tv.grid.columns.add({cellType: 'piechart'});
  pie.settings = {columnNames: columnNames};
}

function addCustomTooltip(table: string) {
  const view = grok.shell.tableView(table);
  view.grid.onCellTooltip(function (cell, x, y) {
    if (cell.isTableCell) {
      const subgroup = cell.tableColumn!.name;
      const value = cell.cell.value;
      //@ts-ignore
      const model = template.subgroup.models.find(model => subgroup.includes(model.name));
      //const model = Object.values(propertiesNew).flatMap((category: any) => category.models).find((m: any) => subgroup.includes(m.name));
      if (!model) return;
      const interpretation = model;
  
      let tooltipContent = '';
      //@ts-ignore
      if (interpretation && interpretation.ranges) {
        //@ts-ignore
        for (let rangeKey in interpretation.ranges) {
          const rangeParts = rangeKey.split(' ');
          const rangeType = rangeParts.includes('-') ? '-' : rangeParts[0];
          const rangeStart = rangeParts.includes('-') ? parseFloat(rangeParts[0]) : parseFloat(rangeParts[1]);
          const rangeEnd = parseFloat(rangeParts[2]);
          if ((rangeType === '-' && value >= rangeStart && value <= rangeEnd) ||
              (rangeType === '<' && value < rangeStart) ||
              (rangeType === '>' && value > rangeStart)) {
                //@ts-ignore
            tooltipContent += `${interpretation.ranges[rangeKey as keyof typeof interpretation.ranges]}\n`;
            break;
          }
        }
      } else
        tooltipContent += `${interpretation}\n`;
    
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

  for (let i = 0; i < modelNames.length; ++i) {
    let column: DG.Column = table.columns.byName(modelNames[i]);
    const newColumnName = viewTable.columns.getUnusedName(modelNames[i]);
    column.name = newColumnName;
    column.setTag(DG.TAGS.FORMAT, '0.00');

    for (const key in models) {
      if (modelNames[i].includes(key)) {
        column.setTag(DG.TAGS.DESCRIPTION, models[key]);
        column.setTag(DG.TAGS.UNITS, getModelUnits(key));
        break;
      }
    }

    updatedModelNames[i] = newColumnName;
    viewTable.columns.add(column);
  }

  addColorCoding(viewTable, updatedModelNames);
  //addNormalizedColumns(viewTable, updatedModelNames);
  addCustomTooltip(viewTable.name);
}

export function getModelUnits(modelName: string): string {
  //@ts-ignore
  //const modelCategories = Object.values(propertiesNew['subgroup']);
  //const model = modelCategories.flatMap((category: any) => category.models).find((m: any) => m.name === modelName);
  const model = template.subgroup.models.find(model => modelName.includes(model.name));
  //@ts-ignore
  return model!.units!;
};

export function getModelsSingle(smiles: string): DG.Accordion {
  const acc = ui.accordion('ADME/Tox');

  const update = async (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties[modelName]['models']
      .filter((model: any) => !model.skip)
      .map((model: any) => model.name);

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

  for (const property of Object.keys(properties)) {
    const models = properties[property]['models'];
    const shouldAddProperty = models.some((model: any) => !model.skip);

    if (shouldAddProperty) {
      const result = ui.div();
      acc.addPane(property, () => {
        update(result, property);
        return result;
      }, false);
    }
  }

  return acc;
}

async function openModelsDialog(selected: any, viewTable: DG.DataFrame, onOK: any): Promise<void> {
  //const jsonData = await _package.files.readAsText('template.json');
  const items = (await grok.dapi.files.list('System:AppData/Admetox/templates')).map((file) => file.fileName.split('.')[0]);
  const templates = ui.choiceInput('Template', items[0], items);
  let templateCopy = template; 
  const host = ui.div([]);
  const tabsV = ui.tabControl(null, true);
  templateCopy.subgroup.forEach((subgroup) => {
    subgroup.models.forEach((model) => {
      const inputs = ui.divV([]);
      inputs.classList.add('admetox-input-form');
      const properties = model.properties;
      properties.map((p) => {
        let object: any = {};
        object[p.name] = p.defaultValue;
        const prop = DG.Property.fromOptions(p);
        const input = DG.InputBase.forProperty(prop, object);
        input.enabled = p.enable;
        inputs.appendChild(input.root);
      });
      tabsV.addPane(model.name, () => inputs);
    });
  });

  tabsV.panes.forEach((pane)=> {
    const functionCheck = ui.boolInput('', false//, (v:boolean) => {
     // calculatedFunctions[funcNamesMap[pane.name]] = !!functionCheck.value;
      /*if (!v)
        $(pane.content).find('input').attr('disabled', 'true');
      else
        $(pane.content).find('input').removeAttr('disabled');*/
    /*}*/);
    pane.header.appendChild(functionCheck.root);
    pane.header.classList.add('admetox-pane-header');
  });
  /*templateCopy.subgroup.models.forEach((model) => {
    const inputs = ui.divV([]);
    Object.entries(model.ranges).map(([range, meaning]) => {
      const input = ui.divH([
        ui.stringInput('', range, (value: string) => {
          const newRanges = { ...model.ranges }; // Create a shallow copy of the ranges object
          const meaningValue = newRanges[range as keyof typeof model.ranges]; // Store the meaning value of the old range
          delete newRanges[range as keyof typeof model.ranges]; // Delete the old range
          newRanges[value as keyof typeof model.ranges] = meaningValue; // Add the new range with the old meaning value
          model.ranges = newRanges;
          range = value;
        }).root,
        ui.stringInput('', meaning, (value: string) => {
          const r = range as keyof typeof model.ranges;
          model.ranges[r] = value;
        }).root
      ]);
      inputs.appendChild(input);
    });
    inputs.appendChild(ui.boolInput('Desirable', true).root);
    tabsV.addPane(model.name, () => inputs);
  });*/
  host.appendChild(tabsV.root);

  
  /*const tree = ui.tree();
  tree.root.classList.add('admetox-dialog-tree');

  const groups: { [_: string]: any } = {};
  const items: DG.TreeViewNode[] = [];
  const selectedModels: { [_: string]: string } = {};

  const checkAll = (val: boolean) => {
    for (const g of Object.values(groups))
      g.checked = val;
    for (const i of items)
      i.checked = val;
  };

  const selectAll = ui.label('All', {classes: 'd4-link-label', onClick: () => checkAll(true)});
  selectAll.classList.add('admetox-dialog-select-all');
  const selectNone = ui.label('None', {classes: 'd4-link-label', onClick: () => checkAll(false)});
  const countLabel = ui.label('0 checked');
  countLabel.classList.add('admetox-dialog-count');

  const keys = Object.keys(properties);
  for (const groupName of keys) {
    const group = tree.group(groupName, null, false);
    group.enableCheckBox();
    groups[groupName] = group;

    group.checkBox!.onchange = (_e) => {
      countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
      if (group.checked) 
        selectedModels[group.text] = group.text;
      group.items.filter((i) => {
        if (i.checked) 
          selectedModels[i.text] = group.text;
      })
    };

    for (const property of properties[groupName]['models']) {
      if (property['skip'] === false) {
        const item = group.item(property['name'], property);
        item.enableCheckBox(selected.includes(property['name']));
        items.push(item);
        
        item.checkBox!.onchange = (_e) => {
          countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
          if (item.checked) 
            selectedModels[item.text] = groupName;
        };
      }
    }

    if (group.items.length === 0) 
      group.remove(); 
    
    checkAll(false);
  }

  const saveInputHistory = (): any => {
    let resultHistory: { [_: string]: any } = {};
    const modelNames = Object.keys(selectedModels);
    for (const modelName of modelNames) 
      resultHistory[modelName] = selectedModels[modelName];
    return resultHistory;
  }

  const loadInputHistory = (history: any): void => {
    checkAll(false);
    const keys: string[] = Object.keys(history);
    for (const key of keys) {
      groups[history[key]].items.filter(function (i: any) {
        if (i.text === key) 
          i.checked = true;
      })
      if (key === history[key])
        groups[history[key]].checked = true;
    }
    countLabel.textContent = `${keys.length} checked`;
  }*/
  
  let smilesCol = viewTable.columns.bySemType(DG.SEMTYPE.MOLECULE);
  /*const molInput = ui.columnInput('Molecules', viewTable, smilesCol, async (col: DG.Column) => {
    smilesCol = col;
  }, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE} as ColumnInputOptions);
  molInput.root.classList.add('admetox-mol-input');
  const boolInput = ui.boolInput('Probability', false);
  boolInput.root.classList.add('admetox-bool-input');*/

  const dlg = ui.dialog('ADME/Tox');
  dlg
    .add(templates)
    .add(host)
    //.add(molInput)
    //.add(ui.divH([selectAll, selectNone, countLabel]))
    //.add(tree.root)
    //.add(boolInput)
    .onOK(() => {
      propertiesNew = templateCopy;
      onOK(['Caco2','PPBR'], false, smilesCol);
    })
    .show()
    //.history(
      //() => saveInputHistory(),
      //(x) => loadInputHistory(x) 
    //);
}

async function getSelected() : Promise<any> {
  const str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  return selected;
}