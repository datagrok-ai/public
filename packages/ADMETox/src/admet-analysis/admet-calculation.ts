import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { properties } from './const';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
const _COLUMN_NAME_STORAGE = 'column_names';

/**
 * Sends a POST request to the Admetox server with CSV data and specific query parameters.
 *
 * @param {string} csvString - The CSV data to be sent to the server (structures).
 * @param {string} queryParams - The query parameters to be included in the request URL (model names).
 * @returns {Promise<string | null | undefined>} A promise that resolves to the server response or undefined if an error occurs.
 */
export async function runAdmetox(csvString: string, queryParams: string): Promise<string | null> {
  const admetoxContainer = await grok.dapi.docker.dockerContainers.filter('admetox').first();
  if (admetoxContainer.status !== 'started' && admetoxContainer.status !== 'checking') {
    grok.log.warning('ADMETox container has not started yet.');
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

  const path = `/df_upload?models=${queryParams}&probability=False`;
  try {
    const response = await grok.dapi.docker.dockerContainers.request(admetoxContainer.id, path, params);
    return response;
  } catch (error) {
    console.error('Failed to access the server:', error);
    return null;
  }
}

export async function addCalculations(smilesCol: DG.Column, viewTable: DG.DataFrame) {
  openModelsDialog(await getSelected(), smilesCol, async (selected: any) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    selected = await getSelected();
    //const pi = DG.TaskBarProgressIndicator.create('Prediction...');
    const csvString = DG.DataFrame.fromColumns([smilesCol]).toCsv();
    const admetoxResults = await runAdmetox(csvString, selected.join(','));
    const table = admetoxResults ? DG.DataFrame.fromCsv(admetoxResults) : null;
    table ? addResultColumns(table, viewTable) : grok.log.warning('');
  });
}

/**
 * Applies color coding to the specified columns in the current TableView.
 * The color coding is based on a conditional rule.
 *
 * @param {string[]} columnNames - An array of column names to which color coding will be applied.
 */
export function addColorCoding(columnNames: string[]) {
  const tv = grok.shell.tv;
  for (const columnName of columnNames) {
    tv.grid.col(columnName)!.isTextColorCoded = true;
    tv.grid.col(columnName)!.column!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
    tv.grid.col(columnName)!.column!.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"<=0.5":"#e87c79",">0.5":"#43b579"}`;
  }  
}

/**
 * Adds form-based results from the Admetox server to the specified DataFrame view.
 * The form properties and corresponding models are fetched from the 'properties' object.
 *
 * @param {DG.Column} smilesCol - The column containing SMILES data to be sent to the server.
 * @param {DG.DataFrame} viewTable - The DataFrame view to which the results will be added.
 */
export async function addForm(smilesCol: DG.Column, viewTable: DG.DataFrame) {
  const queryParams = Object.keys(properties)
    .flatMap(property => properties[property]['models'])
    .filter(obj => obj['skip'] !== true)
    .map(obj => obj['name'])
    .join(',');
  const csvString = await runAdmetox(DG.DataFrame.fromColumns([smilesCol]).toCsv(), queryParams);
  const table = csvString ? DG.DataFrame.fromCsv(csvString) : null;
  if (table)
    addResultColumns(table, viewTable);
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame): void {
  //substitute after demo
  const confidenceInterval = 0.5;
  if (table.columns.length > 0) {
    const modelNames: string[] = table.columns.names();
    const updatedModelNames: string[] = [];
    for (let i = 0; i < modelNames.length; ++i) {
      let column: DG.Column = table.columns.byName(modelNames[i]);
      column.name = viewTable.columns.getUnusedName(modelNames[i]);
      updatedModelNames[i] = column.name;
      column = column.convertTo("double");
      viewTable.columns.add(column);
    }
    viewTable.rows.removeWhere((row) => row.get('Y_Pgp-Inhibitor_probability') <= 0.5);
    addColorCoding(updatedModelNames);
  }
}

/**
 * Creates and returns an accordion with results for ADME/Tox models for a single molecule.
 *
 * @param {DG.SemanticValue<string>} smiles - The SMILES representation of the molecule.
 * @returns {DG.Accordion} An accordion control containing results for ADME/Tox models.
 */
export function getModelsSingle(smiles: DG.SemanticValue<string>): DG.Accordion {
  const acc = ui.accordion('ADME/Tox');
  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === 'ADME/Tox') 
      accPanes[i].append(ui.icons.help(() => {window.open('https://github.com/datagrok-ai/public/blob/1ef0f6c050754a432640301139f41fcc26e2b6c3/packages/ADMETox/README.md', '_blank')}));
  }
  const update = async (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties[modelName]['models']
    .filter((model: any) => model.skip !== false)
    .map((model: any) => model['name']);

    if (smiles.value === 'MALFORMED_INPUT_VALUE') {
      result.appendChild(ui.divText('The molecule is possibly malformed'));
    } else {
      result.appendChild(ui.loader());
      runAdmetox(
        `smiles
        ${smiles.value}`,
        queryParams.toString()
      ).catch((e: any) => {
        result.appendChild(ui.divText('Couldn\'t analyse properties'));
        console.log(e);
      }).then((csvString: any) => {
        ui.empty(result);
        const table = DG.DataFrame.fromCsv(csvString);
        const map: { [_: string]: any } = {};
        for (const model of queryParams)
          map[model] = Number(table.col(model)?.get(0)).toFixed(2);
  
          result.appendChild(ui.tableFromMap(map));
      }); 
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

/**
 * Opens a dialog for selecting models to apply on the specified SMILES column in the DataFrame view.
 *
 * @param {any} selected - An array containing the names of the initially selected models.
 * @param {DG.Column} smilesColumn - The column containing SMILES data.
 * @param {any} onOK - A function to be called when the 'OK' button is clicked, with the selected model names as an argument.
 * @returns {Promise<void>} A promise that resolves when the dialog is closed.
 */
async function openModelsDialog(selected: any, smilesColumn: DG.Column, onOK: any): Promise<void> {
  const tree = ui.tree();
  tree.root.style.maxHeight = '400px';
  tree.root.style.width = '200px';

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
  selectAll.style.marginLeft = '6px';
  selectAll.style.marginRight = '12px';
  const selectNone = ui.label('None', {classes: 'd4-link-label', onClick: () => checkAll(false)});

  const countLabel = ui.label('0 checked');
  countLabel.style.marginLeft = '24px';
  countLabel.style.display = 'inline-flex';

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
  }
  
  const df = grok.shell.tv.dataFrame;
  await grok.dapi.userDataStorage.postValue(_COLUMN_NAME_STORAGE, _KEY, smilesColumn.name);
  const molInput = ui.columnInput('Molecules', df, smilesColumn, async (col: DG.Column) => {
    smilesColumn = col;
    await grok.dapi.userDataStorage.postValue(_COLUMN_NAME_STORAGE, _KEY, smilesColumn.name);
    //@ts-ignore
  }, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE && col.getTag(DG.TAGS.UNITS) === DG.UNITS.Molecule.SMILES} as ColumnInputOptions);
  molInput.root.style.marginLeft = '-70px';
  const sliderInput = ui.sliderInput('Confidence interval', 0.6, 0, 1);
  const boolInput = ui.boolInput('Add columns', false);
  sliderInput.root.style.marginLeft = '-20px';
  boolInput.root.style.marginLeft = '-55px';
  const dlg = ui.dialog('ADME/Tox');
  dlg
    .add(molInput)
    .add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .add(sliderInput.root)
    .add(boolInput)
    .onOK(() => onOK(items.filter((i) => i.checked).map((i: any) => i.value['name'])))
    .show()
    .history(
      () => saveInputHistory(),
      (x) => loadInputHistory(x) 
    );
}

/**
 * Retrieves the selected models from the user data storage.
 *
 * @returns {Promise<any>} A promise that resolves to an array containing the selected models.
 */
async function getSelected() : Promise<any> {
  const str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  return selected;
}