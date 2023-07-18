import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { models, properties } from './const';
import { delay } from '@datagrok-libraries/utils/src/test';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
const _COLUMN_NAME_STORAGE = 'column_names';

export async function accessServer(csvString: string, queryParams: string) {
  //@ts-ignore
  const admetDockerfile = await grok.dapi.docker.dockerContainers.filter('admetox').first();
  const params: RequestInit = {
    method: 'POST',
    headers: {
      'Accept': 'text/csv',
      'Content-type': 'text/csv'
    },
    body: csvString
  };
  const path = `/smiles/df_upload/?models=${queryParams}`;
  let response = await grok.dapi.docker.dockerContainers.request(admetDockerfile.id, path, params);
  return response;
}

export function addTooltip() {
  const tableView = grok.shell.tv;

  tableView.grid.onCellTooltip((cell, x, y) => {
    const col = cell.tableColumn?.name;

    if (cell.isTableCell && col && models[col]) {
      const keys = Object.keys(models[col]);
      const rowValue = tableView.dataFrame.get(cell.gridColumn.name, cell.gridRow);

      let val = '';
      keys.some((range, i) => {
        const nextRange = keys[i + 1];
        const isLastRange = i === keys.length - 1;

        if (nextRange && rowValue >= +range && rowValue <= +nextRange) {
          val = models[col][range];
          return true;
        } else if (isLastRange) {
          val = models[col][range];
          return true;
        }
      });

      ui.tooltip.show(ui.divV([ui.div(val)]), x, y);
      return true;
    }
  });
}

export function addColorCoding(columnNames: string[]) {
  const tv = grok.shell.tv;
  for (const columnName of columnNames) {
    //@ts-ignore
    tv.grid.col(columnName)!.isTextColorCoded = true;
    tv.grid.col(columnName)!.column!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
    tv.grid.col(columnName)!.column!.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"<=0.5":"#e87c79",">0.5":"#43b579"}`;
  }  
}

export async function addForm(smilesCol: DG.Column, viewTable: DG.DataFrame) {
  const queryParams = Object.keys(properties)
    .flatMap(property => properties[property]['models'])
    .filter(obj => obj['skip'] !== true)
    .map(obj => obj['name'])
    .join(',');
  let csvString = await accessServer(DG.DataFrame.fromColumns([smilesCol]).toCsv(), queryParams);
  const table = processCsv(csvString);
  addResultColumns(table, viewTable);
}

async function getPredictions(viewTable: DG.DataFrame, column?: DG.Column): Promise<DG.DataFrame> {
  const selected = await getSelected();

  return new Promise<DG.DataFrame>((resolve, reject) => {
    openModelsDialog(selected, async (selected: any) => {
      await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      if (selected.length === 0) {
        grok.shell.warning('No models have been selected!');
        reject('No models selected');
        return;
      }
      const queryParams = selected.join(',');

      const colName = await grok.dapi.userDataStorage.getValue(_COLUMN_NAME_STORAGE, _KEY);
      const smilesCol = viewTable.columns.byName(colName);
      const malformedIndexes = await grok.functions.call('Chem:_getMolSafe', {molecules: smilesCol});
      const smilesColFiltered = DG.Column.fromStrings(smilesCol.name, Array.from(smilesCol.values()).filter((_, index) => !malformedIndexes.includes(index)));

      if (smilesCol.length > 10000) {
        const dialog = ui.dialog({ title: 'Proceed with computations' });
        dialog
          .add(
            ui.divText(
              `Performing computations could potentially consume a significant amount of time.
               Do you want to continue?`
            )
          )
          .addButton('YES', async () => {
            dialog.close();
            const result = await addColumnsAndProcessInBatches(viewTable, selected, smilesColFiltered, queryParams, malformedIndexes);
            resolve(result);
          })
          .show();
      } else {
        const result = await addColumnsAndProcessInBatches(viewTable, selected, smilesColFiltered, queryParams, malformedIndexes);
        resolve(result);
      }
    });
  });
}

export async function addPredictions(viewTable: DG.DataFrame) {
  const df = await getPredictions(viewTable);
  addResultColumns(df, viewTable);
}


async function addColumnsAndProcessInBatches(viewTable: DG.DataFrame, selected: any[], smilesCol: DG.Column, queryParams: string, malformedIndexes: number[]): Promise<DG.DataFrame> {
  return await processColumnInBatches(smilesCol, viewTable, 100, queryParams, selected, malformedIndexes);
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame, malformedIndexes?: any[]): void {
  if (table.columns.length > 0) {
    if (malformedIndexes)
      malformedIndexes.map((index) => table.rows.insertAt(index));
    const modelNames: string[] = table.columns.names()
    for (let i = 0; i < modelNames.length; ++i) {
      let column: DG.Column = table.columns.byName(modelNames[i]);
      column.name = viewTable.columns.getUnusedName(modelNames[i]);
      column = column.convertTo("double");
      viewTable.columns.add(column);
    }
    addColorCoding(modelNames);
    addTooltip();
  }
}

function processCsv(csvString: string | null | undefined): DG.DataFrame {
  csvString = csvString!.replaceAll('"', '');
  const table = DG.DataFrame.fromCsv(csvString);
  table.rows.removeAt(table.rowCount - 1);
  const removeRow = Array.from(table.row(table.rowCount - 1).cells).some((cell: DG.Cell) => cell.value == '');
  if (removeRow)
    table.rows.removeAt(table.rowCount - 1);
  const modelNames: string[] = [];
  const prevColNames = table.columns.names();
  for (let i = 0; i < prevColNames.length; ++i) {
    modelNames[i] = table.get(prevColNames[i], 0);
  }
  table.rows.removeAt(0);
  for (let i = 0; i < prevColNames.length; ++i) {
    const column: DG.Column = table.columns.byName(prevColNames[i]);
    column.name = column.dataFrame.columns.getUnusedName(modelNames[i]);
  }
  return table;
}

export function getModelsSingle(smiles: DG.SemanticValue<string>): DG.Accordion {
  const acc = ui.accordion('ADME/Tox');
  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === 'ADME/Tox') 
      accPanes[i].append(ui.icons.help(() => {window.open('https://github.com/datagrok-ai/public/blob/1ef0f6c050754a432640301139f41fcc26e2b6c3/packages/ADMETox/README.md', '_blank')}));
  }
  const update = async (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties[modelName]['models'].map((model: any) => model['name']);
    if (smiles.value === 'MALFORMED_INPUT_VALUE') {
      result.appendChild(ui.divText('The molecule is possibly malformed'));
    } else {
      result.appendChild(ui.loader());
      accessServer(
        `smiles
        ${smiles.value}`,
        queryParams.toString()
      ).catch((e: any) => {
        console.log(e);
        //grok.shell.warning(e.toString());
      }).then((csvString: any) => {
        removeChildren(result);
        const table = processCsv(csvString);
        const map: { [_: string]: any } = {};
        for (const model of queryParams)
          map[model] = Number(table.col(model)?.get(0)).toFixed(2);
  
          result.appendChild(ui.tableFromMap(map));
      }); 
    }
  };


  for (const property of Object.keys(properties)) {
    const result = ui.div();
    acc.addPane(property, () => {
      update(result, property);
      return result;
    }, false);
  }

  return acc;
}

async function openModelsDialog(selected: any, onOK: any): Promise<void> {
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
  let column = DG.Utils.firstOrNull(df.columns.byTags({'units': 'smiles'}));
  await grok.dapi.userDataStorage.postValue(_COLUMN_NAME_STORAGE, _KEY, column!.name);
  const molInput = ui.columnInput('Molecules', df, column, 
  async (col: DG.Column) => {
    column = col;
    await grok.dapi.userDataStorage.postValue(_COLUMN_NAME_STORAGE, _KEY, column.name);
    //@ts-ignore
  }, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE && col.getTag(DG.TAGS.UNITS) === DG.UNITS.Molecule.SMILES} as ColumnInputOptions);
  molInput.root.style.marginLeft = '-70px';
  let dlg = ui.dialog('ADME/Tox');
  dlg
    .add(molInput)
    .add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .onOK(() => onOK(items.filter((i) => i.checked).map((i: any) => i.value['name'])))
    .show()
    .history(
      () => saveInputHistory(),
      (x) => loadInputHistory(x) 
    );
}

async function getSelected() : Promise<any> {
  const str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  return selected;
}

//description: Removes all children from node
function removeChildren(node: any): void {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}

let resultDf: DG.DataFrame;
let entered: boolean = false;
async function processBatch(batch: any, queryParams: string, viewTable: DG.DataFrame, startIndex: number, batchSize: number, modelNames: any[], malformedIndexes: any[]) {
  const [csvString] = await Promise.all([
    accessServer(DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', batch)]).toCsv(), queryParams),
    new Promise((resolve) => setTimeout(resolve, 1000))]);
  if (!entered) {
    resultDf = processCsv(csvString);
    entered = true;
  } else {
    const table = processCsv(csvString);
    resultDf.append(table, true);
  }
}

class Semaphore {
  concurrency: number;
  current: number;
  queue: (() => void)[];

  constructor(concurrency: number) {
    this.concurrency = concurrency;
    this.current = 0;
    this.queue = [];
  }

  async acquire() {
    if (this.current < this.concurrency) {
      this.current++;
      return Promise.resolve();
    } else {
      return new Promise<void>(resolve => {
        this.queue.push(resolve);
      });
    }
  }

  release() {
    this.current--;
    if (this.queue.length > 0) {
      const resolve = this.queue.shift();
      resolve?.();
    }
  }
}

async function processColumnInBatches(column: DG.Column, viewTable: DG.DataFrame, batchSize = 100, queryParams: string, modelNames: any[], malformedIndexes: any[]): Promise<DG.DataFrame> {
  //@ts-ignore
  const progressIndicator = DG.TaskBarProgressIndicator.create('Evaluating predictions...', {cancelable: true});
  //@ts-ignore
  progressIndicator.onCanceled.subscribe(() => {
    progressIndicator.close();
  });
  const semaphore = new Semaphore(10);

  const totalBatches = Math.ceil(column.length / batchSize);
  let processedBatches = 0;

  for (let batchIndex = 0; batchIndex < totalBatches; batchIndex++) {
    await semaphore.acquire();

    const start = batchIndex * batchSize;
    const end = start + batchSize;
    const batch = Array.from(column.values()).slice(start, end);

    processBatch(batch, queryParams, viewTable, start, batchSize, modelNames, malformedIndexes)
      .then(() => {
        processedBatches++;
        const percent = (processedBatches / totalBatches) * 100;
        progressIndicator.update(percent, `${percent.toFixed(2)}% is evaluated...`);
        semaphore.release();
      })
      .catch(error => {
        console.error('Error processing batch:', error);
        semaphore.release();
      });
  }

  while (processedBatches < totalBatches) {
    await delay(100);
  }

  progressIndicator.close();
  return resultDf;
}