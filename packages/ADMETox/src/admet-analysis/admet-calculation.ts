import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import { models, properties } from './const';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
let _COLOR = 'true';

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
  let response;
  try {
    //@ts-ignore
    response = await grok.dapi.docker.dockerContainers.request(admetDockerfile.id, path, params);
  } catch(e) {
    console.log(e);
    grok.shell.warning('Starting the ADME Docker Container that is required for this computation.<br>\
     Please try again in a few minutes');
  } finally {
    //@ts-ignore
    await grok.dapi.docker.dockerContainers.run(admetDockerfile.id);
  }
  return response;
}

export function addTooltip() {
  const tableView = grok.shell.tv;
  const { grid, dataFrame } = tableView;
  const between = (x: number, min: number, max: number) => x >= min && x <= max;

  grid.onCellTooltip((cell, x, y) => {
    const col = cell.tableColumn!.name;
    if (cell.isTableCell && Object.prototype.hasOwnProperty.call(models, col)) {
      const ranges = Object.keys(models[col]).map(parseFloat);
      const rowValue = dataFrame.get(cell.gridColumn.name, cell.gridRow);
      let val;
      for (const range of ranges) {
        if (between(rowValue, range, ranges[ranges.length - 1])) {
          val = models[col][range];
          break;
        }
        val = models[col][ranges[0]];
      }
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
    tv.grid.col(columnName)!.categoryColors = {
      '<0.5': 0xFFF1B6B4,
      '>0.5': 0xFFB4F1BC
    };
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

export async function addPredictions(smilesCol: DG.Column, viewTable: DG.DataFrame): Promise<void> {
  openModelsDialog(await getSelected(), async (selected: any) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    selected = await getSelected();
    const queryParams = selected.join(',');
    const pi = DG.TaskBarProgressIndicator.create('Evaluating predictions...');
    try {
      const [csvString] = await Promise.all([
      accessServer(DG.DataFrame.fromColumns([smilesCol]).toCsv(), queryParams),
      new Promise((resolve) => setTimeout(resolve, 1000))]);
      
      const table = processCsv(csvString);
      addResultColumns(table, viewTable);
      addTooltip();

      if (_COLOR === 'true') {
        let columnNames = table.columns.names();
        addColorCoding(columnNames);
      }
      
      pi.update(100, 'Evaluation complete');
    } catch (e: any) {
      console.error(e);
      grok.shell.warning(e.toString());
      pi.update(100, 'Evaluation failed');
    } finally {
      pi.close();
    }
  })
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame): void {
  if (table.columns.length > 0) {
    const modelNames: string[] = table.columns.names()
    for (let i = 0; i < modelNames.length; ++i) {
      let column: DG.Column = table.columns.byName(modelNames[i]);
      column.name = viewTable.columns.getUnusedName(modelNames[i]);
      column = column.convertTo("double");
      viewTable.columns.add(column);
    }
  }
}

function processCsv(csvString: string | null): DG.DataFrame {
  csvString = csvString!.replaceAll('"', '');
  const table = DG.DataFrame.fromCsv(csvString);
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

export function getModelsSingle(smiles: string): DG.Accordion {
  const acc = ui.accordion('ADME/Tox');
  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === 'ADME/Tox') 
      accPanes[i].append(ui.icons.help(() => {window.open('https://github.com/datagrok-ai/public/blob/1ef0f6c050754a432640301139f41fcc26e2b6c3/packages/ADMETox/README.md', '_blank')}));
  }
  const update = (result: HTMLDivElement, modelName: string) => {
    const queryParams = properties[modelName]['models'].map((model: any) => model['name']);
    result.appendChild(ui.loader());
    accessServer(
      `smiles
      ${smiles}`,
      queryParams.toString()
    ).catch((e: any) => {
      console.log(e);
      grok.shell.warning(e.toString());
    }).then((csvString: any) => {
      removeChildren(result);
      const table = processCsv(csvString);
      const map: { [_: string]: any } = {};
      for (const model of queryParams)
        map[model] = Number(table.col(model)?.get(0)).toFixed(2);

        result.appendChild(ui.tableFromMap(map));
    });
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

function openModelsDialog(selected: any, onOK: any): void {
  /*_package.files.readAsText('properties.json').then((res) => {
    properties = JSON.parse(res);
  })*/

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
      const item = group.item(property['name'], property);
      item.enableCheckBox(selected.includes(property['name']));
      items.push(item);

      item.checkBox!.onchange = (_e) => {
        countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
        if (item.checked) 
          selectedModels[item.text] = groupName;
      };
    }

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

  let bIColor = ui.boolInput('Color Coding', true);
  tree.root.appendChild(bIColor.root);

  bIColor.onChanged(async () => {
    _COLOR = bIColor.value!.toString();
  });
  

  let dlg = ui.dialog('ADME/Tox');
  dlg.add(ui.divH([selectAll, selectNone, countLabel]))
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
