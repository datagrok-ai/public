import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { properties, models } from './const';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';
let _COLOR = 'true';

export async function accessServer(csvString: string, queryParams: string) {
  const admetDockerfile = await grok.dapi.docker.dockerContainers.filter('admetox').first();
  const params: RequestInit = {
    method: 'POST',
    headers: {
      'Accept': 'text/csv',
      'Content-type': 'text/csv'
    },
    //body: DG.DataFrame.fromColumns([smilesCol]).toCsv()
    body: csvString
  };
  const path = `/smiles/df_upload/?models=${queryParams}`;
  const response = await grok.dapi.docker.dockerContainers.request(admetDockerfile.id, path, params);
  return response;
}

export function addTooltip() {
  const tableView = grok.shell.tv;
  const between = (x: number, min: number, max: number) => {
    return x >= min && x <= max;
  }
  tableView.grid.onCellTooltip(function (cell, x, y) {
    const col = cell.tableColumn!.name;
    const rangeNumbers = [];
    if (cell.isTableCell && Object.keys(models).includes(col)) {
      const keys = Object.keys(models[col])
      for (let i = 0; i < keys.length; ++i) {
        rangeNumbers[i] = keys[i];
      }
      const rowValue = tableView.dataFrame.get(cell.gridColumn.name, cell.gridRow);
      let val = '';
      //let val: any = _.inRange(rowValue) <= 0.5 ? Object.values(models[col])[0] : Object.values(models[col])[1];
      for (let i = 0; i <= rangeNumbers.length; ++i) {
        if (i != rangeNumbers.length - 1 && between(rowValue, parseFloat(rangeNumbers[i]), parseFloat(rangeNumbers[i + 1]))) {
          val = models[col][rangeNumbers[i]];
          break;
        } else if (i === rangeNumbers.length - 1){
          val = models[col][rangeNumbers[i]];
          break;
        } else {
          val = models[col][rangeNumbers[i + 1]];
        }
      }
      ui.tooltip.show(ui.divV([
        ui.div(val)
      ]), x, y);
      return true;
    }
  });
}

export function addColorCoding(columnNames: string[], viewTable: DG.DataFrame) {
  for (const columnName of columnNames) {
    viewTable.col(columnName)!.tags['.color-coding-type'] = 'Conditional';
    viewTable.col(columnName)!.tags['.color-coding-conditional'] = `{"0-0.5":"#f1b6b4","0.5-1":"#b4f1bc"}`;
  }   
}

export async function addForm(smilesCol: DG.Column, viewTable: DG.DataFrame) {
  let queryParams = 'Pgp-Inhibitor,Pgp-Substrate,HIA,F(20%),F(30%),Ames,SkinSen,BBB,CYP1A2-Inhibitor,CYP1A2-Substrate,CYP3A4-Inhibitor,CYP3A4-Substrate,CYP2C19-Inhibitor,CYP2C19-Substrate,CYP2C9-Inhibitor,CYP2C9-Substrate,CYP2D6-Inhibitor,CYP2D6-Substrate';
  let csvString = await accessServer(DG.DataFrame.fromColumns([smilesCol]).toCsv(), queryParams);
  const table = processCsv(csvString);
  addResultColumns(table, viewTable);
  addColorCoding(table.columns.names(), viewTable);
}

export async function addPredictions(smilesCol: DG.Column, viewTable: DG.DataFrame): Promise<void> {
  openModelsDialog(await getSelected(), async (selected: any) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    selected = await getSelected();
    const queryParams = selected.join(',');
    const pi = DG.TaskBarProgressIndicator.create('Evaluating predictions...');
    let csvString: string | null = '';
    try {
      csvString = await accessServer(DG.DataFrame.fromColumns([smilesCol]).toCsv(), queryParams);
    } catch (e: any) {
      console.error(e);
      grok.shell.warning(e.toString());
      pi.update(100, 'Evaluation failed...');
      pi.close();
    }
    const table = processCsv(csvString);
    addResultColumns(table, viewTable);
    addTooltip();
    console.log(_COLOR);
    if (_COLOR === 'true') {
      let columnNames = table.columns.names();
      addColorCoding(columnNames, viewTable);
    }
    pi.close();
  });
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame): void {
  if (table.columns.length > 0) {
    const modelNames: string[] = table.columns.names()
    for (let i = 0; i < modelNames.length; ++i) {
      let column: DG.Column = table.columns.byName(modelNames[i]);
      column.name = viewTable.columns.getUnusedName(modelNames[i]);
      column.setTag(column.name, models[column.name]);
      column = column.convertTo("double");
      viewTable.columns.add(column);
    }
  }
}

function processCsv(csvString: string | null): DG.DataFrame {
  csvString = csvString!.replaceAll('"', '');
  let table = DG.DataFrame.fromCsv(csvString);
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
  acc.root.appendChild(ui.divV([
    ui.link('Detailed info', () => window.open('README.md link', '_blank'))]));
  let accHeader = document.getElementsByClassName('d4-accordion-pane-header')[17] as HTMLElement;
  accHeader.append(ui.icons.help(() => {window.open('README.md link', '_blank')}));
  const update = (result: HTMLDivElement, modelName: string) => {
    let queryParams: string[] = [];
    let model = properties[modelName]['models']
    for (let i = 0; i < model.length; ++i) 
      queryParams[i] = model[i]['name'];
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
        map[model] = table.col(model)?.get(0);

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

  /*widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();*/

  return acc;
}

function openModelsDialog(selected: any, onOK: any): void {
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
