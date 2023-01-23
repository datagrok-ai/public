import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { properties } from './const';

const _STORAGE_NAME = 'admet_models';
const _KEY = 'selected';

async function accessServer(csvString: string, queryParams: string) {
  const dockerId = (await grok.dapi.dockerfiles.filter('admetox').first()).id;
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
  const response = await grok.dapi.dockerfiles.request(dockerId, path, params);
  return response;
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
    } catch (e) {
      console.error(e);
      grok.shell.warning('Dataset contains malformed data!');
      pi.update(100, 'Evaluation failed...');
      pi.close();
    }
    const table = processCsv(csvString);
    addResultColumns(table, viewTable);
    pi.close();
  });
}

function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame): void {
  if (table.columns.length > 0) {
    const models: string[] = table.columns.names()
    for (let i = 0; i < models.length; ++i) {
      const column: DG.Column = table.columns.byName(models[i]);
      column.name = viewTable.columns.getUnusedName(models[i]);
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

export function getModelsSingle(smiles: string): DG.Widget {
  const widget = new DG.Widget(ui.div());
  const result = ui.div();
  const selectButton = ui.bigButton('SELECT', async () => {
    openModelsDialog(await getSelected(), async (selected: any) => {
      await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      update();    
    });
  });
  selectButton.style.marginTop = '20px';
  
  const update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    getSelected().then((selected) => {
      const queryParams = selected.join(',');
      accessServer(
       `smiles
        ${smiles}`,
        queryParams
      ).catch((e) => {
        console.log(e);
        grok.shell.warning('Dataset contains malformed data!');
      }).then((csvString: any) => {
        removeChildren(result);
        const table = processCsv(csvString);
        const map: { [_: string]: any } = {};
        for (const model of selected)
          map[model] = table.col(model)?.get(0);

        result.appendChild(ui.tableFromMap(map));
      });
    });
  };

  widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();

  return widget;
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

  let dlg = ui.dialog('ADME/Tox');
  dlg.add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .addButton('Form', async () => {
      const pi = DG.TaskBarProgressIndicator.create('Creating a form...');
      dlg.close();
      let queryParams = 'Pgp-Inhibitor,Pgp-Substrate,HIA,F(20%),F(30%),BBB,CYP1A2-Inhibitor,CYP1A2-Substrate,CYP3A4-Inhibitor,CYP3A4-Substrate,CYP2C19-Inhibitor,CYP2C19-Substrate,CYP2C9-Inhibitor,CYP2C9-Substrate,CYP2D6-Inhibitor,CYP2D6-Substrate,Ames,SkinSen';
      const df = grok.shell.tv.grid.dataFrame;
      const col = df.columns.bySemType(DG.SEMTYPE.MOLECULE)
      if (col) {
        const csvString = await accessServer(DG.DataFrame.fromColumns([col]).toCsv(), queryParams);
        const table = processCsv(csvString);
        addResultColumns(table, df);
        pi.close();
        grok.shell.tv.addViewer(DG.Viewer.fromType('Form', df));
      }
    })
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
