import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import { getDescriptorsTree, getDescriptorsPy } from '../scripts-api';

const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';

/**
 * Adds descriptors to table
 * @export
 * @param {DG.Column} smiles column with smiles.
 * @param {DG.DataFrame} viewTable current view table.
 */
export async function addDescriptors(smiles: DG.Column, viewTable: DG.DataFrame): Promise<void> {

  openDescriptorsDialog(await getSelected(), async (selected: any) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    getSelected().then(selected => {
      let pi = DG.TaskBarProgressIndicator.create('Calculating descriptors');
      //grok.chem.descriptors(DG.DataFrame.fromColumns([smiles]), smiles.name, selected).then((table: any) => {
      getDescriptorsPy(smiles.name, DG.DataFrame.fromColumns([smiles]), 'selected', DG.DataFrame.fromColumns([DG.Column.fromList('string', 'selected', selected)])
      ).then((table: any) => {
        addResultColumns(table, viewTable);
      });
      pi.close();
    });
  });
}

/**
 * calculates descriptors for single entry
 * @export
 * @param {DG.Column} smiles column with smiles.
 */
export function getDescriptorsSingle(smiles: string) {
  let widget = new DG.Widget(ui.div());
  let result = ui.div();
  let selectButton = ui.bigButton('SELECT', async () => {
    openDescriptorsDialog(await getSelected(), async (selected: any) => {
      await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      update();
    });
  });
  selectButton.style.marginTop = '20px';

  let update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    getSelected().then(selected => {
      getDescriptorsPy('smiles', DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'selected', DG.DataFrame.fromColumns([DG.Column.fromList('string', 'selected', selected)])
      ).then((table: any) => {
      //grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then((table: any) => {
        removeChildren(result);
        let map: { [_: string]: any } = {};
        for (let descriptor of selected)
          map[descriptor] = table.col(descriptor).get(0);
        result.appendChild(ui.tableFromMap(map));
      });
    });
  }

  widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();

  return widget;
}

/**
 * demo application for descriptor calculation
 * @export
 */
export function getDescriptorsApp() {
  let defaultSmiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
  let sketcherValue = defaultSmiles;

  let windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  let table = DG.DataFrame.create();
  table.name = 'Descriptors';
  let view = grok.shell.addTableView(table);

  let dsDiv = ui.divV([], 'grok-prop-panel');
  dsDiv.appendChild(getDescriptorsSingle(defaultSmiles).root);

  let sketcher = grok.chem.sketcher((smiles: string, molfile: string) => {
    sketcherValue = smiles;
    removeChildren(dsDiv);
    dsDiv.appendChild(getDescriptorsSingle(smiles).root);
  }, defaultSmiles);
  let addButton = ui.bigButton('ADD', async () => {
    getSelected().then(selected => {
      getDescriptorsPy('smiles', DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'selected', DG.DataFrame.fromColumns([DG.Column.fromList('string', 'selected', selected)])
      ).then(t => {
      //grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
        let columnNames = table.columns.names();
        if ((table.columns.length !== selected.length + 1) || selected.some((s: any) => !columnNames.includes(s))) {
          table = DG.DataFrame.create();
          table.name = 'Descriptors';
          view.dataFrame = table;
          for (let col of t.columns.toList())
            table.columns.addNew(col.name, col.type);
        }
        table.rows.addNew(t.columns.toList().map((c: any) => c.get(0)));
      });
    });
  });
  addButton.style.marginTop = '12px';
  let skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

  let skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
  view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

  grok.events.onViewRemoved.subscribe((v: any) => {
    if (v.name === view.name) {
      windows.showToolbox = true;
      windows.showHelp = true;
      windows.showProperties = true;
    }
  });
}

//description: Open descriptors selection dialog
function openDescriptorsDialog(selected: any, onOK: any) {
  //grok.chem.descriptorsTree().then((descriptors: { [_: string]: any }) => {
  getDescriptorsTree().then((descriptors: { [_: string]: any }) => {
    let tree = ui.tree();
    tree.root.style.maxHeight = '400px';

    let groups: { [_: string]: any } = {};
    let items: DG.TreeViewNode[] = [];

    const checkAll = (val: boolean) => {
      for (let g in groups) groups[g].checked = val;
      for (let i of items) i.checked = val;
    };

    let selectAll = ui.label('All', { classes: 'd4-link-label', onClick: () => checkAll(true) });
    selectAll.style.marginLeft = '6px'
    selectAll.style.marginRight = '12px'
    let selectNone = ui.label('None', { classes: 'd4-link-label', onClick: () => checkAll(false) });

    let countLabel = ui.label('0 checked')
    countLabel.style.marginLeft = '24px';
    countLabel.style.display = 'inline-flex';

    let keys = Object.keys(descriptors);
    for (let groupName in keys) {
      let group = tree.group(keys[groupName], null, false);
      group.enableCheckBox();
      groups[keys[groupName]] = group;

      group.checkBox!.onchange = (e) => {
        countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
      };

      for (let descriptor of descriptors[keys[groupName]]['descriptors']) {
        let item = group.item(descriptor['name'], descriptor);
        item.enableCheckBox(selected.includes(descriptor['name']));
        items.push(item);

        item.checkBox!.onchange = (e) => {
          countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
        };
      }

      checkAll(false);
    }

    ui.dialog('Chem Descriptors')
      .add(ui.divH([selectAll, selectNone, countLabel]))
      .add(tree.root)
      .onOK(() => onOK(items.filter(i => i.checked).map((i: any) => i.value['name'])))
      .show();
  });
}

//description: Get selected descriptors
async function getSelected() {
  let str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  if (selected.length === 0) {

    //selected = (await grok.chem.descriptorsTree() as any)['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    selected = (await getDescriptorsTree() as any)['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Removes all children from node
function removeChildren(node: any) {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}

//description: add columns into table.
function addResultColumns(table: DG.DataFrame, viewTable: DG.DataFrame): void {

  if (table.columns.length > 0) {
    let descriptors: string[] = table.columns.names();

    for (let i = 1; i < descriptors.length; i++) {
      let column: DG.Column = table.columns.byName(descriptors[i]);
      column.name = getName(column.name, viewTable.columns.names());
      viewTable.columns.add(column);
    }
  }
}

function getName(initialName: string, existingNames: string[]) {
  if (!existingNames.includes(initialName)) {
    return initialName;
  } else {
    let counter: number = 1;
    let newName: string = (' ' + initialName + '_' + counter).slice(1);
    while (existingNames.includes(newName)) {
      counter++;
      newName = (' ' + initialName + '_' + counter).slice(1);
    }

    return newName;
  }
}
