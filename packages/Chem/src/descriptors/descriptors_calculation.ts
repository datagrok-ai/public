import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';

export async function getDescriptors(smiles: DG.Column): Promise<void> {



  openDescriptorsDialog(await getSelected(), async (selected: any) => {
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    getSelected().then(selected => {
      grok.chem.descriptors(DG.DataFrame.fromColumns([smiles]), smiles.name, selected).then((table: any) => {

        let map: { [_: string]: any } = {};
        for (let descriptor of selected)
          map[descriptor] = table.col(descriptor).get(0);
 
      });
    });
  });


  // let t2 = DG.DataFrame.fromColumns([
  //   DG.Column.fromList('int', 'smi', smiles.toList()),
  // ]);


}

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
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then((table: any) => {
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
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
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
  grok.chem.descriptorsTree().then((descriptors: { [_: string]: any }) => {
    let tree = ui.tree();
    tree.root.style.maxHeight = '400px';

    let groups: { [_: string]: any } = {};
    let items: DG.TreeViewNode[] = [];

    for (let groupName in descriptors) {
      let group = tree.group(groupName, null, false);
      group.enableCheckBox();
      groups[groupName] = group;

      for (let descriptor of descriptors[groupName]['descriptors']) {
        let item = group.item(descriptor['name'], descriptor);
        item.enableCheckBox(selected.includes(descriptor['name']));
        items.push(item);
      }
    }

    let clear = ui.button('NONE', () => {
      for (let g in groups) groups[g].checked = false;
      for (let i of items) i.checked = false;
    });

    ui.dialog('Chem Descriptors')
      .add(clear)
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
    selected = (await grok.chem.descriptorsTree() as any)['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Removes all children from node
function removeChildren(node: any) {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}