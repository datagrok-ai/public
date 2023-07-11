import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {_package} from '../package';
import {addCopyIcon} from '../utils/ui-utils';
import {MESSAGE_MALFORMED} from '../constants';

const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';
const MOLECULE_NAME = 'smiles';
let descriptors: {[_: string]: any};

// gets a widget with calculated descriptors for a cell with molecule
export function getDescriptorsSingle(smiles: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  smiles = _convertMolNotation(smiles, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  if (smiles === MESSAGE_MALFORMED)
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  const molecule = DG.chem.isMolBlock(smiles) ? `\"${smiles}\"` : smiles;
  const widget = new DG.Widget(ui.div());
  const result = ui.div();
  const selectButton = ui.bigButton('SELECT', async () => {
    openDescriptorsDialog(await getSelected(), async (selected: string []) => {
      await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      update();
    });
  });
  selectButton.style.marginTop = '20px';

  const update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    getSelected().then((selected) => {
      grok.chem.descriptors(DG.DataFrame.fromCsv(`${MOLECULE_NAME}\n${molecule}`), MOLECULE_NAME, selected)
        .then((table: any) => {
          removeChildren(result);
          const map: { [_: string]: any } = {};
          for (const descriptor of selected)
            map[descriptor] = table.col(descriptor).get(0);
          result.appendChild(ui.tableFromMap(map));
        });
    });
  };

  addCopyIcon(result, 'Descriptors');
  widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();

  return widget;
}

//description: Open descriptors selection dialog
async function openDescriptorsDialog(selected: any, onOK: any): Promise<void> {
  descriptors = await grok.chem.descriptorsTree();
  const tree = ui.tree();
  tree.root.style.maxHeight = '400px';

  const groups: { [_: string]: any } = {};
  const items: DG.TreeViewNode[] = [];
  const selectedDescriptors: { [_: string]: string } = {};

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

  const keys = Object.keys(descriptors);
  for (const groupName of keys) {
    const group = tree.group(groupName, null, false);
    group.enableCheckBox();
    groups[groupName] = group;

    group.checkBox!.onchange = (_e) => {
      countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
      if (group.checked)
        selectedDescriptors[group.text] = group.text;
      group.items.filter((i) => {
        if (i.checked)
          selectedDescriptors[i.text] = group.text;
      });
    };

    for (const descriptor of descriptors[groupName]['descriptors']) {
      const item = group.item(descriptor['name'], descriptor);
      item.enableCheckBox(selected.includes(descriptor['name']));
      items.push(item);

      item.checkBox!.onchange = (_e) => {
        countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
        if (item.checked)
          selectedDescriptors[item.text] = groupName;
      };
    }

    checkAll(false);
  }

  const saveInputHistory = (): any => {
    const resultHistory: { [_: string]: any } = {};
    const descriptorNames = Object.keys(selectedDescriptors);
    for (const descriptorName of descriptorNames)
      resultHistory[descriptorName] = selectedDescriptors[descriptorName];
    return resultHistory;
  };

  const loadInputHistory = (history: any): void => {
    checkAll(false);
    const keys: string[] = Object.keys(history);
    for (const key of keys) {
      groups[history[key]].items.filter(function(i: any) {
        if (i.text === key)
          i.checked = true;
      });
      if (key === history[key])
        groups[history[key]].checked = true;
    }
    countLabel.textContent = `${keys.length} checked`;
  };

  ui.dialog('Chem Descriptors')
    .add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .onOK(() => onOK(items.filter((i) => i.checked).map((i: any) => i.value['name'])))
    .show()
    .history(
      () => saveInputHistory(),
      (x) => loadInputHistory(x),
    );
}

//description: Get selected descriptors
async function getSelected() : Promise<string []> {
  const str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  if (selected.length === 0) {
    selected = descriptors['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Removes all children from node
function removeChildren(node: any): void {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}
