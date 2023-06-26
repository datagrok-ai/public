import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {ChemPropNames, IDescriptorTree} from './types';
import {ChemPropertyGroupMap} from './consts';


export async function chemDescriptorsDialog(onOk: (selected: {[_: string]: string[]}) => void, onCancel: () => void) {
  // select all, none buttons and count label
  const selectAll = ui.label('All', {classes: 'd4-link-label', onClick: () => toggleAll(true)});
  selectAll.style.marginLeft = '6px';
  selectAll.style.marginRight = '12px';
  const selectNone = ui.label('None', {classes: 'd4-link-label', onClick: () => toggleAll(false)});
  const countLabel = ui.label('0 checked');
  countLabel.style.marginLeft = '24px';
  countLabel.style.display = 'inline-flex';

  // tree groups
  const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
  const tree = ui.tree();
  tree.root.style.maxHeight = '400px';
  tree.root.style.width = '300px';

  const descriptorGroups: { [_: string]: DG.TreeViewGroup } = {};
  const descriptorItems: DG.TreeViewNode[] = [];
  const toxicityItems: DG.TreeViewNode[] = [];
  const alertsItems: DG.TreeViewNode[] = [];
  const propertiesItems: DG.TreeViewNode[] = [];
  const allItemsMap: { [_ in ChemPropNames]: DG.TreeViewNode[] } = {
    'Descriptors': descriptorItems,
    'Toxicity Risks': toxicityItems,
    'Structural Alerts': alertsItems,
    'Chemical Properties': propertiesItems,
  };

  function onSelectionChange() {
    countLabel.textContent = `${
      Object.values(allItemsMap)
        .map((items) => items.filter((i) => i.checked).length).reduce((a, b) => a + b, 0)
    } checked`;
  }
  function createTreeGroup(name: string, treeNode: DG.TreeViewGroup): DG.TreeViewGroup {
    const res = treeNode.group(name, null, false);
    res.enableCheckBox(false);
    res.checkBox!.onchange = (_e) => onSelectionChange();
    return res;
  };
  const toxicityGroup = createTreeGroup('Toxicity Risks', tree);
  const alertsGroup = createTreeGroup('Structural Alerts', tree);
  const propertiesGroup = createTreeGroup('Chemical Properties', tree);
  const descriptorsGroup = createTreeGroup('Descripors', tree);

  const allGroupsMap: { [_ in ChemPropNames]: DG.TreeViewGroup } = {
    'Descriptors': descriptorsGroup,
    'Toxicity Risks': toxicityGroup,
    'Structural Alerts': alertsGroup,
    'Chemical Properties': propertiesGroup,
  };
  const toggleAll = (val: boolean) => {
    for (const g of Object.values(descriptorGroups))
      g.checked = val;
    Object.values(allGroupsMap).forEach((g) => g.checked = val);
    Object.values(allItemsMap).forEach((items) => {
      for (const i of items)
        i.checked = val;
    });
    onSelectionChange();
  };
  // processing chem descriptprs tree
  const keys = Object.keys(descriptorsTree);
  for (const groupName of keys) {
    const group = createTreeGroup(groupName, descriptorsGroup);
    descriptorGroups[groupName] = group;

    for (const descriptor of descriptorsTree[groupName].descriptors) {
      const item = group.item(descriptor.name, descriptor);
      descriptorItems.push(item);
      item.enableCheckBox(false);
        item.checkBox!.onchange = (_e) => onSelectionChange();
    }
    toggleAll(false);
  };

  //adding toxicity risks, alerts and properties to the tree

  Object.values(ChemPropertyGroupMap).forEach((group) => {
    const groupName = group.name;
    const curGroup = allGroupsMap[groupName];
    group.values.forEach((prop) => {
      const item = curGroup.item(prop.name, prop);
      item.enableCheckBox(false);
      item.checkBox!.onchange = (_e) => onSelectionChange();
      allItemsMap[groupName].push(item);
    });
  });

  function onOkProxy() {
    const res: {[_: string]: string[]} = {};
    Object.keys(allItemsMap).forEach((groupName) => {
      res[groupName] = allItemsMap[groupName as ChemPropNames].filter((i) => i.checked)
        .map((i) => i.value[groupName === 'Descriptors' ? 'name' : 'propertyName']);
    });
    onOk(res);
  }

  // enable toxicity risks and structural alerts by default
  toxicityGroup.checked = true;
  alertsGroup.checked = true;
  toxicityItems.forEach((i) => i.checked = true);
  alertsItems.forEach((i) => i.checked = true);
  onSelectionChange();

  ui.dialog('Select descriptors')
    .add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .onOK(() => onOkProxy())
    .onCancel(onCancel)
    .show();
};
