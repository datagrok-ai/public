import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

type IDescriptorTree = {
    [key: string]: {
      descriptors: Array<Descriptor>,
    } & Descriptor;
  }

  type Descriptor = {
    name: string,
    description: string,
  };

export async function chemDescriptorsDialog(onOk: (selected: string[]) => void, onCancel: () => void) {
  const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
  const tree = ui.tree();
  const groups: { [_: string]: DG.TreeViewGroup } = {};
  const items: DG.TreeViewNode[] = [];
  const selectedDescriptors: { [_: string]: string } = {};
  tree.root.style.maxHeight = '400px';

  const toggleAll = (val: boolean) => {
    for (const g of Object.values(groups))
      g.checked = val;
    for (const i of items)
      i.checked = val;
  };
  const selectAll = ui.label('All', {classes: 'd4-link-label', onClick: () => toggleAll(true)});
  selectAll.style.marginLeft = '6px';
  selectAll.style.marginRight = '12px';
  const selectNone = ui.label('None', {classes: 'd4-link-label', onClick: () => toggleAll(false)});
  const countLabel = ui.label('0 checked');
  countLabel.style.marginLeft = '24px';
  countLabel.style.display = 'inline-flex';

  const keys = Object.keys(descriptorsTree);
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

    for (const descriptor of descriptorsTree[groupName].descriptors) {
      const item = group.item(descriptor['name'], descriptor);
      items.push(item);
      item.enableCheckBox(false);
        item.checkBox!.onchange = (_e) => {
          countLabel.textContent = `${items.filter((i) => i.checked).length} checked`;
          if (item.checked)
            selectedDescriptors[item.text] = groupName;
        };
    }

    toggleAll(false);
  }

  ui.dialog('Select descriptors')
    .add(ui.divH([selectAll, selectNone, countLabel]))
    .add(tree.root)
    .onOK(() => onOk(items.filter((i) => i.checked).map((i: any) => i.value['name'])))
    .onCancel(onCancel)
    .show();
};
