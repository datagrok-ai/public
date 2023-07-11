import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {IComputeDialogResult, IDescriptorTree} from '../types';
import '../../../css/hit-triage.css';

export async function chemDescriptorsDialog(onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  functions: {packageName: string, name: string}[], useDescriptors: boolean = true) {
  // tree groups
  const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
  const descriptorsGroup = ui.tree();
  const host = ui.div();
  host.style.maxHeight = '500px';
  host.style.maxWidth = '500px';
  descriptorsGroup.root.style.marginLeft = '16px';
  const descriptorItems: DG.TreeViewNode[] = [];

  function createTreeGroup(name: string, treeNode: DG.TreeViewGroup): DG.TreeViewGroup {
    const res = treeNode.group(name, null, false);
    res.enableCheckBox(false);
    return res;
  };

  // const descriptorsGroup = createTreeGroup('Descriptors', tree);
  const keys = Object.keys(descriptorsTree);
  for (const groupName of keys) {
    const group = createTreeGroup(groupName, descriptorsGroup);

    for (const descriptor of descriptorsTree[groupName].descriptors) {
      const item = group.item(descriptor.name, descriptor);
      descriptorItems.push(item);
      item.enableCheckBox(false);
    }
  };

  const descriptorsName = 'Descriptors';
  const funcNamesMap: {[key: string]: string} = {[descriptorsName]: descriptorsName};
  const tabControlArgs: {[key: string]: HTMLElement} =
  useDescriptors ? {[descriptorsName]: descriptorsGroup.root} : {};
  const funcInputsMap: {[funcName: string]: DG.FuncCall} = {};
  const calculatedFunctions: {[key: string]: boolean} = {[descriptorsName]: true};
  for (const func of functions) {
    const f = DG.Func.find({package: func.packageName, name: func.name})[0];
    const funcCall = f.prepare();
    const keyName = `${func.packageName}:${func.name}`;
    funcInputsMap[keyName] = funcCall;
    const editor = await funcCall.getEditor();
    editor.style.overflowY = 'scroll';
    tabControlArgs[f.friendlyName ?? f.name] = editor;
    funcNamesMap[f.friendlyName ?? f.name] = keyName;
    calculatedFunctions[keyName] = true;
    (editor.children[0] as HTMLElement).style.display = 'none'; // table input
    (editor.children[1] as HTMLElement).style.display = 'none'; // column input
  }
  const tc = ui.tabControl(tabControlArgs, true);
  host.appendChild(tc.root);
  // add checkboxes to each hader
  tc.panes.forEach((pane)=> {
    const functionCheck = ui.boolInput('', calculatedFunctions[funcNamesMap[pane.name]], () => {
      calculatedFunctions[funcNamesMap[pane.name]] = !!functionCheck.value;
    });
    functionCheck.setTooltip('Toggle calculation of this function');
    pane.header.appendChild(functionCheck.root);
    pane.header.classList.add('hit-triage-compute-dialog-pane-header');
  });
  function onOkProxy() {
    const res: IComputeDialogResult = {descriptors: [], externals: {}};
    res.descriptors =calculatedFunctions[descriptorsName] ?
      descriptorItems.filter((i) => i.checked).map((i) => i.value.name) : [];
    Object.entries(funcInputsMap).filter(([name, _]) => calculatedFunctions[name])
      .forEach(([name, value]) => {
        res.externals[name] = {};
        const inputs = value.inputs;
        Object.entries(inputs).slice(2).forEach(([inputName, inputVal]) => {
          res.externals[name][inputName] = inputVal;
        });
      });
    onOk(res);
  }

  ui.dialog('Compute')
    .add(host)
    .onOK(() => onOkProxy())
    .onCancel(onCancel)
    .show();
};
