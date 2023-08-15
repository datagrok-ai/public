import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';
import {IChemFunctionsDialogResult, IComputeDialogResult, IDescriptorTree, ITemplate} from '../types';
import '../../../css/hit-triage.css';

export async function chemFunctionsDialog(onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  template: ITemplate, dialog?: boolean,
): Promise<IChemFunctionsDialogResult> {
  const functions = template?.compute?.functions ?? [];
  const useDescriptors = !!template?.compute?.descriptors?.enabled;

  // tree groups
  const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
  const descriptorsGroup = ui.tree();
  descriptorsGroup.root.classList.add('hit-triage-compute-dialog-descriptors-group');
  const host = ui.div([], {classes: 'hit-triage-compute-dialog-host'});
  const descriptorItems: DG.TreeViewNode[] = [];

  function createTreeGroup(name: string, treeNode: DG.TreeViewGroup): DG.TreeViewGroup {
    const res = treeNode.group(name, null, false);
    res.enableCheckBox(false);
    return res;
  };

  // const descriptorsGroup = createTreeGroup('Descriptors', tree);
  const keys = Object.keys(descriptorsTree);
  const preselectedDescriptors: string[] = template?.compute?.descriptors?.args ?? [];
  for (const groupName of keys) {
    const group = createTreeGroup(groupName, descriptorsGroup);

    for (const descriptor of descriptorsTree[groupName].descriptors) {
      const item = group.item(descriptor.name, descriptor);
      descriptorItems.push(item);
      item.enableCheckBox(preselectedDescriptors.includes(descriptor.name));
    }
  };

  const descriptorsName = 'Descriptors';
  const funcNamesMap: {[key: string]: string} = {[descriptorsName]: descriptorsName};
  const tabControlArgs: {[key: string]: HTMLElement} =
  useDescriptors ? {[descriptorsName]: descriptorsGroup.root} : {};
  const funcInputsMap: {[funcName: string]: DG.FuncCall} = {};
  const calculatedFunctions: {[key: string]: boolean} = {[descriptorsName]: true};
  for (const func of functions) {
    const f = DG.Func.find({package: func.package, name: func.name})[0];
    const funcCall = f.prepare(func.args);
    const keyName = `${func.package}:${func.name}`;
    funcInputsMap[keyName] = funcCall;
    const editor = ui.div();
    const inputs = await funcCall.buildEditor(editor, {condensed: false});
    editor.classList.add('oy-scroll');
    editor.style.marginLeft = '15px';
    tabControlArgs[f.friendlyName ?? f.name] = editor;
    funcNamesMap[f.friendlyName ?? f.name] = keyName;
    calculatedFunctions[keyName] = true;
    (editor.children[0] as HTMLElement).style.display = 'none'; // table input
    (editor.children[1] as HTMLElement).style.display = 'none'; // column input
    inputs.forEach((input) => {
      if (input.property?.name && Object.keys(func.args).includes(input.property?.name))
        input.value = func.args[input.property.name];
      input.fireChanged();
    });
  }

  const tc = ui.tabControl(tabControlArgs, true);
  host.appendChild(tc.root);
  // add checkboxes to each hader
  tc.panes.forEach((pane)=> {
    const functionCheck = ui.boolInput('', calculatedFunctions[funcNamesMap[pane.name]], (v:boolean) => {
      calculatedFunctions[funcNamesMap[pane.name]] = !!functionCheck.value;
      if (!v){
        $(pane.content).find('input').attr('disabled', 'true');
      } else {
        $(pane.content).find('input').removeAttr('disabled');
      }
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

  if (dialog) {
    ui.dialog('Compute')
      .add(host)
      .onOK(() => onOkProxy())
      .onCancel(onCancel)
      .show();
  }

  return {
    root: host,
    okProxy: onOkProxy,
  };
};
