import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {IComputeDialogResult, IDescriptorTree, IFunctionInput, IFunctionProps} from './types';

function getArgsFromFuncName(packageName: string, name: string) {
  const func = DG.Func.find({package: packageName, name: name})[0];
  const funcObj: IFunctionProps = {
    name: func.name,
    package: packageName,
    description: func.description,
    calculate: false,
    displayName: func.friendlyName ?? func.name,
    inputs: func.inputs.slice(2)// first two args are table and column
      .map(({caption, choices, defaultValue, description, name, propertyType}): IFunctionInput => ({
        name,
        displayName: caption ?? name,
        description,
        type: propertyType,
        defaultValue,
        choices,
        value: defaultValue,
      })),
    onlyCheckboxes: func.inputs.slice(2).every((input) => input.propertyType === DG.TYPE.BOOL),
  };
  return funcObj;
}

function getInput(inputProps: IFunctionInput, dataframe: DG.DataFrame, onChange: (val: any) => void) {
  if (inputProps.choices)
    return ui.choiceInput(inputProps.displayName, inputProps.defaultValue, inputProps.choices, onChange);

  switch (inputProps.type) {
  case DG.TYPE.BOOL:
    return ui.boolInput(inputProps.displayName, inputProps.defaultValue ?? false, onChange);
  case DG.TYPE.STRING:
    return ui.stringInput(inputProps.displayName, inputProps.defaultValue ?? null, onChange);
  case DG.TYPE.FLOAT:
    return ui.floatInput(inputProps.displayName, inputProps.defaultValue ?? null, onChange);
  case DG.TYPE.INT:
  case DG.TYPE.BIG_INT:
    return ui.intInput(inputProps.displayName, inputProps.defaultValue ?? null, onChange);
  case DG.TYPE.COLUMN:
    return ui.columnInput(inputProps.displayName, dataframe, inputProps.defaultValue ?? null, onChange);
  case DG.TYPE.COLUMN_LIST:
    return ui.columnsInput(inputProps.displayName, dataframe, onChange);
  default:
    throw new Error(`Unsupported input type: ${inputProps.type}`);
  }
}

function assembleCheckboxTree(treeNode: DG.TreeViewGroup, functionView: IFunctionProps) {
  function onSelectionChange() {
    let isAnyChecked = false;
    treeNode.items.forEach((child, i) => {
      functionView.inputs[i].value = child.checked;
      isAnyChecked ||= child.checked;
    });
    functionView.calculate = isAnyChecked || treeNode.checked;
  }

  treeNode.enableCheckBox(false);
  treeNode.checkBox!.onchange = () => {
    treeNode.items.forEach((child) => child.checked = treeNode.checked); onSelectionChange();
  };

  functionView.inputs.forEach((prop) => {
    const item = treeNode.item(prop.displayName, prop);
    item.enableCheckBox(false);
    item.checkBox!.onchange = () => onSelectionChange();
    prop.description && ui.tooltip.bind(item.root, prop.description);
  });
}

function assembleFunctionTree(treeNode: DG.TreeViewGroup, functionView: IFunctionProps, dataframe: DG.DataFrame): void {
  const group = treeNode.group(functionView.displayName, null, false);

  if (functionView.onlyCheckboxes) {
    assembleCheckboxTree(group, functionView);
    return;
  }
  group.enableCheckBox(false);
  group.checkBox!.onchange = (_e) => {functionView.calculate = group.checked;};
  const groupLabel = group.root.getElementsByClassName('d4-tree-view-node')[0] as HTMLElement;
  functionView.description && groupLabel && ui.tooltip.bind(groupLabel, functionView.description);
  functionView.inputs.forEach((prop, i) => {
    const inputView = getInput(prop, dataframe,
      (val: any) => {
        functionView.inputs[i].value = val;
        functionView.calculate = true;
        group.checked = true;
      });
    group.root.getElementsByClassName('d4-tree-view-group-host')[0].appendChild(inputView.root);
  });
}


export async function chemDescriptorsDialog(onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  dataframe: DG.DataFrame, functions: {packageName: string, name: string}[]) {
  // tree groups
  const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
  const tree = ui.tree();
  tree.root.style.maxHeight = '500px';
  tree.root.style.width = '400px';

  const descriptorItems: DG.TreeViewNode[] = [];

  function createTreeGroup(name: string, treeNode: DG.TreeViewGroup): DG.TreeViewGroup {
    const res = treeNode.group(name, null, false);
    res.enableCheckBox(false);
    return res;
  };

  const descriptorsGroup = createTreeGroup('Descriptors', tree);
  const keys = Object.keys(descriptorsTree);
  for (const groupName of keys) {
    const group = createTreeGroup(groupName, descriptorsGroup);

    for (const descriptor of descriptorsTree[groupName].descriptors) {
      const item = group.item(descriptor.name, descriptor);
      descriptorItems.push(item);
      item.enableCheckBox(false);
    }
  };
  const externalFunctions = functions.map((func) => {
    const args = getArgsFromFuncName(func.packageName, func.name);
    assembleFunctionTree(tree, args, dataframe);
    return args;
  });

  function onOkProxy() {
    const res: IComputeDialogResult = {descriptors: [], externals: {}};
    res.descriptors = descriptorItems.filter((i) => i.checked).map((i) => i.value.name);
    externalFunctions.filter((func) => func.calculate).forEach((func) => {
      const funcName = func.name;
      const packageName = func.package;
      const funcArgs = func.inputs.map((i) => ({name: i.name, value: i.value}));
      res.externals[`${packageName}:${funcName}`] = {};
      funcArgs.forEach((arg) => res.externals[`${packageName}:${funcName}`][arg.name] = arg.value);
    });
    onOk(res);
  }

  ui.dialog('Calculate')
    .add(tree.root)
    .onOK(() => onOkProxy())
    .onCancel(onCancel)
    .show();
};
