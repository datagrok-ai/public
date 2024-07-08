import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';
import {IChemFunctionsDialogResult, IComputeDialogResult, IDescriptorTree,
  HitTriageTemplate, HitTriageTemplateFunction, HitTriageTemplateScript} from '../types';
import '../../../css/hit-triage.css';
import {HTQueryPrefix, HTScriptPrefix} from '../consts';
import {HitAppBase} from '../hit-app-base';

export async function chemFunctionsDialog(app: HitAppBase<any>,
  onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  template: Omit<HitTriageTemplate, 'dataSourceType'>, dialog?: boolean,
): Promise<IChemFunctionsDialogResult> {
  // if compute is in dialog form, we need to show all the functions
  const computeFunctions = await app.computeFunctions;
  const functions = computeFunctions.functions
    .filter(({inputs: i}) => i.length > 2 && i[0].propertyType === 'dataframe' && i[1].propertyType === 'column')
    .map((f): HitTriageTemplateFunction => ({package: f.package.name, name: f.name, args:
      template?.compute?.functions?.find((tf) => tf.name === f.name && tf.package === f.package.name)?.args ?? {}}));

  const scripts: HitTriageTemplateScript[] = computeFunctions.scripts
    .filter(({inputs: i}) => i.length >= 2 && i[0].propertyType === 'dataframe' && i[1].propertyType === 'column')
    .map((s) => ({name: s.name, id: s.id, args: template?.compute?.scripts?.find((ts) => ts.id === s.id)?.args ?? {}}));

  const queries: HitTriageTemplateScript[] = computeFunctions.queries
    .filter((q) => q.inputs.length > 0 && q.inputs[0].propertyType === 'list')// for now only lists...
    .map((q) => {
      return {
        name: q.friendlyName ?? q.name,
        id: q.id,
        args: template?.compute?.queries?.find((ts) => ts.id === q.id)?.args ?? {},
      };
    })
    ;

  const useDescriptors = !!template?.compute?.descriptors?.enabled || dialog;

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
  // if compute is in dialog form, we need to show all the functions
  const tabControlArgs: {[key: string]: HTMLElement} =
  useDescriptors ? {[descriptorsName]: descriptorsGroup.root} : {};
  const funcInputsMap: {[funcName: string]: DG.FuncCall} = {};
  const calculatedFunctions: {[key: string]: boolean} = {[descriptorsName]: false};
  calculatedFunctions[descriptorsName] =
    !!template?.compute?.descriptors?.enabled && template?.compute?.descriptors?.args?.length > 0;
  // handling package functions
  for (const func of functions) {
    try {
      const f = DG.Func.find({package: func.package, name: func.name})[0];
      const funcCall = f.prepare(func.args);
      const keyName = `${func.package}:${func.name}`;
      funcInputsMap[keyName] = funcCall;
      const editor = ui.div();
      const inputs = await funcCall.buildEditor(editor, {condensed: false});
      editor.classList.add('oy-scroll');
      editor.style.marginLeft = '15px';
      editor.style.removeProperty('max-width');
      tabControlArgs[f.friendlyName ?? f.name] = editor;
      funcNamesMap[f.friendlyName ?? f.name] = keyName;
      calculatedFunctions[keyName] = template?.compute?.functions?.some(
        (f) => f.name === func.name && f.package === func.package,
      ) ?? false;
      (editor.children[0] as HTMLElement).style.display = 'none'; // table input
      (editor.children[1] as HTMLElement).style.display = 'none'; // column input
      inputs.forEach((input) => {
        if (input.property?.name && Object.keys(func.args).includes(input.property?.name))
          input.value = func.args[input.property.name];
        input.fireChanged();
      });
    } catch (e) {
      console.error(e);
      continue;
    }
  }

  // handling scripts
  for (const script of scripts) {
    try {
      const f = await grok.dapi.scripts.find(script.id);
      const funcCall = f.prepare(script.args);
      const keyName = `${HTScriptPrefix}:${script.name ?? ''}:${script.id}`;
      funcInputsMap[keyName] = funcCall;
      const editor = ui.div();
      const inputs = await funcCall.buildEditor(editor, {condensed: false});
      editor.classList.add('oy-scroll');
      editor.style.marginLeft = '15px';
      editor.style.removeProperty('max-width');
      tabControlArgs[f.friendlyName ?? f.name] = editor;
      funcNamesMap[f.friendlyName ?? f.name] = keyName;
      calculatedFunctions[keyName] = template?.compute?.scripts?.some(
        (s) => s.id === script.id,
      ) ?? false;
      (editor.children[0] as HTMLElement).style.display = 'none'; // table input
      (editor.children[1] as HTMLElement).style.display = 'none'; // column input
      inputs.forEach((input) => {
        if (input.property?.name && Object.keys(script.args).includes(input.property?.name))
          input.value = script.args[input.property.name];
        input.fireChanged();
      });
    } catch (e) {
      console.error(e);
      continue;
    }
  }

  // handling queries
  for (const query of queries) {
    try {
      const f = await grok.dapi.queries.find(query.id);
      const funcCall = f.prepare(query.args);
      const keyName = `${HTQueryPrefix}:${query.name ?? ''}:${query.id}`;
      funcInputsMap[keyName] = funcCall;
      const editor = ui.div();
      const inputs = await funcCall.buildEditor(editor, {condensed: false});
      editor.classList.add('oy-scroll');
      editor.style.marginLeft = '15px';
      editor.style.removeProperty('max-width');
      tabControlArgs[f.friendlyName ?? f.name] = editor;
      funcNamesMap[f.friendlyName ?? f.name] = keyName;
      calculatedFunctions[keyName] = template?.compute?.queries?.some(
        (s) => s.id === query.id,
      ) ?? false;
      (editor.children[0] as HTMLElement).style.display = 'none'; // list of molecules input
      inputs.forEach((input) => {
        if (input.property?.name && Object.keys(query.args).includes(input.property?.name))
          input.value = query.args[input.property.name];
        input.fireChanged();
      });
    } catch (e) {
      console.error(e);
      continue;
    }
  }

  const tc = ui.tabControl(tabControlArgs, true);
  tc.onTabChanged.subscribe(() => {try {tc.currentPane.content.style.removeProperty('max-width');} catch (e) {}});
  tc.header.style.overflow = 'scroll';
  tc.root.style.width = '100%';
  tc.root.style.minWidth = '350px';
  host.appendChild(tc.root);
  // add checkboxes to each hader
  tc.panes.forEach((pane)=> {
    const functionCheck = ui.input.bool('', {value: calculatedFunctions[funcNamesMap[pane.name]], onValueChanged: (input) => {
      calculatedFunctions[funcNamesMap[pane.name]] = !!functionCheck.value;
      if (!input.value)
        $(pane.content).find('input').attr('disabled', 'true');
      else
        $(pane.content).find('input').removeAttr('disabled');
    }});
    functionCheck.setTooltip('Toggle calculation of this function');
    pane.header.appendChild(functionCheck.root);
    pane.header.classList.add('hit-triage-compute-dialog-pane-header');
  });
  function onOkProxy() {
    const res: IComputeDialogResult = {descriptors: [], externals: {}, scripts: {}, queries: {}};
    res.descriptors = calculatedFunctions[descriptorsName] ?
      descriptorItems.filter((i) => i.checked).map((i) => i.value.name) : [];
    const filteredFuncs = Object.entries(funcInputsMap).filter(([name, _]) => calculatedFunctions[name]);
    // handling package functions
    filteredFuncs.filter(([name, _]) => !name.startsWith(HTScriptPrefix) && !name.startsWith(HTQueryPrefix))
      .forEach(([name, value]) => {
        res.externals[name] = {};
        const inputs = value.inputs;
        Object.entries(inputs).slice(2).forEach(([inputName, inputVal]) => {
          res.externals[name][inputName] = inputVal;
        });
      });
    // handling scripts
    filteredFuncs.filter(([name, _]) => name.startsWith(HTScriptPrefix))
      .forEach(([name, value]) => {
        res.scripts![name] = {};
        const inputs = value.inputs;
        Object.entries(inputs).slice(2).forEach(([inputName, inputVal]) => {
          res.scripts![name][inputName] = inputVal;
        });
      });
    // handling queries
    filteredFuncs.filter(([name, _]) => name.startsWith(HTQueryPrefix))
      .forEach(([name, value]) => {
        res.queries![name] = {};
        const inputs = value.inputs;
        Object.entries(inputs).slice(1).forEach(([inputName, inputVal]) => {
          res.queries![name][inputName] = inputVal;
        });
      });
    console.log(res);
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
