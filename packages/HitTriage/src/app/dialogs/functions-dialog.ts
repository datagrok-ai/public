import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';
import {IChemFunctionsDialogResult, IComputeDialogResult, IDescriptorTree,
  HitTriageTemplate,
  IFunctionArgs} from '../types';
import '../../../css/hit-triage.css';
import {funcTypeNames, HTQueryPrefix, HTScriptPrefix} from '../consts';
import {HitAppBase} from '../hit-app-base';
import { FunctionOrdering, getReorderingInput, getSavedFunctionOrdering } from '../utils';

export async function chemFunctionsDialog(app: HitAppBase<any>,
  onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  template: Omit<HitTriageTemplate, 'dataSourceType'>, dialog?: boolean,
): Promise<IChemFunctionsDialogResult> {
  // if compute is in dialog form, we need to show all the functions

  const computeFunctions = app.computeFunctions;

  const allFunctions: {func: DG.Func | DG.Script | DG.DataQuery, args: IFunctionArgs}[] = [];
  computeFunctions.functions
    .filter(({inputs: i}) => i.length >= 2 && i[0].propertyType === 'dataframe' && i[1].propertyType === 'column')
    .forEach((f) => {
      allFunctions.push({
        func: f,
        args: template?.compute?.functions?.find((tf) =>
          tf.name === f.name && tf.package === f.package.name)?.args ?? {}});
    });
  computeFunctions.scripts
    .filter(({inputs: i}) => i.length >= 2 && i[0].propertyType === 'dataframe' && i[1].propertyType === 'column')
    .forEach((s) => {
      allFunctions.push({
        func: s,
        args: template?.compute?.scripts?.find((ts) => ts.id === s.id)?.args ?? {}});
    });
  computeFunctions.queries
    .filter((q) => q.inputs.length > 0 && q.inputs[0].propertyType === 'list')// for now only lists...
    .forEach((q) => {
      allFunctions.push({
        func: q,
        args: template?.compute?.queries?.find((ts) => ts.id === q.id)?.args ?? {}});
    });

  const host = ui.div([], {classes: 'hit-triage-compute-dialog-host'});
  const descriptorItems: DG.TreeViewNode[] = [];
  const descriptorsGroup = ui.tree();
  const descriptorsGroupHost = ui.wait(async () => {
    try {
      // tree groups
        const descriptorsTree = (await grok.chem.descriptorsTree()) as IDescriptorTree;
        descriptorsGroup.root.classList.add('hit-triage-compute-dialog-descriptors-group');
    
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
        return descriptorsGroup.root;
      } catch (e) {
        console.error(e);
      }
      return ui.divText('Failed to load descriptors');
  })
  descriptorsGroupHost.classList.add('oy-scroll');
  
  const descriptorsName = 'Descriptors';
  const funcNamesMap: {[key: string]: string} = {[descriptorsName]: descriptorsName};
  // if compute is in dialog form, we need to show all the functions
  const tabControlArgs: {[key: string]: HTMLElement} = {[descriptorsName]: descriptorsGroupHost};
  const funcInputsMap: {[funcName: string]: DG.FuncCall} = {};
  const calculatedFunctions: {[key: string]: boolean} = {[descriptorsName]: false};
  calculatedFunctions[descriptorsName] =
    !!template?.compute?.descriptors?.enabled && template?.compute?.descriptors?.args?.length > 0;
  // handling package functions
  for (const func of allFunctions) {
    try {
      const f = func.func;
      const funcCall = f.prepare(func.args);
      const keyName = f.type === funcTypeNames.function ?
        `${f.package.name}:${f.name}`: f.type === funcTypeNames.script ?
          `${HTScriptPrefix}:${f.name ?? ''}:${f.id}` : `${HTQueryPrefix}:${f.name ?? ''}:${f.id}`;
      funcInputsMap[keyName] = funcCall;
      const editor = ui.div();
      ui.setUpdateIndicator(editor, true);
      // deffered action to not wait for the editor to be built
      funcCall.buildEditor(editor, {condensed: false}).then((inputs) => {
        inputs?.[0]?.root && (inputs[0].root.style.display = 'none'); // table input
        // for queries we don't have column input, we only have one input for the list of molecules
        inputs?.[1]?.root && f.type !== funcTypeNames.query && (inputs[1].root.style.display = 'none'); // column input

        inputs.forEach((input) => {
          if (input.property?.name && Object.keys(func.args).includes(input.property?.name))
            input.value = func.args[input.property.name];
          input.fireChanged();
        });
      }).catch((e) => {
        console.error(e);
        ui.setUpdateIndicator(editor, false);
        editor.appendChild(ui.divText('Failed to build function editor'));
      }).finally(() => ui.setUpdateIndicator(editor, false));
      editor.classList.add('oy-scroll');
      editor.style.marginLeft = '15px';
      editor.style.removeProperty('max-width');
      tabControlArgs[f.friendlyName ?? f.name] = editor;
      funcNamesMap[f.friendlyName ?? f.name] = keyName;
      calculatedFunctions[keyName] = (template?.compute?.functions?.some(
        (f) => func.func.type === funcTypeNames.function &&
        f.name === func.func.name && f.package === func.func.package?.name) ||
          template?.compute?.scripts?.some((s) => s.id === f.id) ||
          template?.compute?.queries?.some((s) => s.id === f.id)) ?? false;
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
  function addCheckboxesToPaneHeaders () {
    tc.panes.forEach((pane)=> {
      const functionCheck =
        ui.input.bool('', {value: calculatedFunctions[funcNamesMap[pane.name]], onValueChanged: (value) => {
          calculatedFunctions[funcNamesMap[pane.name]] = !!value;
          if (!value)
            $(pane.content).find('input')?.attr('disabled', 'true');
          else
            $(pane.content).find('input')?.removeAttr('disabled');
        }});
      functionCheck.setTooltip('Toggle calculation of this function');
      pane.header.appendChild(functionCheck.root);
      pane.header.classList.add('hit-triage-compute-dialog-pane-header');
    });
  }

  function initTabs(order: FunctionOrdering) {
    tc.clear();
    let actOrder = order.order ?? Object.keys(tabControlArgs);
    if (actOrder.length === 0)
      actOrder = Object.keys(tabControlArgs);
    
    actOrder.forEach((n) => {
      if (tabControlArgs[n])
        tc.addPane(n, () => tabControlArgs[n]);
    });
    // after adding ordered panes, we also need to add new functions, not included in the order
    Object.keys(tabControlArgs).forEach((n) => {
      if (!actOrder.includes(n) && tabControlArgs[n] && !(order.hidden ?? []).includes(n)) {
        tc.addPane(n, () => tabControlArgs[n]);
      }
    });
    // make sure that hidden functions are not checked;
    (order.hidden ?? []).forEach((n) => {
      calculatedFunctions[funcNamesMap[n]] = false;
    });
    addCheckboxesToPaneHeaders();
    if (tc.panes.length > 0)
      tc.currentPane = tc.panes[0];
  }

  initTabs(getSavedFunctionOrdering());
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

  tc.root.appendChild(getReorderingInput(Object.keys(tabControlArgs), (newOrder) => {
    initTabs(newOrder);
  }));
  if (dialog) {
    ui.dialog('Compute')
      .add(host)
      .onOK(() => onOkProxy())
      .onCancel(onCancel)
      .show({resizable: true});
  }

  return {
    root: host,
    okProxy: onOkProxy,
  };
};
