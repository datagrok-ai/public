/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './../../css/biochem-functions-dialog-widget.css';

interface FunctionState {
  selected: boolean;
  funcCall: DG.FuncCall;
  editor: HTMLElement;
  navItem: HTMLElement;
  selectionCheckbox: DG.InputBase<boolean | null>;
  paramInputs: {[key: string]: DG.InputBase};
  tableInputName: string;
  moleculesInputName: string;
}

interface MethodInfo {
  name: string;
  package: string;
  authors: string;
  year: string;
  description: string;
  github?: string;
  citation?: string;
}


export async function biochemicalPropertiesDialog(): Promise<void> {
  const calculatorFuncs = await DG.Func.find({meta: {function_family: 'biochem-calculator'}});

  if (calculatorFuncs.length === 0) {
    grok.shell.warning('No biochemical calculators found.');
    return;
  }

  const dialog = ui.dialog({title: 'Biochemical Properties'});
  dialog.root.classList.add('biochem-calc-dialog');

  let table = grok.shell.t;
  let moleculesInput: DG.InputBase<DG.Column | null>;
  const moleculesInputDiv = ui.div();
  const onTableChanged = () => {
    moleculesInput = ui.input.column('Molecules', {table: table!, filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE});

    if (!moleculesInput.value && table) {
      const firstMoleculeCol = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
      if (firstMoleculeCol)
        moleculesInput.value = firstMoleculeCol;
    }

    ui.empty(moleculesInputDiv);
    moleculesInputDiv.appendChild(moleculesInput.root);
  };
  const tableInput = ui.input.table('Table', {value: table, items: grok.shell.tables, onValueChanged: (t) => {table = t!; onTableChanged();}});
  onTableChanged();
  const topInputContainer = ui.divH([tableInput.root, moleculesInputDiv], 'biochem-calc-top-inputs');
  const inputSection = ui.div([topInputContainer], 'biochem-calc-input-section');

  const functionState = new Map<string, FunctionState>();

  const searchInput = ui.input.string('', {placeholder: 'Search...'});
  const searchIcon = ui.iconFA('search', null, 'Search');
  searchIcon.classList.add('biochem-calc-search-icon');
  const navHeader = ui.div([searchIcon, searchInput.root], 'biochem-calc-nav-header');
  const navList = ui.div([], 'biochem-calc-nav-list');
  const navPanel = ui.div([navHeader, navList], 'biochem-calc-nav-panel');
  const editorPanel = ui.div([], 'biochem-calc-editor-panel');
  const methodFooter = ui.div([], 'biochem-calc-method-footer');
  const mainContent = ui.div([navPanel, editorPanel], 'biochem-calc-main-content');
  let firstFuncKey = '';

  const sanitizeValue = (value: any): any => {
    if (value instanceof DG.Column) return {_type: 'column', name: value.name};
    if (value instanceof DG.DataFrame) return {_type: 'dataframe', name: value.name};
    return value;
  };

  const desanitizeValue = (value: any): any => {
    if (value && value._type) {
      if (value._type === 'column' && table) return table.col(value.name);
      if (value._type === 'dataframe') return grok.shell.tableByName(value.name);
    }
    return value;
  };


  const getStringInput = (): string => {
    const stateToSave: {[key: string]: any} = {};
    stateToSave.globalInputs = {
      table: table?.name || null,
      moleculesColumn: moleculesInput?.value?.name || null,
    };
    stateToSave.functions = {};
    for (const [funcName, state] of functionState.entries()) {
      const currentParamValues: {[key: string]: any} = {};
      for (const [paramName, input] of Object.entries(state.paramInputs))
        currentParamValues[paramName] = sanitizeValue(input.value);
      stateToSave.functions[funcName] = {
        selected: state.selectionCheckbox.value,
        params: currentParamValues,
      };
    }
    return JSON.stringify(stateToSave);
  };

  const applyStringInput = (input: string): void => {
    try {
      const stateToRestore = JSON.parse(input);
      if (stateToRestore.globalInputs) {
        const globalInputs = stateToRestore.globalInputs;
        if (globalInputs.table) {
          const targetTable = grok.shell.tables.find((t) => t.name === globalInputs.table);
          if (targetTable) {
            tableInput.value = targetTable;
            table = targetTable;
            onTableChanged();
          }
        }
        if (globalInputs.moleculesColumn && table) {
          const targetColumn = table.col(globalInputs.moleculesColumn);
          if (targetColumn)
            moleculesInput.value = targetColumn;
        }
      }
      if (stateToRestore.functions) {
        for (const [funcName, savedState] of Object.entries(stateToRestore.functions)) {
          const state = functionState.get(funcName);
          if (state && savedState) {
            const typedSavedState = savedState as any;
            state.selectionCheckbox.value = typedSavedState.selected;
            state.selected = typedSavedState.selected;
            const savedParams = typedSavedState.params;
            for (const [paramName, paramValue] of Object.entries(savedParams)) {
              if (state.paramInputs[paramName]) {
                const restoredValue = desanitizeValue(paramValue);
                state.paramInputs[paramName].value = restoredValue;
              }
            }
          }
        }
      }
    } catch (e) {
      console.error('Failed to restore state from history:', e);
      grok.shell.error('Failed to restore state from history');
    }
  };


  const getMethodInfoFromMeta = (func: DG.Func): MethodInfo => {
    const opts = func.options || {};
    return {
      name: func.friendlyName || func.name, package: opts['method_info.package'] || 'N/A',
      authors: opts['method_info.author'] || 'N/A', year: opts['method_info.year'] || '',
      description: func.description || '', github: opts['method_info.github'], citation: opts['method_info.citation'],
    };
  };

  const buildFunctionEditors = async () => {
    for (const func of calculatorFuncs) {
      if (!func.inputs || func.inputs.length < 2) {
        console.warn(`Biochem-Widget: Skipping function "${func.name}" because it has fewer than two inputs.`);
        continue;
      }
      const tableInputDef = func.inputs[0];
      const moleculesInputDef = func.inputs[1];

      if (!tableInputDef || !moleculesInputDef) {
        console.warn(`Biochem-Widget: Skipping function "${func.name}" because its first or second input is malformed (null).`);
        continue;
      }
      if (tableInputDef.propertyType !== 'dataframe' || moleculesInputDef.propertyType !== 'column') {
        console.warn(`Biochem-Widget: Skipping function "${func.name}" because its first two inputs are not a dataframe and a column.`);
        continue;
      }

      const tableInputName = tableInputDef.name;
      const moleculesInputName = moleculesInputDef.name;

      const funcName = (func.friendlyName ?? func.name).replace(/^Calculate\s+/i, '');
      if (!firstFuncKey) firstFuncKey = funcName;
      const funcCall = func.prepare({});
      const paramInputs: {[key: string]: DG.InputBase} = {};
      const editorContainer = ui.divV([], 'ui-form');
      editorContainer.style.padding = '15px';

      const inputs = await funcCall.buildEditor(editorContainer);
      inputs.forEach((input) => {
        if (input.property)
          paramInputs[input.property.name] = input;
      });

      if (inputs[0]) inputs[0].root.style.display = 'none';
      if (inputs[1]) inputs[1].root.style.display = 'none';

      const checkbox = ui.input.bool('', {value: false});

      inputs.forEach((input) => {
        input.onChanged.subscribe(() => {
          if (checkbox.value === false)
            checkbox.value = true;
        });
      });

      const navItem = ui.div([checkbox.root, ui.span([funcName])], 'biochem-calc-nav-item');

      const state: FunctionState = {
        selected: false, funcCall, editor: editorContainer, navItem,
        selectionCheckbox: checkbox, paramInputs,
        tableInputName: tableInputName,
        moleculesInputName: moleculesInputName,
      };

      checkbox.onChanged.subscribe((v) => {state.selected = v;});
      functionState.set(funcName, state);
      navList.appendChild(navItem);
      navItem.addEventListener('click', (e) => {if (e.target !== checkbox.input) setActiveFunction(funcName);});
    }
  }; ;

  await buildFunctionEditors();

  const setActiveFunction = (funcName: string) => {
    functionState.forEach((state, name) => state.navItem.classList.toggle('active', name === funcName));
    ui.empty(editorPanel);
    const activeState = functionState.get(funcName);
    if (activeState) {
      editorPanel.appendChild(activeState.editor);
      updateMethodInfoFooter(funcName, methodFooter, calculatorFuncs);
    }
  };

  if (firstFuncKey) setActiveFunction(firstFuncKey);

  const searchInFunctionMetadata = (func: DG.Func, searchTerm: string): boolean => {
    searchTerm = searchTerm.toLowerCase();
    if (func.name.toLowerCase().includes(searchTerm) || func.friendlyName?.toLowerCase().includes(searchTerm) || func.description?.toLowerCase().includes(searchTerm))
      return true;
    for (const [key, value] of Object.entries(func.options || {})) {
      if (typeof value === 'string' && value.toLowerCase().includes(searchTerm)) return true;
      if (key.toLowerCase().includes(searchTerm)) return true;
    }
    return false;
  };

  searchInput.onChanged.subscribe(() => {
    const searchTerm = searchInput.value;
    for (const [funcName, state] of functionState.entries()) {
      const func = calculatorFuncs.find((f) => (f.friendlyName ?? f.name).replace(/^Calculate\s+/i, '') === funcName)!;
      state.navItem.style.display = searchInFunctionMetadata(func, searchTerm) ? 'flex' : 'none';
    }
  });

  dialog.add(ui.divV([inputSection, mainContent, methodFooter]));

  dialog.onOK(async () => {
    const selectedFunctions = Array.from(functionState.entries()).filter(([, state]) => state.selected);
    if (selectedFunctions.length === 0) return grok.shell.warning('Please select at least one calculation.');
    if (!moleculesInput.value) return grok.shell.warning('Please select a molecule column.');
    await executeSelectedFunctions(selectedFunctions, table!, moleculesInput.value);
  });

  const getHistoryInput = () => {
    const selectedFunctions = Array.from(functionState.entries())
      .filter(([, state]) => state.selected)
      .map(([funcName]) => funcName);

    const infoString = selectedFunctions.length > 0 ? selectedFunctions.join(', ') : 'No calculations';
    return {
      info: infoString,
      editorSettings: getStringInput(),
    };
  };

  dialog.history(getHistoryInput, (x: any) => applyStringInput(x['editorSettings']));

  dialog.show({width: 900, height: 680, resizable: true});


  async function executeSelectedFunctions(
    selectedFunctions: [string, FunctionState][],
    table: DG.DataFrame,
    molecules: DG.Column,
  ) {
    const progress = DG.TaskBarProgressIndicator.create('Calculating properties...');
    let completedCount = 0;
    try {
      for (const [funcName, state] of selectedFunctions) {
        progress.update((completedCount / selectedFunctions.length) * 100, `Calculating ${funcName}...`);
        state.funcCall.inputs[state.tableInputName] = table;
        state.funcCall.inputs[state.moleculesInputName] = molecules;

        await state.funcCall.call(undefined, undefined, {processed: false});
        completedCount++;
      }
    } catch (error: unknown) {
      grok.shell.error(`Calculation failed: ${error instanceof Error ? error.message : String(error)}`);
    } finally {
      progress.close();
      if (completedCount > 0) grok.shell.info(`Successfully calculated ${completedCount} properties.`);
    }
  }

  function updateMethodInfoFooter(paneName: string, container: HTMLElement, funcs: DG.Func[]) {
    ui.empty(container);
    const func = funcs.find((f) => (f.friendlyName ?? f.name).replace(/^Calculate\s+/i, '') === paneName);
    if (!func) return;

    const methodInfo = getMethodInfoFromMeta(func);
    const gridItems = [
      ui.div([ui.span(['Package: '], 'font-weight-bold'), ui.span([methodInfo.package])]),
      ui.div([
        ui.span(['Authors: '], 'font-weight-bold'),
        ui.span([methodInfo.year ? `${methodInfo.authors} (${methodInfo.year})` : methodInfo.authors]),
      ]),
    ];

    if (methodInfo.github) {
      gridItems.push(ui.div([
        ui.span(['Source: '], 'font-weight-bold'),
        ui.link('GitHub', () => window.open(methodInfo.github!, '_blank'), 'Link to source code repository'),
      ]));
    }
    const infoGrid = ui.div(gridItems, 'biochem-calc-method-footer-grid');
    const desc = ui.divText(func.description || methodInfo.description, {style: {marginTop: '5px'}});
    container.appendChild(infoGrid);
    container.appendChild(desc);
  }
}
