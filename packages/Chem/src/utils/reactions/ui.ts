/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {renderReactionToCanvas} from '../../rendering/rdkit-reaction-renderer';
import {getRdKitModule} from '../chem-common-rdkit';
import {resolveReactionVariables, runTransformationReaction, runTwoComponentReaction} from './reactions';
import {ReactionBrowser} from './reaction-browser';
import {openReactionEditor} from './reaction-editor';
import {NamedReaction} from './types';

// ==========================================
//  Shared component types
// ==========================================

type TransformationComponents = {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  columnContainer: HTMLElement;
  saltsInput: DG.InputBase<boolean | null>;
  outputNameInput: DG.InputBase<string>;
  browser: ReactionBrowser;
  newReactionButton: HTMLElement;
  previewCanvas: HTMLCanvasElement;
  variablesDiv: HTMLElement;
  exampleGrid: ExampleGridInfo;
  getColumnInput: () => DG.InputBase<DG.Column | null>;
  getResolvedSmarts: () => string | null;
  runOnExistingTable: () => Promise<void>;
  runToNewDataFrame: () => Promise<DG.DataFrame | null>;
};

type TwoComponentComponents = {
  table1Input: DG.InputBase<DG.DataFrame | null>;
  col1Container: HTMLElement;
  table2Input: DG.InputBase<DG.DataFrame | null>;
  col2Container: HTMLElement;
  modeInput: DG.InputBase<string | null>;
  saltsInput: DG.InputBase<boolean | null>;
  browser: ReactionBrowser;
  newReactionButton: HTMLElement;
  previewCanvas: HTMLCanvasElement;
  variablesDiv: HTMLElement;
  exampleGrid: ExampleGridInfo;
  getColumn1Input: () => DG.InputBase<DG.Column | null>;
  getColumn2Input: () => DG.InputBase<DG.Column | null>;
  getResolvedSmarts: () => string | null;
  runOnExistingTable: () => Promise<void>;
  runToNewDataFrame: () => Promise<DG.DataFrame | null>;
};

type ExampleGridInfo = {root: HTMLElement; grid: DG.Grid | null; df: DG.DataFrame};

/** Create a column input that gracefully handles null table (no table loaded yet). */
function createColumnInput(
  label: string, table: DG.DataFrame | null | undefined,
  column: DG.Column | null | undefined, tooltip: string,
): DG.InputBase<DG.Column | null> {
  const opts: any = {
    value: column ?? undefined,
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
  };
  if (table) opts.table = table;
  const inp = ui.input.column(label, opts);
  ui.tooltip.bind(inp.input, tooltip);
  return inp;
}

// ==========================================
//  Builder: transformation (single column)
// ==========================================

function buildTransformationComponents(
  rdkit: RDModule, table?: DG.DataFrame | null, column?: DG.Column | null,
): TransformationComponents {
  table = table ?? grok.shell.t ?? null;
  column = column ?? table?.columns?.bySemType(DG.SEMTYPE.MOLECULE) ?? null;

  // ---- Inputs ----
  const tableInput = ui.input.table('Table', {value: table ?? undefined, nullable: true});
  ui.tooltip.bind(tableInput.input, 'Table containing the molecule column to transform');

  const columnContainer = ui.div();
  let columnInput = createColumnInput('Molecules', table, column,
    'Column with molecules to run the reaction on');
  columnContainer.append(columnInput.root);

  const saltsInput = ui.input.bool('Remove salts and water', {value: true});
  ui.tooltip.bind(saltsInput.input, 'Remove common salts and water from molecules before running the reaction');

  const outputNameInput = ui.input.string('Output column', {
    value: table?.columns?.getUnusedName(`Reacted(${column?.name ?? 'molecules'})`) ?? 'Products',
    nullable: false,
  });
  ui.tooltip.bind(outputNameInput.input, 'Name of the new column that will contain the reaction products');
  outputNameInput.root.style.flexGrow = '1';
  outputNameInput.root.style.minWidth = '0';
  outputNameInput.root.style.maxWidth = '399px';

  tableInput.onChanged.subscribe(async () => {
    const t = tableInput.value;
    await t?.meta?.detectSemanticTypes();
    await DG.delay(100); // wait for semantic type detection to finish and UI to update
    const molCol = t?.columns?.bySemType(DG.SEMTYPE.MOLECULE) ?? null;
    columnInput = createColumnInput('Molecules', t, molCol,
      'Column with molecules to run the reaction on');
    columnContainer.replaceChildren(columnInput.root);
  });

  // ---- Browser ----
  const browser = new ReactionBrowser(rdkit, {
    modeFilter: 'transformation',
    showDetails: true,
    allowAdd: true,
    allowEdit: true,
  });

  // ---- Variables ----
  const variablesDiv = ui.divV([], {style: {paddingTop: '4px'}});
  let currentVariableInputs: DG.InputBase[] = [];

  // ---- Preview ----
  const r = window.devicePixelRatio;
  const previewCanvas = ui.canvas(400 * r, 100 * r);
  previewCanvas.style.width = '400px';
  previewCanvas.style.maxWidth = '100%';
  previewCanvas.style.height = '100px';
  previewCanvas.style.border = '1px solid var(--grey-2)';
  previewCanvas.style.borderRadius = '4px';
  previewCanvas.style.cursor = 'pointer';
  previewCanvas.addEventListener('click', () => openFullscreenPreview(rdkit, browser.selected));

  const exampleGrid = createExampleGrid();

  function getResolvedVariables(): Record<string, any> {
    const reaction = browser.selected;
    const vars: Record<string, any> = {};
    if (reaction?.variables) {
      Object.keys(reaction.variables).forEach((key, i) => {
        if (currentVariableInputs[i])
          vars[key] = currentVariableInputs[i].value;
      });
    }
    return vars;
  }

  function getResolvedSmarts(): string | null {
    const reaction = browser.selected;
    if (!reaction) return null;
    const vars = getResolvedVariables();
    return reaction.variables ?
      resolveReactionVariables(reaction.reactionSmarts, vars) :
      reaction.reactionSmarts;
  }

  function updatePreview() {
    const reaction = browser.selected;
    if (!reaction) return;
    const smarts = getResolvedSmarts();
    if (!smarts) return;

    const ctx = previewCanvas.getContext('2d')!;
    ctx.clearRect(0, 0, previewCanvas.width, previewCanvas.height);
    try {
      renderReactionToCanvas(rdkit, previewCanvas, smarts, previewCanvas.width, previewCanvas.height);
    } catch {/* preview error */}

    const col = columnInput.value;
    if (col) {
      const examples = col.categories.slice(0, 10).filter((s: string) => !!s);
      if (examples.length > 0) {
        runTransformationReaction(rdkit, smarts, examples, {removeSaltsAndWater: saltsInput.value}).then((result) => {
          updateExampleGrid(exampleGrid, examples, result.products);
        }).catch(() => {/* silently fail preview */});
      }
    }
  }

  // Wire browser selection → variables + preview
  browser.onSelectionChanged.subscribe((reaction) => {
    variablesDiv.innerHTML = '';
    currentVariableInputs = [];
    if (!reaction) return;

    if (reaction.variables) {
      for (const [_key, varDef] of Object.entries(reaction.variables)) {
        const inputCtor = varDef.type === 'string' ? ui.input.string :
          varDef.type === 'int' ? ui.input.int : ui.input.float;
        const inp = inputCtor(varDef.name, {
          value: varDef.defaultValue,
          nullable: false,
          onValueChanged: () => updatePreview(),
        });
        if (varDef.description)
          ui.tooltip.bind(inp.input, varDef.description);
        currentVariableInputs.push(inp);
        variablesDiv.append(inp.root);
      }
    }
    updatePreview();
  });

  const newReactionButton = ui.button('+ New Reaction', () => {
    openReactionEditor(rdkit, {
      onSave: async () => {await browser.reload();},
    });
  }, 'Create a new custom reaction');

  // ---- Run: add column to existing table ----
  async function runOnExistingTable(): Promise<void> {
    const reaction = browser.selected;
    if (!reaction) {grok.shell.warning('No reaction selected.'); return;}
    const col = columnInput.value;
    const t = tableInput.value;
    if (!col || !t) {grok.shell.warning('No table or molecule column selected.'); return;}

    const smarts = getResolvedSmarts()!;
    const pg = DG.TaskBarProgressIndicator.create('Running reaction...');
    try {
      const result = await runTransformationReaction(rdkit, smarts, col.toList(), {
        removeSaltsAndWater: saltsInput.value,
        showInfo: true,
      });
      const newCol = t.columns.addNewString(t.columns.getUnusedName(outputNameInput.value || `Reacted(${col.name})`));
      newCol.semType = DG.SEMTYPE.MOLECULE;
      newCol.setTag('cell.renderer', 'Molecule');
      newCol.init((i) => result.products[i]);
    } catch (e: any) {
      grok.shell.error(`Reaction failed: ${e?.message ?? e}`);
    } finally {
      pg.close();
    }
  }

  // ---- Run: create new DataFrame ----
  async function runToNewDataFrame(): Promise<DG.DataFrame | null> {
    const reaction = browser.selected;
    if (!reaction) {grok.shell.warning('No reaction selected.'); return null;}
    const col = columnInput.value;
    if (!col) {grok.shell.warning('No molecule column selected.'); return null;}

    const smarts = getResolvedSmarts()!;
    const pg = DG.TaskBarProgressIndicator.create('Running reaction...');
    try {
      const inputs = col.toList();
      const result = await runTransformationReaction(rdkit, smarts, inputs, {
        removeSaltsAndWater: saltsInput.value,
        showInfo: true,
      });
      const df = DG.DataFrame.create(inputs.length);
      df.name = `${reaction.name} Products`;
      const inputCol = df.columns.addNewString('Input');
      inputCol.semType = DG.SEMTYPE.MOLECULE;
      inputCol.setTag('cell.renderer', 'Molecule');
      inputCol.init((i) => inputs[i] ?? '');
      const prodCol = df.columns.addNewString('Product');
      prodCol.semType = DG.SEMTYPE.MOLECULE;
      prodCol.setTag('cell.renderer', 'Molecule');
      prodCol.init((i) => result.products[i] ?? '');
      return df;
    } catch (e: any) {
      grok.shell.error(`Reaction failed: ${e?.message ?? e}`);
      return null;
    } finally {
      pg.close();
    }
  }

  return {
    tableInput, columnContainer, saltsInput, outputNameInput,
    browser, newReactionButton,
    previewCanvas, variablesDiv, exampleGrid,
    getColumnInput: () => columnInput,
    getResolvedSmarts,
    runOnExistingTable, runToNewDataFrame,
  };
}

// ==========================================
//  Builder: two-component (two columns)
// ==========================================

function buildTwoComponentComponents(
  rdkit: RDModule, table?: DG.DataFrame | null,
): TwoComponentComponents {
  table = table ?? grok.shell.t ?? null;
  const molColumns = table?.columns?.toList().filter((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE) ?? [];

  // ---- Inputs ----
  const table1Input = ui.input.table('Table 1', {value: table ?? undefined, nullable: true});
  ui.tooltip.bind(table1Input.input, 'Table containing the first reactant column');
  const col1Container = ui.div();
  let column1Input = createColumnInput('Reactant 1', table, molColumns[0],
    'First reactant molecule column');
  col1Container.append(column1Input.root);

  const table2Input = ui.input.table('Table 2', {value: table ?? undefined, nullable: true});
  ui.tooltip.bind(table2Input.input, 'Table containing the second reactant column');
  const col2Container = ui.div();
  let column2Input = createColumnInput('Reactant 2', table,
    molColumns.length > 1 ? molColumns[1] : molColumns[0],
    'Second reactant molecule column');
  col2Container.append(column2Input.root);

  table1Input.onChanged.subscribe(async () => {
    const t = table1Input.value;
    await t?.meta?.detectSemanticTypes();
    await DG.delay(100); // wait for semantic type detection to finish and UI to update
    const molCol = t?.columns?.bySemType(DG.SEMTYPE.MOLECULE) ?? null;
    column1Input = createColumnInput('Reactant 1', t, molCol,
      'First reactant molecule column');
    col1Container.replaceChildren(column1Input.root);
  });
  table2Input.onChanged.subscribe(async () => {
    const t = table2Input.value;
    await t?.meta?.detectSemanticTypes();
    await DG.delay(100); // wait for semantic type detection to finish and UI to update
    const molCol = t?.columns?.bySemType(DG.SEMTYPE.MOLECULE) ?? null;
    column2Input = createColumnInput('Reactant 2', t, molCol,
      'Second reactant molecule column');
    col2Container.replaceChildren(column2Input.root);
  });

  const modeInput = ui.input.choice('Combination Mode', {
    value: 'pairwise',
    items: ['pairwise', 'matrix'],
    nullable: false,
  });
  ui.tooltip.bind(modeInput.input, 'Pairwise: row 1 + row 1, row 2 + row 2, etc. Matrix: all N*M combinations.');

  const saltsInput = ui.input.bool('Remove salts and water', {value: true});
  ui.tooltip.bind(saltsInput.input, 'Remove common salts and water before reacting');

  // ---- Browser ----
  const browser = new ReactionBrowser(rdkit, {
    modeFilter: 'two-component',
    showDetails: true,
    allowAdd: true,
    allowEdit: true,
  });

  // ---- Variables ----
  const variablesDiv = ui.divV([], {style: {paddingTop: '4px'}});
  let currentVariableInputs: DG.InputBase[] = [];

  // ---- Preview ----
  const r = window.devicePixelRatio;
  const previewCanvas = ui.canvas(400 * r, 120 * r);
  previewCanvas.style.width = '400px';
  previewCanvas.style.maxWidth = '100%';
  previewCanvas.style.height = '120px';
  previewCanvas.style.border = '1px solid var(--grey-2)';
  previewCanvas.style.borderRadius = '4px';
  previewCanvas.style.cursor = 'pointer';
  previewCanvas.addEventListener('click', () => openFullscreenPreview(rdkit, browser.selected));

  const exampleGrid = createTwoCompExampleGrid();

  function getResolvedVariables(): Record<string, any> {
    const reaction = browser.selected;
    const vars: Record<string, any> = {};
    if (reaction?.variables) {
      Object.keys(reaction.variables).forEach((key, i) => {
        if (currentVariableInputs[i])
          vars[key] = currentVariableInputs[i].value;
      });
    }
    return vars;
  }

  function getResolvedSmarts(): string | null {
    const reaction = browser.selected;
    if (!reaction) return null;
    const vars = getResolvedVariables();
    return reaction.variables ?
      resolveReactionVariables(reaction.reactionSmarts, vars) :
      reaction.reactionSmarts;
  }

  function updatePreview(reaction: NamedReaction | null) {
    const ctx = previewCanvas.getContext('2d')!;
    ctx.clearRect(0, 0, previewCanvas.width, previewCanvas.height);
    if (!reaction) return;

    let smarts = reaction.reactionSmarts;
    if (reaction.variables) {
      const defaults: Record<string, any> = {};
      for (const [key, varDef] of Object.entries(reaction.variables))
        defaults[key] = varDef.defaultValue;
      smarts = resolveReactionVariables(smarts, defaults);
    }

    try {
      renderReactionToCanvas(rdkit, previewCanvas, smarts, previewCanvas.width, previewCanvas.height);
    } catch {/* preview error */}

    const col1 = column1Input.value;
    const col2 = column2Input.value;
    if (col1 && col2) {
      const examples1 = col1.categories.slice(0, 10).filter((s: string) => !!s);
      const examples2 = col2.categories.slice(0, 10).filter((s: string) => !!s);
      const len = Math.min(examples1.length, examples2.length);
      if (len > 0) {
        runTwoComponentReaction(rdkit, smarts, examples1.slice(0, len), examples2.slice(0, len), {
          mode: 'pairwise',
          removeSaltsAndWater: saltsInput.value,
        }).then((result) => {
          updateTwoCompExampleGrid(exampleGrid, examples1.slice(0, len), examples2.slice(0, len), result.products);
        }).catch(() => {/* silently fail preview */});
      }
    }
  }

  // Wire browser selection → variables + preview
  browser.onSelectionChanged.subscribe((reaction) => {
    variablesDiv.innerHTML = '';
    currentVariableInputs = [];
    if (!reaction) {updatePreview(null); return;}

    if (reaction.variables) {
      for (const [_key, varDef] of Object.entries(reaction.variables)) {
        const inputCtor = varDef.type === 'string' ? ui.input.string :
          varDef.type === 'int' ? ui.input.int : ui.input.float;
        const inp = inputCtor(varDef.name, {
          value: varDef.defaultValue,
          nullable: false,
          onValueChanged: () => updatePreview(browser.selected),
        });
        if (varDef.description)
          ui.tooltip.bind(inp.input, varDef.description);
        currentVariableInputs.push(inp);
        variablesDiv.append(inp.root);
      }
    }
    updatePreview(reaction);
  });

  const newReactionButton = ui.button('+ New Reaction', () => {
    openReactionEditor(rdkit, {
      defaultMode: 'two-component',
      onSave: async () => {await browser.reload();},
    });
  }, 'Create a new custom reaction');

  // ---- Run: add column to existing table / create new table (matrix) ----
  async function runOnExistingTable(): Promise<void> {
    const reaction = browser.selected;
    if (!reaction) {grok.shell.warning('No reaction selected.'); return;}
    const col1 = column1Input.value;
    const col2 = column2Input.value;
    const t1 = table1Input.value;
    if (!col1 || !col2 || !t1) {grok.shell.warning('Please select both reactant columns.'); return;}

    const pg = DG.TaskBarProgressIndicator.create('Running two-component reaction...');
    try {
      const result = await runTwoComponentReaction(rdkit, reaction.reactionSmarts, col1.toList(), col2.toList(), {
        mode: modeInput.value as 'pairwise' | 'matrix',
        removeSaltsAndWater: saltsInput.value,
        showInfo: true,
      });

      if (modeInput.value === 'pairwise') {
        const newCol = t1.columns.addNewString(t1.columns.getUnusedName(`Product(${col1.name}+${col2.name})`));
        newCol.semType = DG.SEMTYPE.MOLECULE;
        newCol.setTag('cell.renderer', 'Molecule');
        newCol.init((i) => result.products[i] ?? '');
      } else {
        const df = DG.DataFrame.create(result.products.length);
        df.name = `${reaction.name} Products`;
        const r1Col = df.columns.addNewString('Reactant 1');
        r1Col.semType = DG.SEMTYPE.MOLECULE;
        r1Col.setTag('cell.renderer', 'Molecule');
        r1Col.init((i) => result.reactants1[i] ?? '');
        const r2Col = df.columns.addNewString('Reactant 2');
        r2Col.semType = DG.SEMTYPE.MOLECULE;
        r2Col.setTag('cell.renderer', 'Molecule');
        r2Col.init((i) => result.reactants2[i] ?? '');
        const pCol = df.columns.addNewString('Product');
        pCol.semType = DG.SEMTYPE.MOLECULE;
        pCol.setTag('cell.renderer', 'Molecule');
        pCol.init((i) => result.products[i] ?? '');
        grok.shell.addTableView(df);
      }
    } catch (e: any) {
      grok.shell.error(`Reaction failed: ${e?.message ?? e}`);
    } finally {
      pg.close();
    }
  }

  // ---- Run: always create new DataFrame ----
  async function runToNewDataFrame(): Promise<DG.DataFrame | null> {
    const reaction = browser.selected;
    if (!reaction) {grok.shell.warning('No reaction selected.'); return null;}
    const col1 = column1Input.value;
    const col2 = column2Input.value;
    if (!col1 || !col2) {grok.shell.warning('Please select both reactant columns.'); return null;}

    const pg = DG.TaskBarProgressIndicator.create('Running two-component reaction...');
    try {
      const result = await runTwoComponentReaction(rdkit, reaction.reactionSmarts, col1.toList(), col2.toList(), {
        mode: modeInput.value as 'pairwise' | 'matrix',
        removeSaltsAndWater: saltsInput.value,
        showInfo: true,
      });

      const df = DG.DataFrame.create(result.products.length);
      df.name = `${reaction.name} Products`;
      const r1Col = df.columns.addNewString('Reactant 1');
      r1Col.semType = DG.SEMTYPE.MOLECULE;
      r1Col.setTag('cell.renderer', 'Molecule');
      r1Col.init((i) => result.reactants1[i] ?? '');
      const r2Col = df.columns.addNewString('Reactant 2');
      r2Col.semType = DG.SEMTYPE.MOLECULE;
      r2Col.setTag('cell.renderer', 'Molecule');
      r2Col.init((i) => result.reactants2[i] ?? '');
      const pCol = df.columns.addNewString('Product');
      pCol.semType = DG.SEMTYPE.MOLECULE;
      pCol.setTag('cell.renderer', 'Molecule');
      pCol.init((i) => result.products[i] ?? '');
      return df;
    } catch (e: any) {
      grok.shell.error(`Reaction failed: ${e?.message ?? e}`);
      return null;
    } finally {
      pg.close();
    }
  }

  return {
    table1Input, col1Container, table2Input, col2Container,
    modeInput, saltsInput,
    browser, newReactionButton,
    previewCanvas, variablesDiv, exampleGrid,
    getColumn1Input: () => column1Input,
    getColumn2Input: () => column2Input,
    getResolvedSmarts,
    runOnExistingTable, runToNewDataFrame,
  };
}

// ==========================================
//  DIALOG: transformation (single column)
// ==========================================

/** Main entry point for the single-column transformation reaction dialog. */
export async function transformationReactionsUI(table?: DG.DataFrame, column?: DG.Column | null): Promise<void> {
  const rdkit = getRdKitModule();
  table ??= grok.shell.t;
  column ??= table?.columns?.bySemType(DG.SEMTYPE.MOLECULE);

  if (!table || !column) {
    grok.shell.warning('No table or molecule column found.');
    return;
  }

  const c = buildTransformationComponents(rdkit, table, column);

  const previewSection = ui.divV([
    ui.h3('Reaction Preview'),
    c.previewCanvas,
    c.variablesDiv,
    c.exampleGrid.root,
  ], {style: {minWidth: '480px', flexShrink: '0'}});

  const leftPanel = ui.divV([
    ui.divH([c.tableInput.root, c.columnContainer], {style: {gap: '8px'}}),
    ui.divH([c.outputNameInput.root, c.saltsInput.root], {style: {gap: '8px'}}),
    c.browser.root,
    c.newReactionButton,
  ], {style: {flex: '1 1 auto', overflow: 'auto', maxWidth: 'min(calc(100vw - 600px), 840px)'}});

  const mainContent = ui.divH([leftPanel, previewSection], {style: {gap: '16px', alignItems: 'stretch'}});

  const dialog = ui.dialog('Run Reaction')
    .add(mainContent);

  dialog.onOK(async () => {await c.runOnExistingTable();});
  dialog.show({resizable: true});
}

// ==========================================
//  DIALOG: two-component (two columns)
// ==========================================

/** Main entry point for the two-component reaction dialog. */
export async function twoComponentReactionUI(table?: DG.DataFrame): Promise<void> {
  const rdkit = getRdKitModule();
  table ??= grok.shell.t;

  if (!table) {
    grok.shell.warning('No table found.');
    return;
  }

  const c = buildTwoComponentComponents(rdkit, table);

  const rightPanel = ui.divV([
    ui.h3('Reaction Preview'),
    c.previewCanvas,
    c.variablesDiv,
    c.exampleGrid.root,
  ], {style: {minWidth: '480px', flexShrink: '0'}});

  const leftPanel = ui.divV([
    ui.divH([c.table1Input.root, c.col1Container], {style: {gap: '8px'}}),
    ui.divH([c.table2Input.root, c.col2Container], {style: {gap: '8px'}}),
    ui.divH([c.modeInput.root, c.saltsInput.root], {style: {gap: '8px'}}),
    c.browser.root,
    c.newReactionButton,
  ], {style: {flex: '1 1 auto', overflow: 'auto', maxWidth: 'min(calc(100vw - 600px), 840px)'}});

  const mainContent = ui.divH([leftPanel, rightPanel], {style: {gap: '16px', alignItems: 'stretch'}});

  const dialog = ui.dialog('Two-Component Reaction')
    .add(mainContent);

  dialog.onOK(async () => {await c.runOnExistingTable();});
  dialog.show({resizable: true});
}

// ==========================================
//  VIEW: transformation (single column)
// ==========================================

/** Returns a DG.ViewBase with splitH layout for transformation reactions.
 *  Left: reaction browser with filters. Right: inputs, preview, example grid, run button. */
export async function transformationReactionsView(
  table?: DG.DataFrame, column?: DG.Column | null,
): Promise<DG.ViewBase> {
  const rdkit = getRdKitModule();
  table ??= grok.shell.t;
  column ??= table?.columns?.bySemType(DG.SEMTYPE.MOLECULE);

  const c = buildTransformationComponents(rdkit, table, column);

  // ---- Left: browser (cards fill available space) ----
  c.browser.root.style.height = '100%';
  c.browser.root.style.display = 'flex';
  c.browser.root.style.flexDirection = 'column';
  const cardsEl = c.browser.root.children[1] as HTMLElement;
  if (cardsEl)
    cardsEl.style.maxHeight = 'none';
    // cardsEl.style.flex = '1 1 auto';

  const leftSide = ui.divV([c.browser.root], {style: {
    height: '100%', overflow: 'hidden',
  }});

  // ---- Right: inputs + preview + grid + run button ----
  c.outputNameInput.input.style.width = '60%';
  const inputsRow = ui.divV([
    ui.divH([c.tableInput.root, c.columnContainer], {style: {gap: '8px', flexWrap: 'wrap'}}),
    ui.divH([c.outputNameInput.root, c.saltsInput.root], {style: {gap: '8px', flexWrap: 'wrap'}}),
  ], {style: {paddingBottom: '8px'}});
  styleViewInputs(inputsRow);

  const runButton = ui.button('Run Reaction', async () => {
    const df = await c.runToNewDataFrame();
    if (df)
      grok.shell.addTableView(df);
  });
  const runRow = ui.divH([runButton], {style: {
    justifyContent: 'flex-end', marginTop: '8px', paddingTop: '8px',
    borderTop: '1px solid var(--grey-2)',
  }});

  const rightSide = ui.divV([
    inputsRow,
    ui.h3('Reaction Preview'),
    c.previewCanvas,
    c.variablesDiv,
    c.exampleGrid.root,
    runRow,
  ], {style: {height: '100%', overflow: 'auto', padding: '8px'}});

  // ---- View ----
  const view = DG.View.create();
  view.name = 'Transformation Reactions';
  const split = ui.splitH([leftSide, rightSide], null, true);
  split.style.width = '100%';
  split.style.height = '100%';
  view.append(split);

  return view;
}

// ==========================================
//  VIEW: two-component (two columns)
// ==========================================

/** Returns a DG.ViewBase with splitH layout for two-component reactions.
 *  Left: reaction browser with filters. Right: inputs, preview, example grid, run button. */
export async function twoComponentReactionsView(table?: DG.DataFrame): Promise<DG.ViewBase> {
  const rdkit = getRdKitModule();
  table ??= grok.shell.t;

  const c = buildTwoComponentComponents(rdkit, table);

  // ---- Left: browser (cards fill available space) ----
  c.browser.root.style.height = '100%';
  c.browser.root.style.display = 'flex';
  c.browser.root.style.flexDirection = 'column';
  const cardsEl = c.browser.root.children[1] as HTMLElement;
  if (cardsEl)
    cardsEl.style.maxHeight = 'none';
    // cardsEl.style.flex = '1 1 auto';

  const leftSide = ui.divV([c.browser.root], {style: {
    height: '100%', overflow: 'hidden',
  }});

  // ---- Right: inputs + preview + grid + run button ----
  const inputsRow = ui.divV([
    ui.divH([c.table1Input.root, c.col1Container], {style: {gap: '8px', flexWrap: 'wrap'}}),
    ui.divH([c.table2Input.root, c.col2Container], {style: {gap: '8px', flexWrap: 'wrap'}}),
    ui.divH([c.modeInput.root, c.saltsInput.root], {style: {gap: '8px', flexWrap: 'wrap'}}),
  ], {style: {paddingBottom: '8px'}});
  styleViewInputs(inputsRow);

  const runButton = ui.button('Run Reaction', async () => {
    const df = await c.runToNewDataFrame();
    if (df)
      grok.shell.addTableView(df);
  });
  const runRow = ui.divH([runButton], {style: {
    justifyContent: 'flex-end', marginTop: '8px', paddingTop: '8px',
    borderTop: '1px solid var(--grey-2)',
  }});

  const rightSide = ui.divV([
    inputsRow,
    ui.h3('Reaction Preview'),
    c.previewCanvas,
    c.variablesDiv,
    c.exampleGrid.root,
    runRow,
  ], {style: {height: '100%', overflow: 'auto', padding: '8px'}});

  // ---- View ----
  const view = DG.View.create();
  view.name = 'Two-Component Reactions';
  const split = ui.splitH([leftSide, rightSide], null, true);
  split.style.width = '100%';
  split.style.height = '100%';
  view.append(split);

  return view;
}

// ==========================================
//  Helpers
// ==========================================

/** Make d4-input-root elements inside a container expand to fill available width. */
function styleViewInputs(container: HTMLElement): void {
  container.querySelectorAll('.d4-input-root').forEach((el) => {
    const inputRoot = el as HTMLElement;
    inputRoot.style.flexGrow = '1';
    inputRoot.style.minWidth = '120px';
  });
}

function createExampleGrid(): ExampleGridInfo {
  const df = DG.DataFrame.create(0);
  const molCol = df.columns.addNew('Input', DG.TYPE.STRING);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.setTag('cell.renderer', 'Molecule');
  const prodCol = df.columns.addNew('Product', DG.TYPE.STRING);
  prodCol.semType = DG.SEMTYPE.MOLECULE;
  prodCol.setTag('cell.renderer', 'Molecule');
  const grid = df.plot.grid();
  grid.root.style.width = '100%';
  grid.root.style.flexGrow = '1';
  grid.props.rowHeight = 100;
  return {root: grid.root, grid, df};
}

function updateExampleGrid(
  gridInfo: ExampleGridInfo,
  inputs: string[], products: string[],
): void {
  const len = Math.min(inputs.length, products.length, 10);
  const df = DG.DataFrame.create(len);
  const molCol = df.columns.addNew('Input', DG.TYPE.STRING);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.setTag('cell.renderer', 'Molecule');
  molCol.init((i) => inputs[i]);
  const prodCol = df.columns.addNew('Product', DG.TYPE.STRING);
  prodCol.semType = DG.SEMTYPE.MOLECULE;
  prodCol.setTag('cell.renderer', 'Molecule');
  prodCol.init((i) => products[i]);

  if (gridInfo.grid) {
    gridInfo.grid.dataFrame = df;
    gridInfo.grid.props.rowHeight = 100;
    const inpGridCol = gridInfo.grid.col('Input');
    const prodGridCol = gridInfo.grid.col('Product');
    if (inpGridCol) inpGridCol.width = 230;
    if (prodGridCol) prodGridCol.width = 230;
  }
}

function createTwoCompExampleGrid(): ExampleGridInfo {
  const df = DG.DataFrame.create(0);
  const r1Col = df.columns.addNew('Reactant 1', DG.TYPE.STRING);
  r1Col.semType = DG.SEMTYPE.MOLECULE;
  r1Col.setTag('cell.renderer', 'Molecule');
  const r2Col = df.columns.addNew('Reactant 2', DG.TYPE.STRING);
  r2Col.semType = DG.SEMTYPE.MOLECULE;
  r2Col.setTag('cell.renderer', 'Molecule');
  const prodCol = df.columns.addNew('Product', DG.TYPE.STRING);
  prodCol.semType = DG.SEMTYPE.MOLECULE;
  prodCol.setTag('cell.renderer', 'Molecule');
  const grid = df.plot.grid();
  grid.root.style.width = '100%';
  grid.root.style.flexGrow = '1';
  grid.props.rowHeight = 100;
  return {root: grid.root, grid, df};
}

function updateTwoCompExampleGrid(
  gridInfo: ExampleGridInfo,
  reactants1: string[], reactants2: string[], products: string[],
): void {
  const len = Math.min(reactants1.length, reactants2.length, products.length, 10);
  const df = DG.DataFrame.create(len);
  const r1Col = df.columns.addNew('Reactant 1', DG.TYPE.STRING);
  r1Col.semType = DG.SEMTYPE.MOLECULE;
  r1Col.setTag('cell.renderer', 'Molecule');
  r1Col.init((i) => reactants1[i]);
  const r2Col = df.columns.addNew('Reactant 2', DG.TYPE.STRING);
  r2Col.semType = DG.SEMTYPE.MOLECULE;
  r2Col.setTag('cell.renderer', 'Molecule');
  r2Col.init((i) => reactants2[i]);
  const prodCol = df.columns.addNew('Product', DG.TYPE.STRING);
  prodCol.semType = DG.SEMTYPE.MOLECULE;
  prodCol.setTag('cell.renderer', 'Molecule');
  prodCol.init((i) => products[i]);

  if (gridInfo.grid) {
    gridInfo.grid.dataFrame = df;
    gridInfo.grid.props.rowHeight = 100;
    const r1GridCol = gridInfo.grid.col('Reactant 1');
    const r2GridCol = gridInfo.grid.col('Reactant 2');
    const prodGridCol = gridInfo.grid.col('Product');
    if (r1GridCol) r1GridCol.width = 160;
    if (r2GridCol) r2GridCol.width = 160;
    if (prodGridCol) prodGridCol.width = 160;
  }
}

/** Opens a full-screen dialog showing the selected reaction rendered large. */
function openFullscreenPreview(rdkit: RDModule, reaction: NamedReaction | null): void {
  if (!reaction) return;

  const r = window.devicePixelRatio;
  const w = Math.min(window.innerWidth - 80, 1200);
  const h = Math.min(window.innerHeight - 160, 600);
  const canvas = ui.canvas(w * r, h * r);
  canvas.style.width = `${w}px`;
  canvas.style.height = `${h}px`;
  canvas.style.backgroundColor = 'white';
  canvas.style.borderRadius = '4px';

  let smarts = reaction.reactionSmarts;
  if (reaction.variables) {
    const defaults: Record<string, any> = {};
    for (const [key, varDef] of Object.entries(reaction.variables))
      defaults[key] = varDef.defaultValue;
    smarts = resolveReactionVariables(smarts, defaults);
  }

  try {
    renderReactionToCanvas(rdkit, canvas, smarts, canvas.width, canvas.height);
  } catch {
    const ctx = canvas.getContext('2d')!;
    ctx.save();
    ctx.scale(r, r);
    ctx.fillStyle = '#cc0000';
    ctx.font = '16px Roboto, sans-serif';
    ctx.fillText('Failed to render reaction preview', 20, h / 2);
    ctx.restore();
  }

  const dialog = ui.dialog(reaction.name)
    .add(canvas);
  dialog.show({resizable: true});
}
