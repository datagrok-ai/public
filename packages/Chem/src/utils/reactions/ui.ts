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
//  TRANSFORMATION REACTION UI (single col)
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

  // ---- Inputs ----
  const tableInput = ui.input.table('Table', {value: table});
  ui.tooltip.bind(tableInput.input, 'Table containing the molecule column to transform');
  const columnContainer = ui.div();
  let columnInput = ui.input.column('Molecules', {
    value: column,
    table: table,
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
  });
  ui.tooltip.bind(columnInput.input, 'Column with molecules to run the reaction on');
  columnContainer.append(columnInput.root);

  const saltsInput = ui.input.bool('Remove salts and water', {value: true});
  ui.tooltip.bind(saltsInput.input, 'Remove common salts and water from molecules before running the reaction');
  const outputNameInput = ui.input.string('Output column', {
    value: table.columns.getUnusedName(`Reacted(${column.name})`),
    nullable: false,
  });
  ui.tooltip.bind(outputNameInput.input, 'Name of the new column that will contain the reaction products');
  outputNameInput.root.style.flexGrow = '1';
  outputNameInput.root.style.minWidth = '0';
  outputNameInput.root.style.maxWidth = '399px';


  tableInput.onChanged.subscribe(() => {
    const t = tableInput.value;
    if (t) {
      const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
      columnInput = ui.input.column('Molecules', {
        value: molCol ?? undefined,
        table: t,
        filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      });
      ui.tooltip.bind(columnInput.input, 'Column with molecules to run the reaction on');
      columnContainer.replaceChildren(columnInput.root);
    }
  });

  // ---- Reaction browser (filtered to transformation mode) ----
  const browser = new ReactionBrowser(rdkit, {
    modeFilter: 'transformation',
    showDetails: true,
    allowAdd: true,
    allowEdit: true,
  });

  // ---- Variables section ----
  const variablesDiv = ui.divV([], {style: {paddingTop: '4px'}});

  // ---- Preview section (right panel) ----
  const r = window.devicePixelRatio;
  const previewCanvas = ui.canvas(400 * r, 100 * r);
  previewCanvas.style.width = '100%';
  previewCanvas.style.height = '100px';
  previewCanvas.style.border = '1px solid var(--grey-2)';
  previewCanvas.style.borderRadius = '4px';
  previewCanvas.style.cursor = 'pointer';
  previewCanvas.addEventListener('click', () => openFullscreenPreview(rdkit, browser.selected));

  const exampleGrid = createExampleGrid();

  const previewSection = ui.divV([
    ui.h3('Reaction Preview'),
    previewCanvas,
    variablesDiv,
    exampleGrid.root,
  ], {style: {minWidth: '480px', flexShrink: '0'}});

  let currentVariableInputs: DG.InputBase[] = [];

  function updatePreview() {
    const reaction = browser.selected;
    if (!reaction) return;

    // Resolve variables
    const resolvedVars: Record<string, any> = {};
    if (reaction.variables) {
      const varKeys = Object.keys(reaction.variables);
      varKeys.forEach((key, i) => {
        if (currentVariableInputs[i])
          resolvedVars[key] = currentVariableInputs[i].value;
      });
    }

    const smarts = reaction.variables ? resolveReactionVariables(reaction.reactionSmarts, resolvedVars) : reaction.reactionSmarts;

    // Render preview
    const ctx = previewCanvas.getContext('2d')!;
    ctx.clearRect(0, 0, previewCanvas.width, previewCanvas.height);
    try {
      renderReactionToCanvas(rdkit, previewCanvas, smarts, previewCanvas.width, previewCanvas.height);
    } catch {/* preview error */}

    // Update example grid
    const col = columnInput.value;
    if (col) {
      const examples = col.categories.slice(0, 10).filter((s) => !!s);
      if (examples.length > 0) {
        runTransformationReaction(rdkit, smarts, examples, {removeSaltsAndWater: saltsInput.value}).then((result) => {
          updateExampleGrid(exampleGrid, examples, result.products);
        }).catch(() => {/* silently fail preview */});
      }
    }
  }

  // React to reaction selection
  browser.onSelectionChanged.subscribe((reaction) => {
    variablesDiv.innerHTML = '';
    currentVariableInputs = [];

    if (!reaction) return;

    // Create variable inputs
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

  // ---- New Reaction button ----
  const newReactionButton = ui.button('+ New Reaction', () => {
    openReactionEditor(rdkit, {
      onSave: async () => {
        await browser.reload();
      },
    });
  }, 'Create a new custom reaction');

  // ---- Layout: left panel (inputs + browser) / right panel (preview) ----
  const leftPanel = ui.divV([
    ui.divH([tableInput.root, columnContainer], {style: {gap: '8px'}}),
    ui.divH([outputNameInput.root, saltsInput.root], {style: {gap: '8px'}}),
    browser.root,
    newReactionButton,
  ], {style: {flex: '1 1 auto', overflow: 'auto', maxWidth: 'min(calc(100vw - 600px), 840px)'}});

  const mainContent = ui.divH([leftPanel, previewSection], {style: {gap: '16px', alignItems: 'stretch'}});

  // ---- Dialog ----
  const dialog = ui.dialog('Run Reaction')
    .add(mainContent);

  dialog.onOK(async () => {
    const reaction = browser.selected;
    if (!reaction) {
      grok.shell.warning('No reaction selected.');
      return;
    }
    const col = columnInput.value;
    const t = tableInput.value;
    if (!col || !t) {
      grok.shell.warning('No table or molecule column selected.');
      return;
    }

    // Resolve variables
    const resolvedVars: Record<string, any> = {};
    if (reaction.variables) {
      Object.keys(reaction.variables).forEach((key, i) => {
        if (currentVariableInputs[i])
          resolvedVars[key] = currentVariableInputs[i].value;
      });
    }

    const smarts = reaction.variables ? resolveReactionVariables(reaction.reactionSmarts, resolvedVars) : reaction.reactionSmarts;

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
  });

  dialog.show({resizable: true, center: true});
}

// ==========================================
//  TWO-COMPONENT REACTION UI (two columns)
// ==========================================

/** Main entry point for the two-component reaction dialog. */
export async function twoComponentReactionUI(table?: DG.DataFrame): Promise<void> {
  const rdkit = getRdKitModule();
  table ??= grok.shell.t;

  if (!table) {
    grok.shell.warning('No table found.');
    return;
  }

  const molColumns = table.columns.toList().filter((c) => c.semType === DG.SEMTYPE.MOLECULE);

  // ---- Inputs ----
  const table1Input = ui.input.table('Table 1', {value: table});
  ui.tooltip.bind(table1Input.input, 'Table containing the first reactant column');
  const col1Container = ui.div();
  let column1Input = ui.input.column('Reactant 1', {
    value: molColumns[0] ?? undefined,
    table: table,
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
  });
  ui.tooltip.bind(column1Input.input, 'First reactant molecule column');
  col1Container.append(column1Input.root);

  const table2Input = ui.input.table('Table 2', {value: table});
  ui.tooltip.bind(table2Input.input, 'Table containing the second reactant column');
  const col2Container = ui.div();
  let column2Input = ui.input.column('Reactant 2', {
    value: molColumns.length > 1 ? molColumns[1] : molColumns[0] ?? undefined,
    table: table,
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
  });
  ui.tooltip.bind(column2Input.input, 'Second reactant molecule column');
  col2Container.append(column2Input.root);

  table1Input.onChanged.subscribe(() => {
    const t = table1Input.value;
    if (t) {
      const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
      column1Input = ui.input.column('Reactant 1', {
        value: molCol ?? undefined,
        table: t,
        filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      });
      ui.tooltip.bind(column1Input.input, 'First reactant molecule column');
      col1Container.replaceChildren(column1Input.root);
    }
  });
  table2Input.onChanged.subscribe(() => {
    const t = table2Input.value;
    if (t) {
      const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
      column2Input = ui.input.column('Reactant 2', {
        value: molCol ?? undefined,
        table: t,
        filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      });
      ui.tooltip.bind(column2Input.input, 'Second reactant molecule column');
      col2Container.replaceChildren(column2Input.root);
    }
  });

  const modeInput = ui.input.choice('Combination Mode', {
    value: 'pairwise',
    items: ['pairwise', 'matrix'],
    nullable: false,
  });
  ui.tooltip.bind(modeInput.input, 'Pairwise: row 1 + row 1, row 2 + row 2, etc. Matrix: all N*M combinations.');

  const saltsInput = ui.input.bool('Remove salts and water', {value: true});
  ui.tooltip.bind(saltsInput.input, 'Remove common salts and water before reacting');

  // ---- Reaction browser (filtered to two-component mode) ----
  const browser = new ReactionBrowser(rdkit, {
    modeFilter: 'two-component',
    showDetails: true,
    allowAdd: true,
    allowEdit: true,
  });

  // ---- Preview (right panel) ----
  const r = window.devicePixelRatio;
  const previewCanvas = ui.canvas(400 * r, 120 * r);
  previewCanvas.style.width = '100%';
  previewCanvas.style.height = '120px';
  previewCanvas.style.border = '1px solid var(--grey-2)';
  previewCanvas.style.borderRadius = '4px';
  previewCanvas.style.cursor = 'pointer';
  previewCanvas.addEventListener('click', () => openFullscreenPreview(rdkit, browser.selected));

  const exampleGrid = createTwoCompExampleGrid();

  function updateTwoCompPreview(reaction: NamedReaction | null) {
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

    // Run example reactions using both columns
    const col1 = column1Input.value;
    const col2 = column2Input.value;
    if (col1 && col2) {
      const examples1 = col1.categories.slice(0, 10).filter((s) => !!s);
      const examples2 = col2.categories.slice(0, 10).filter((s) => !!s);
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

  browser.onSelectionChanged.subscribe((reaction) => updateTwoCompPreview(reaction));

  // ---- New Reaction button ----
  const newReactionButton = ui.button('+ New Reaction', () => {
    openReactionEditor(rdkit, {
      defaultMode: 'two-component',
      onSave: async () => {
        await browser.reload();
      },
    });
  }, 'Create a new custom reaction');

  // ---- Layout: left (inputs + browser) / right (preview) ----
  const leftPanel = ui.divV([
    ui.divH([table1Input.root, col1Container], {style: {gap: '8px'}}),
    ui.divH([table2Input.root, col2Container], {style: {gap: '8px'}}),
    ui.divH([modeInput.root, saltsInput.root], {style: {gap: '8px'}}),
    browser.root,
    newReactionButton,
  ], {style: {flex: '1 1 auto', overflow: 'auto', maxWidth: 'min(calc(100vw - 600px), 840px)'}});

  const rightPanel = ui.divV([
    ui.h3('Reaction Preview'),
    previewCanvas,
    exampleGrid.root,
  ], {style: {minWidth: '480px', flexShrink: '0'}});

  const mainContent = ui.divH([leftPanel, rightPanel], {style: {gap: '16px', alignItems: 'stretch'}});

  // ---- Dialog ----
  const dialog = ui.dialog('Two-Component Reaction')
    .add(mainContent);

  dialog.onOK(async () => {
    const reaction = browser.selected;
    if (!reaction) {
      grok.shell.warning('No reaction selected.');
      return;
    }
    const col1 = column1Input.value;
    const col2 = column2Input.value;
    const t1 = table1Input.value;
    if (!col1 || !col2 || !t1) {
      grok.shell.warning('Please select both reactant columns.');
      return;
    }

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
  });

  dialog.show({resizable: true, center: true});
}

// ==========================================
//  Helpers
// ==========================================

function createExampleGrid(): {root: HTMLElement; grid: DG.Grid | null; df: DG.DataFrame} {
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
  gridInfo: {root: HTMLElement; grid: DG.Grid | null; df: DG.DataFrame},
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

function createTwoCompExampleGrid(): {root: HTMLElement; grid: DG.Grid | null; df: DG.DataFrame} {
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
  gridInfo: {root: HTMLElement; grid: DG.Grid | null; df: DG.DataFrame},
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
  dialog.show({resizable: true, center: true});
}
