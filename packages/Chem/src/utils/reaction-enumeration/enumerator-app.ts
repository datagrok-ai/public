import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {openConfigDialog} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {enumerate, EnumerationProgress, OutputRow, TemplateInput} from './enumerate';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// Drop-in replacement for the standalone package's rdkit-helper. By the time the app function
// runs, Chem's `init` handler has already initialized the RDKit module.
async function getRdKit() {
  return getRdKitModule();
}

const PREVIEW_TARGET_ROWS = 20;
const PREVIEW_MAX_COMBOS_PER_TEMPLATE = 3;
const PREVIEW_MAX_ROUNDS = 2;
const PREVIEW_DEBOUNCE_MS = 400;

async function loadBundledCsv(name: string): Promise<DG.DataFrame | null> {
  try {
    const text = await _package.files.readAsText(name);
    const df = DG.DataFrame.fromCsv(text);
    df.name = name.replace(/\.csv$/i, '');
    return df;
  } catch (e) {
    console.warn(`Could not load bundled file ${name}: ${e}`);
    return null;
  }
}

function pickFile(accept: string): Promise<File | null> {
  return new Promise((resolve) => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = accept;
    input.style.display = 'none';
    document.body.appendChild(input);
    input.onchange = () => {
      const f = input.files?.[0] ?? null;
      document.body.removeChild(input);
      resolve(f);
    };
    input.click();
  });
}

function downloadText(text: string, filename: string, mime = 'text/plain'): void {
  const blob = new Blob([text], {type: mime});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function getStringColumn(df: DG.DataFrame, name: string): string[] {
  const col = df.col(name);
  if (!col) throw new Error(`Column "${name}" not found in "${df.name}". Available: ${df.columns.names().join(', ')}`);
  const out: string[] = new Array(col.length);
  for (let i = 0; i < col.length; i++) {
    const v = col.get(i);
    out[i] = v == null ? '' : String(v);
  }
  return out;
}

function buildResultDataFrame(rows: OutputRow[], name = 'Enumeration result'): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('product', rows.map((r) => r.product)),
    DG.Column.fromStrings('route', rows.map((r) => r.route)),
    DG.Column.fromStrings('template', rows.map((r) => r.template)),
    DG.Column.fromStrings('reaction_name', rows.map((r) => r.reaction_name)),
    DG.Column.fromInt32Array('round', new Int32Array(rows.map((r) => r.round))),
    DG.Column.fromInt32Array('n_routes', new Int32Array(rows.map((r) => r.n_routes))),
  ]);
  df.name = name;
  df.col('product')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('route')!.semType = 'ChemicalReaction';
  df.col('template')!.semType = 'ChemicalReaction';
  return df;
}

// Shuffle helper.
function shuffleInPlace<T>(arr: T[]): T[] {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
  return arr;
}

function routeStepCount(row: OutputRow): number {
  return row.route ? Math.max(0, row.route.split('>>').length - 1) : 0;
}

// Pick samples for the preview, biased toward multi-step (more interesting visually).
// Strategy: aim for ~70% multi-step when available; fall back to single-step otherwise. Within
// each pool, shuffle randomly so consecutive previews show different examples.
function pickPreviewSamples(rows: OutputRow[], n: number): OutputRow[] {
  if (rows.length <= n) return shuffleInPlace(rows.slice());
  const multi = shuffleInPlace(rows.filter((r) => routeStepCount(r) > 1));
  const single = shuffleInPlace(rows.filter((r) => routeStepCount(r) <= 1));
  const targetMulti = Math.min(multi.length, Math.ceil(n * 0.7));
  const fromMulti = multi.slice(0, targetMulti);
  const remaining = n - fromMulti.length;
  const fromSingle = single.slice(0, remaining);
  // If we still don't have enough (e.g. we asked for 70% multi-step but only had a few), top up
  // from whichever pool still has rows.
  const combined = [...fromMulti, ...fromSingle];
  if (combined.length < n) {
    const leftoverMulti = multi.slice(fromMulti.length);
    const leftoverSingle = single.slice(fromSingle.length);
    combined.push(...leftoverMulti, ...leftoverSingle);
  }
  // Shuffle the final selection so multi-step and single-step rows interleave.
  return shuffleInPlace(combined.slice(0, n));
}

interface BuiltInputs {
  templates: TemplateInput[];
  buildingBlocks: string[];
  exclusionSmarts: string[];
}

const isStringCol = (c: DG.Column) => c.type === DG.COLUMN_TYPE.STRING;

// Build a column input bound to a table. If the preferred column name exists and matches the
// filter, it's pre-selected; otherwise the first matching column is used (Datagrok default).
function makeColInput(
  label: string, table: DG.DataFrame | null, preferredName: string,
  filter: (c: DG.Column) => boolean, tooltip: string,
): DG.InputBase<DG.Column | null> {
  const opts: any = {filter};
  if (table) {
    opts.table = table;
    const c = table.col(preferredName);
    if (c && filter(c)) opts.value = c;
  }
  const inp = ui.input.column(label, opts);
  inp.setTooltip(tooltip);
  return inp;
}

// When the parent table changes, swap the column input's table in-place (preserves DOM identity
// so form layout and onChanged subscribers keep working). Re-selects the preferred column name in
// the new table if it exists.
function bindColToTable(
  colInput: DG.InputBase<DG.Column | null>, tableInput: DG.InputBase<DG.DataFrame | null>,
  getPreferredName: () => string, filter: (c: DG.Column) => boolean,
): void {
  tableInput.onChanged.subscribe(() => {
    const t = tableInput.value;
    if (!t) return;
    ui.input.setColumnInputTable(colInput, t, filter);
    const wanted = getPreferredName();
    if (wanted) {
      const c = t.col(wanted);
      if (c && filter(c)) colInput.value = c;
    }
  });
}

function buildInputs(
  config: EnumeratorConfig, tDf: DG.DataFrame, bDf: DG.DataFrame, xDf: DG.DataFrame | null,
): BuiltInputs {
  const smartsList = getStringColumn(tDf, config.enumeration.smarts_col);
  const blockingCol = config.enumeration.reactant_blocking_groups_per_template_column;
  const rxnNameCol = config.enumeration.reaction_name_col;
  const blockingListRaw = tDf.col(blockingCol) ? getStringColumn(tDf, blockingCol) : null;
  const rxnNameList = tDf.col(rxnNameCol) ? getStringColumn(tDf, rxnNameCol) : null;

  const templates: TemplateInput[] = [];
  for (let i = 0; i < smartsList.length; i++) {
    const smarts = (smartsList[i] ?? '').trim();
    if (!smarts) continue;
    const blockingRaw = blockingListRaw ? blockingListRaw[i] : '';
    const blockingSmartsList = blockingRaw ?
      blockingRaw.split(/[;|]/).map((s) => s.trim()).filter((s) => s.length > 0) : [];
    templates.push({smarts, blockingSmartsList, reactionName: rxnNameList?.[i] ?? ''});
  }

  const buildingBlocks = getStringColumn(bDf, config.enumeration.bb_smiles_column)
    .filter((s) => s.trim().length > 0);

  let exclusionSmarts: string[] = [];
  if (xDf) {
    const excCol = config.products_specs.exclusion_smarts_products_file_smarts_col;
    if (xDf.col(excCol))
      exclusionSmarts = getStringColumn(xDf, excCol).filter((s) => s.trim().length > 0);
  }
  return {templates, buildingBlocks, exclusionSmarts};
}

export async function buildEnumeratorView(): Promise<DG.ViewBase> {
  const view = DG.View.create();
  view.name = 'Chemical Enumerator';
  view.box = true;

  let config: EnumeratorConfig = cloneConfig(DEFAULT_CONFIG);

  const templatesDf = await loadBundledCsv(BUNDLED_TEMPLATES);
  const bbsDf = await loadBundledCsv(BUNDLED_BBS);
  const exclusionDf = await loadBundledCsv(BUNDLED_EXCLUSION);

  // ---- DATA inputs ----
  const templatesInput = ui.input.table('Templates', {
    value: templatesDf ?? undefined, nullable: false,
    tooltipText: 'Table with reaction SMARTS templates. Pick an open workspace table or upload a CSV via the file icon.',
  });
  const smartsColInput = makeColInput('SMARTS column', templatesDf, config.enumeration.smarts_col, isStringCol,
    'Column in the Templates table that contains the reaction SMARTS strings.');

  const blockingColInput = makeColInput('Blocking SMARTS column', templatesDf,
    config.enumeration.reactant_blocking_groups_per_template_column, isStringCol,
    'Optional column whose values are SMARTS patterns (separated by ";" or "|"). A building block ' +
    'matching any blocking SMARTS for a given template is excluded from that template only.');

  const rxnNameColInput = makeColInput('Reaction name column', templatesDf,
    config.enumeration.reaction_name_col, isStringCol,
    'Optional column with a friendly name for each reaction template. Surfaces in the output ' +
    '"reaction_name" column.');

  const bbsInput = ui.input.table('Building blocks', {
    value: bbsDf ?? undefined, nullable: false,
    tooltipText: 'Table with the building-block library (SMILES). Pick a workspace table or upload a CSV.',
  });
  const bbColInput = makeColInput('SMILES column', bbsDf, config.enumeration.bb_smiles_column, isStringCol,
    'Column in the Building blocks table that contains SMILES.');

  const exclusionInput = ui.input.table('Exclusion SMARTS (optional)', {
    value: exclusionDf ?? undefined, nullable: true,
    tooltipText: 'Optional table of SMARTS patterns. Any product matching one of these is rejected.',
  });
  const exclusionColInput = makeColInput('Exclusion SMARTS column', exclusionDf,
    config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol,
    'Column in the Exclusion table that contains the SMARTS strings.');

  // Re-bind column inputs whenever the parent table changes. Preserves the column input identity
  // (so form layout and onChanged subscribers stay valid) and re-selects the preferred column
  // name in the new table when present.
  bindColToTable(smartsColInput, templatesInput, () => config.enumeration.smarts_col, isStringCol);
  bindColToTable(blockingColInput, templatesInput,
    () => config.enumeration.reactant_blocking_groups_per_template_column, isStringCol);
  bindColToTable(rxnNameColInput, templatesInput, () => config.enumeration.reaction_name_col, isStringCol);
  bindColToTable(bbColInput, bbsInput, () => config.enumeration.bb_smiles_column, isStringCol);
  bindColToTable(exclusionColInput, exclusionInput,
    () => config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol);

  // ---- CONFIG inputs ----
  const numRoundsInput = ui.input.int('Number of rounds', {value: config.enumeration.num_rounds, min: 1});
  numRoundsInput.setTooltip(
    'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
    'products and (in depth-first mode) reacts each one with original BBs. Increase for deeper libraries.');

  const depthFirstInput = ui.input.bool('Depth first', {value: config.enumeration.depth_first});
  depthFirstInput.setTooltip(
    'When checked, each round-r > 1 step must combine EXACTLY ONE round-(r-1) product with original ' +
    'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
    'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

  const minCInput = ui.input.int('Min carbon atoms', {value: config.products_specs.min_num_carbon_atoms});
  minCInput.setTooltip('Minimum number of carbon atoms a product must have to pass. Set -1 to disable.');

  const maxCInput = ui.input.int('Max carbon atoms', {value: config.products_specs.max_num_carbon_atoms});
  maxCInput.setTooltip('Maximum number of carbon atoms a product may have. Set -1 to disable.');

  const maxCombosInput = ui.input.int('Max combinations / template',
    {value: config.max_num_combinations_per_template});
  maxCombosInput.setTooltip(
    'Per template, per round: cap on how many reactant combinations are actually run. If the ' +
    'cartesian product of matching BBs across slots exceeds this, the enumerator runs the first ' +
    'N combos and stops. Set -1 to disable.');

  const maxComponentsInput = ui.input.int('Max # components', {value: config.max_num_components, min: 1});
  maxComponentsInput.setTooltip(
    'Templates whose reaction SMARTS have more reactant slots than this are skipped. Default 4 ' +
    'covers all templates in the bundled set.');

  // True while we're pushing config → inputs (syncConfigToQuickInputs). Each setAndFire fires
  // onChanged, which triggers schedulePreview → refreshTooltip → syncQuickInputsToConfig — and
  // that read-back happens MID-loop while only some inputs have been updated, overwriting config
  // with the stale values of the inputs we haven't reached yet (e.g. max # components reverts).
  // The flag short-circuits the read-back during the sync; refreshTooltip just reads the live
  // `config` instead.
  let pushingConfigToInputs = false;

  // ---- Header (i) tooltip with full YAML ----
  const infoIcon = ui.iconFA('info-circle', () => {});
  infoIcon.style.marginLeft = '8px';
  infoIcon.style.color = 'var(--blue-2)';
  const refreshTooltip = () => {
    if (!pushingConfigToInputs) syncQuickInputsToConfig();
    const yaml = configToYaml(config);
    const pre = document.createElement('pre');
    pre.style.fontSize = '11px';
    pre.style.maxHeight = '400px';
    pre.style.maxWidth = '500px';
    pre.style.overflow = 'auto';
    pre.style.margin = '0';
    pre.textContent = yaml;
    ui.tooltip.bind(infoIcon, () => pre);
  };

  const syncQuickInputsToConfig = () => {
    config.enumeration.num_rounds = numRoundsInput.value ?? config.enumeration.num_rounds;
    config.enumeration.depth_first = !!depthFirstInput.value;
    config.products_specs.min_num_carbon_atoms = minCInput.value ?? -1;
    config.products_specs.max_num_carbon_atoms = maxCInput.value ?? -1;
    config.max_num_combinations_per_template = maxCombosInput.value ?? config.max_num_combinations_per_template;
    config.max_num_components = maxComponentsInput.value ?? config.max_num_components;
    // Column inputs hold a Column object; persist its name in the YAML config. If the input has no
    // selection, keep the previous config value so YAML round-trip stays stable.
    config.enumeration.smarts_col = smartsColInput.value?.name ?? config.enumeration.smarts_col;
    config.enumeration.reactant_blocking_groups_per_template_column =
      blockingColInput.value?.name ?? config.enumeration.reactant_blocking_groups_per_template_column;
    config.enumeration.reaction_name_col =
      rxnNameColInput.value?.name ?? config.enumeration.reaction_name_col;
    config.enumeration.bb_smiles_column = bbColInput.value?.name ?? config.enumeration.bb_smiles_column;
    config.products_specs.exclusion_smarts_products_file_smarts_col =
      exclusionColInput.value?.name ?? config.products_specs.exclusion_smarts_products_file_smarts_col;
  };

  // Setting `input.value = X` updates the model but, on int/float/string inputs, does NOT
  // always refresh the rendered <input type="number"> text — the Dart-side widget skips the
  // re-render when the value comes via the API rather than user typing. This is what causes the
  // dialog→main "max # components doesn't update" symptom. The fix is to also push the new value
  // into the underlying HTMLInputElement so its visible text matches the model.
  const setAndFire = <T>(input: DG.InputBase<T>, v: T) => {
    input.value = v;
    try {
      const el = input.input as HTMLInputElement | undefined;
      if (el?.tagName === 'INPUT' && el.type !== 'checkbox') {
        const desired = v == null ? '' : String(v);
        if (el.value !== desired) el.value = desired;
      }
    } catch {/* ignore — non-textual inputs (column/table/bool/etc.) */}
    try {input.fireChanged();} catch {/* ignore — older API versions */}
  };

  const syncConfigToQuickInputs = () => {
    pushingConfigToInputs = true;
    try {
      setAndFire(numRoundsInput, config.enumeration.num_rounds);
      setAndFire(depthFirstInput, config.enumeration.depth_first);
      setAndFire(minCInput, config.products_specs.min_num_carbon_atoms);
      setAndFire(maxCInput, config.products_specs.max_num_carbon_atoms);
      setAndFire(maxCombosInput, config.max_num_combinations_per_template);
      setAndFire(maxComponentsInput, config.max_num_components);
      // For column inputs, look up the column object by name on the currently-selected table.
      const tDf = templatesInput.value;
      if (tDf) {
        const sc = tDf.col(config.enumeration.smarts_col);
        if (sc) setAndFire(smartsColInput, sc);
        const bc = tDf.col(config.enumeration.reactant_blocking_groups_per_template_column);
        if (bc) setAndFire(blockingColInput, bc);
        const rc = tDf.col(config.enumeration.reaction_name_col);
        if (rc) setAndFire(rxnNameColInput, rc);
      }
      const bDf = bbsInput.value;
      if (bDf) {
        const c = bDf.col(config.enumeration.bb_smiles_column);
        if (c) setAndFire(bbColInput, c);
      }
      const xDf = exclusionInput.value;
      if (xDf) {
        const c = xDf.col(config.products_specs.exclusion_smarts_products_file_smarts_col);
        if (c) setAndFire(exclusionColInput, c);
      }
    } finally {
      pushingConfigToInputs = false;
    }
  };

  // ---- Validation ----
  function validate(): string | null {
    syncQuickInputsToConfig();
    const tDf = templatesInput.value;
    if (!tDf) return 'Select a Templates table.';
    if (!smartsColInput.value) return 'Select a SMARTS column.';

    const bDf = bbsInput.value;
    if (!bDf) return 'Select a Building blocks table.';
    if (!bbColInput.value) return 'Select a SMILES column.';

    const xDf = exclusionInput.value;
    if (xDf && !exclusionColInput.value)
      return 'Select an Exclusion SMARTS column or clear the Exclusion table.';

    const rounds = numRoundsInput.value ?? 0;
    if (rounds < 1) return 'Number of rounds must be at least 1.';

    const maxComp = maxComponentsInput.value ?? 0;
    if (maxComp < 1) return 'Max # components must be at least 1.';

    return null;
  }

  // ---- Preview area ----
  const validationDiv = ui.divText('', {style: {color: 'var(--red-3)', fontSize: '12px', marginBottom: '6px', minHeight: '16px'}});
  const previewStatus = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-5)', marginLeft: '8px'}});
  // Preview grid host: flex-fills the remaining vertical space so the Run row stays pinned at the bottom.
  // min-height: 0 is critical inside a flex column, otherwise the child grows to its content height.
  const previewHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', border: '1px solid var(--grey-2)', borderRadius: '4px', overflow: 'hidden', background: 'var(--white)'}});

  let previewRunId = 0;
  let previewTimer: any = null;
  function showInPreview(content: HTMLElement | null): void {
    previewHost.innerHTML = '';
    if (content) previewHost.appendChild(content);
  }

  async function updatePreview(): Promise<void> {
    const myRunId = ++previewRunId;
    const err = validate();
    validationDiv.textContent = err ?? '';
    runBtn.disabled = err != null;
    if (err) {
      previewStatus.textContent = '';
      const msg = ui.divText('Fix the validation error above to see a preview.',
        {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}});
      showInPreview(msg);
      return;
    }
    previewStatus.textContent = 'updating preview…';
    const tDf = templatesInput.value!;
    const bDf = bbsInput.value!;
    const xDf = exclusionInput.value;
    let inputs: BuiltInputs;
    try {
      inputs = buildInputs(config, tDf, bDf, xDf);
    } catch (e) {
      previewStatus.textContent = '';
      validationDiv.textContent = e instanceof Error ? e.message : String(e);
      runBtn.disabled = true;
      showInPreview(null);
      return;
    }

    let rdkit: any;
    try {rdkit = await getRdKit();} catch (e) {
      previewStatus.textContent = '';
      validationDiv.textContent = `Could not load RDKit: ${e instanceof Error ? e.message : String(e)}`;
      runBtn.disabled = true;
      showInPreview(null);
      return;
    }
    if (myRunId !== previewRunId) return;

    const previewConfig = cloneConfig(config);
    previewConfig.enumeration.num_rounds = Math.min(previewConfig.enumeration.num_rounds, PREVIEW_MAX_ROUNDS);
    previewConfig.max_num_combinations_per_template =
      previewConfig.max_num_combinations_per_template < 0 ?
        PREVIEW_MAX_COMBOS_PER_TEMPLATE :
        Math.min(previewConfig.max_num_combinations_per_template, PREVIEW_MAX_COMBOS_PER_TEMPLATE);
    previewConfig.max_num_routes_per_compound = 1;
    previewConfig.keep_building_blocks_in_final_output = false;

    let rows: OutputRow[] = [];
    try {
      const result = await enumerate({
        rdkit, config: previewConfig, ...inputs,
        isCancelled: () => myRunId !== previewRunId,
      });
      rows = result.rows;
    } catch (e) {
      if (myRunId !== previewRunId) return;
      previewStatus.textContent = '';
      validationDiv.textContent = `Preview failed: ${e instanceof Error ? e.message : String(e)}`;
      showInPreview(null);
      return;
    }
    if (myRunId !== previewRunId) return;

    if (rows.length === 0) {
      previewStatus.textContent = '';
      const msg = ui.divText('No products were produced for these inputs in the preview budget. ' +
        'Check that templates and building blocks are compatible, or relax the product filters.',
      {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}});
      showInPreview(msg);
      return;
    }

    const samples = pickPreviewSamples(rows, PREVIEW_TARGET_ROWS);
    const df = buildResultDataFrame(samples, 'Preview');
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    // Make rows tall enough that the reaction renderer has room for multi-step images.
    try {
      grid.props.rowHeight = 110;
      const colByName = (n: string) => grid.col(n);
      const setW = (n: string, w: number) => {const c = colByName(n); if (c) c.width = w;};
      setW('product', 180);
      setW('route', 460);
      setW('reaction_name', 140);
      setW('round', 56);
      setW('n_routes', 70);
      // 'template' column shows the SMARTS of the last step's template — useful for debugging
      // a surprising product; set a fixed width so the route column stays readable.
      setW('template', 280);
    } catch (e) {
      console.warn('Preview grid styling failed:', e);
    }
    showInPreview(grid.root);
    previewStatus.textContent =
      `${samples.length} samples (of ${rows.length} preview rows; rounds=1..${previewConfig.enumeration.num_rounds}, combos≤${PREVIEW_MAX_COMBOS_PER_TEMPLATE}/template)`;
  }

  function schedulePreview(): void {
    refreshTooltip();
    if (previewTimer) clearTimeout(previewTimer);
    previewTimer = setTimeout(() => {previewTimer = null; updatePreview();}, PREVIEW_DEBOUNCE_MS);
  }
  refreshTooltip();
  const wireOnChange = <T>(input: DG.InputBase<T>) => input.onChanged.subscribe(() => schedulePreview());
  wireOnChange(templatesInput);
  wireOnChange(smartsColInput);
  wireOnChange(blockingColInput);
  wireOnChange(rxnNameColInput);
  wireOnChange(bbsInput);
  wireOnChange(bbColInput);
  wireOnChange(exclusionInput);
  wireOnChange(exclusionColInput);
  wireOnChange(numRoundsInput);
  wireOnChange(depthFirstInput);
  wireOnChange(minCInput);
  wireOnChange(maxCInput);
  wireOnChange(maxCombosInput);
  wireOnChange(maxComponentsInput);

  // ---- Buttons ----
  const editConfigBtn = ui.button('Edit full config…', async () => {
    syncQuickInputsToConfig();
    const updated = await openConfigDialog(config);
    if (updated) {config = updated; syncConfigToQuickInputs(); refreshTooltip(); schedulePreview();}
  });
  editConfigBtn.title = 'Open the full config form (every YAML field).';

  const loadYamlBtn = ui.button('Load YAML…', async () => {
    const f = await pickFile('.yaml,.yml');
    if (!f) return;
    try {
      const text = await f.text();
      config = configFromYaml(text);
      syncConfigToQuickInputs();
      refreshTooltip();
      schedulePreview();
      grok.shell.info(`Loaded config from ${f.name}.`);
    } catch (e) {
      grok.shell.error(`Could not load YAML: ${e instanceof Error ? e.message : String(e)}`);
    }
  });
  loadYamlBtn.title = 'Load a YAML config file from disk and apply it to the form.';

  const saveYamlBtn = ui.button('Save YAML', () => {
    syncQuickInputsToConfig();
    downloadText(configToYaml(config), 'enumerator-config.yaml', 'text/yaml');
  });
  saveYamlBtn.title = 'Download the current config as a YAML file.';

  // ---- Run / Cancel ----
  const progressLabel = ui.divText('', {style: {fontSize: '12px', color: 'var(--grey-5)'}});
  let cancelled = false;
  const runBtn = ui.bigButton('Run enumeration', async () => {
    if (validate() != null) return;
    syncQuickInputsToConfig();
    cancelled = false;
    runBtn.disabled = true;
    cancelBtn.style.display = '';
    progressLabel.textContent = 'Initializing…';
    try {
      await runFullEnumeration();
    } catch (e) {
      grok.shell.error(`Enumeration failed: ${e instanceof Error ? e.message : String(e)}`);
      console.error(e);
    } finally {
      runBtn.disabled = false;
      cancelBtn.style.display = 'none';
      progressLabel.textContent = '';
    }
  });
  runBtn.title = 'Run the full enumeration with the current config and add the result to the workspace.';

  const cancelBtn = ui.button('Cancel', () => {cancelled = true;});
  cancelBtn.style.display = 'none';

  async function runFullEnumeration(): Promise<void> {
    progressLabel.textContent = 'Loading RDKit…';
    const rdkit = await getRdKit();
    const tDf = templatesInput.value!;
    const bDf = bbsInput.value!;
    const xDf = exclusionInput.value;
    const inputs = buildInputs(config, tDf, bDf, xDf);

    progressLabel.textContent =
      `Running: ${inputs.templates.length} templates × ${inputs.buildingBlocks.length} BBs × ${config.enumeration.num_rounds} rounds`;
    const onProgress = (p: EnumerationProgress) => {
      const combo = p.combosTotal && p.combosTotal > 0 ?
        `, combos ${p.combosDone}/${p.combosTotal}` : '';
      progressLabel.textContent =
        `Round ${p.round}/${p.numRounds}, template ${p.templateIndex + 1}/${p.numTemplates}${combo}, products: ${p.productsSoFar}`;
    };

    const start = performance.now();
    const {rows, warnings} = await enumerate({
      rdkit, config, ...inputs, onProgress, isCancelled: () => cancelled,
    });
    const elapsed = ((performance.now() - start) / 1000).toFixed(1);

    if (cancelled)
      grok.shell.warning(`Enumeration cancelled. Partial results: ${rows.length} rows.`);
    else
      grok.shell.info(`Enumeration done in ${elapsed}s — ${rows.length} rows.`);
    if (warnings.length > 0) {
      console.warn('Enumeration warnings:', warnings);
      grok.shell.warning(`${warnings.length} warning(s); see console for details.`);
    }
    if (rows.length > 0) grok.shell.addTableView(buildResultDataFrame(rows));
  }

  // ---- Layout ----
  // Sized so every section is visible without page scroll: header + data/config row + action bar
  // are at their natural height; the preview grid (flex: 1) absorbs whatever vertical space is left;
  // run row is pinned at the bottom.
  const headerRow = ui.divH([
    ui.divText('Chemical library enumeration', {style: {fontSize: '16px', fontWeight: 'bold'}}),
    infoIcon,
  ], {style: {alignItems: 'center', flex: '0 0 auto'}});

  const dataSection = ui.divV([
    ui.h2('Data'),
    ui.form([
      templatesInput, smartsColInput, blockingColInput, rxnNameColInput,
      bbsInput, bbColInput,
      exclusionInput, exclusionColInput,
    ]),
  ], {style: {flex: '1 1 0', minWidth: '320px'}});

  const configSection = ui.divV([
    ui.h2('Config'),
    ui.form([numRoundsInput, depthFirstInput, maxComponentsInput, minCInput, maxCInput, maxCombosInput]),
  ], {style: {flex: '1 1 0', minWidth: '280px'}});

  const sectionsRow = ui.divH([dataSection, configSection],
    {style: {gap: '24px', alignItems: 'flex-start', flex: '0 0 auto'}});

  // Combined action bar: config buttons on the left, preview header + status on the right, on one line.
  const previewLabel = ui.divText('Preview — sample products',
    {style: {fontWeight: 'bold', fontSize: '13px'}});
  const spacer = ui.div([], {style: {flex: '1 1 auto'}});
  const actionBar = ui.divH(
    [editConfigBtn, loadYamlBtn, saveYamlBtn, spacer, previewLabel, previewStatus],
    {style: {gap: '8px', alignItems: 'center', marginTop: '8px', marginBottom: '6px', flex: '0 0 auto'}});

  // Validation message — only takes height when there's an error; flex: 0 0 auto so it never grows.
  validationDiv.style.flex = '0 0 auto';

  const runRow = ui.divH([runBtn, cancelBtn, progressLabel],
    {style: {gap: '12px', alignItems: 'center', marginTop: '12px', flex: '0 0 auto'}});

  const root = ui.divV([
    headerRow,
    sectionsRow,
    actionBar,
    validationDiv,
    previewHost,
    runRow,
  ], {style: {padding: '16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'}});

  view.append(root);

  // Kick off the initial preview after the view is in the DOM.
  schedulePreview();
  return view;
}
