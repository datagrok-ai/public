/* eslint-disable max-len */
import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, PackageFunctions} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {buildCombinationLimitFields, buildProductFilterFields, fixNullableIntStepper} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {tryGetRxn} from './enumerate';
import {MountedViewerRegistry} from './viewer-mount';
import {ChipEl, EnumeratorNav} from './enumerator-nav';
import {
  combinationLimitsChanged, detectChemSemTypes, estimateProductCount, MAX_ROUNDS, Mode, MODE_LABEL,
  productFiltersChangedCount, roundsLabel,
} from './shared';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

async function loadBundledCsv(name: string): Promise<DG.DataFrame | null> {
  try {
    const text = await _package.files.readAsText(name);
    const df = DG.DataFrame.fromCsv(text);
    df.name = name.replace(/\.csv$/i, '');
    await detectChemSemTypes(df);
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
    let done = false;
    const cleanup = (f: File | null): void => {
      if (done) return;
      done = true;
      window.removeEventListener('focus', onFocus);
      input.remove();
      resolve(f);
    };
    input.onchange = () => cleanup(input.files?.[0] ?? null);
    input.addEventListener('cancel', () => cleanup(null));
    // Fallback where 'cancel' is unsupported. The delay must outlast 'change' firing, or the `done`
    // guard drops a real selection that arrives after the timeout.
    const onFocus = (): void => {
      setTimeout(() => cleanup(input.files?.[0] ?? null), 1500);
    };
    if (!('oncancel' in input)) window.addEventListener('focus', onFocus);
    input.click();
  });
}

const isStringCol = (c: DG.Column) => c.type === DG.COLUMN_TYPE.STRING;

function makeColInput(
  label: string, table: DG.DataFrame | null, preferredName: string,
  filter: (c: DG.Column) => boolean, tooltip: string, nullable: boolean,
): DG.InputBase<DG.Column | null> {
  // `nullable` must apply even when `table` is already set, or the platform defaults to required.
  const opts: any = {filter, nullable};
  if (table) {
    opts.table = table;
    const c = table.col(preferredName);
    if (c && filter(c)) opts.value = c;
  }
  const inp = ui.input.column(label, opts);
  inp.setTooltip(tooltip);
  return inp;
}

// Swap the column input's table in-place on parent-table change (preserves DOM identity so
// form layout stays valid), re-selecting the preferred column name if it exists in the new table.
function bindColToTable(
  colInput: DG.InputBase<DG.Column | null>, tableInput: DG.InputBase<DG.DataFrame | null>,
  getPreferredName: () => string, filter: (c: DG.Column) => boolean,
) {
  return tableInput.onChanged.subscribe(() => {
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

// A selected column may still hold the wrong data (e.g. "Name" picked where "SMILES" was detected),
// so sample it with the same sampler Chem's own detectors use. validate() reruns on every tracked
// input's onChanged, so cache per Column object: picking a different column or file yields a new
// Column, which makes the cache self-invalidating.
const columnContentValidCache = new WeakMap<DG.Column, boolean>();
function cachedSampleValid(col: DG.Column, isValid: (s: string) => boolean): boolean {
  let result = columnContentValidCache.get(col);
  if (result === undefined) {
    result = DG.Detector.sampleCategories(col, isValid, 1, 5, 0.8);
    columnContentValidCache.set(col, result);
  }
  return result;
}

function isValidReactionSmarts(s: string): boolean {
  const rxn = tryGetRxn(getRdKitModule(), s);
  const ok = !!rxn;
  rxn?.delete();
  return ok;
}

// validateMolecule returns '' (not null) for a valid molecule — only a non-empty string is an error.
function isValidSmiles(s: string): boolean {
  return !PackageFunctions.validateMolecule(s);
}

export interface RibbonOverrideFlags {
  templatesOverride: boolean;
  bbsOverride: boolean;
  reagentsOverride: boolean;
}

export interface EnumeratorConfigFormDeps {
  view: DG.View;
  viewerHost: MountedViewerRegistry;
  refreshValidation: () => void;
  openAccPaneAndSyncTab: (pane: DG.AccordionPane) => void;
  // Lazy: PreviewPanel and RunControls are constructed after this form, since they need its inputs.
  getPreviewRecapCard: () => HTMLElement;
  getPreviewEnumerateBtnWrap: () => HTMLElement;
  // Late-bound — needs the DataPanel instances. Only decides whether Save YAML warns that per-round
  // subsets aren't included.
  hasAnyPerRoundOverride: () => boolean;
}

/** Owns `config`, every data/quick-config input, YAML load/save, and validation. Composes an
 * `EnumeratorNav` for the layout side (strategy cards, left-nav accordion, ribbon chips) and
 * re-exports its panes/chips, so the orchestrator sees one surface. */
export class EnumeratorConfigForm {
  private config: EnumeratorConfig;

  readonly templatesInput: DG.InputBase<DG.DataFrame | null>;
  readonly bbsInput: DG.InputBase<DG.DataFrame | null>;
  readonly reagentsInput: DG.InputBase<DG.DataFrame | null>;
  readonly exclusionInput: DG.InputBase<DG.DataFrame | null>;
  readonly numRoundsInput: DG.InputBase<number | null>;
  readonly depthFirstInput: DG.InputBase<boolean>;

  private readonly smartsColInput: DG.InputBase<DG.Column | null>;
  private readonly blockingColInput: DG.InputBase<DG.Column | null>;
  private readonly rxnNameColInput: DG.InputBase<DG.Column | null>;
  private readonly bbColInput: DG.InputBase<DG.Column | null>;
  private readonly reagentsColInput: DG.InputBase<DG.Column | null>;
  private readonly exclusionColInput: DG.InputBase<DG.Column | null>;
  private readonly maxComponentsInput: DG.InputBase<number | null>;
  private readonly maxRoutesInput: DG.InputBase<number | null>;

  private combinationLimitFields: ReturnType<typeof buildCombinationLimitFields>;
  private productFilterFields: ReturnType<typeof buildProductFilterFields>;
  // Kept out of view.subs: these two field groups are rebuilt on every YAML load, and view.subs
  // only drains at view close, so routing them through it would pile up dead entries per reload.
  private limitFieldSubs: Subscription[] = [];
  private pushingConfigToInputs = false;

  readonly appInfoIcon: HTMLElement;
  private readonly configInfoIcon: HTMLElement;

  readonly loadYamlBtn: HTMLElement;
  readonly saveYamlBtn: HTMLElement;

  private readonly nav: EnumeratorNav;

  // Pass-throughs onto `nav`.
  readonly accReactionsPane: DG.AccordionPane;
  readonly accBbsPane: DG.AccordionPane;
  readonly accExtrasPane: DG.AccordionPane;
  readonly accCombinePane: DG.AccordionPane;
  readonly accPreviewPane: DG.AccordionPane;
  readonly accPanes: DG.AccordionPane[];
  readonly leftPane: HTMLElement;
  readonly chipReactionsC: ChipEl;
  readonly chipBbsC: ChipEl;
  readonly chipExtrasC: ChipEl;
  readonly chipCombineC: ChipEl;
  readonly cfgEstEl: HTMLElement;

  private constructor(
    private readonly deps: EnumeratorConfigFormDeps,
    templatesDf: DG.DataFrame | null, bbsDf: DG.DataFrame | null, exclusionDf: DG.DataFrame | null,
  ) {
    this.config = cloneConfig(DEFAULT_CONFIG);

    // Built eagerly (syncQuickInputsToConfig needs .inputs before the accordion exists); only their
    // ui.form() wrapper is lazy — see lazyFilterForm in enumerator-nav.ts.
    this.combinationLimitFields = buildCombinationLimitFields(this.config);
    this.productFilterFields = buildProductFilterFields(this.config);

    // ---- DATA inputs ----
    this.templatesInput = ui.input.table('Reaction templates file', {
      value: templatesDf ?? undefined, nullable: false,
      tooltipText: 'Table with reaction SMARTS templates. Pick an open workspace table or upload a CSV via the file icon.',
    });
    this.smartsColInput = makeColInput('Reaction SMARTS column', templatesDf, this.config.enumeration.smarts_col, isStringCol,
      'Column in the reaction templates file that contains the reaction SMARTS strings.', false);

    this.blockingColInput = makeColInput('Blocking SMARTS (optional)', templatesDf,
      this.config.enumeration.reactant_blocking_groups_per_template_column, isStringCol,
      'Optional column whose values are SMARTS patterns (separated by ";" or "|"). Excludes building ' +
      'blocks with functional groups incompatible with this template — a building block matching any ' +
      'of them is skipped for this template only and stays available to all others.', true);

    this.rxnNameColInput = makeColInput('Reaction name (optional)', templatesDf,
      this.config.enumeration.reaction_name_col, isStringCol,
      'Optional column with a friendly name for each reaction template. Surfaces in the output ' +
      '"reaction_name" column.', true);

    this.bbsInput = ui.input.table('Building blocks file', {
      value: bbsDf ?? undefined, nullable: false,
      tooltipText: 'Table with the building-block library (SMILES). Pick a workspace table or upload a CSV.',
    });
    this.bbColInput = makeColInput('SMILES column', bbsDf, this.config.enumeration.bb_smiles_column, isStringCol,
      'Column in the building blocks file that contains SMILES.', false);

    this.reagentsInput = ui.input.table('Reagents file (optional)', {
      value: undefined, nullable: true,
      tooltipText: 'Optional table of reagent SMILES. When set, switches to reagents mode: every ' +
        'round uses exactly one building block (or product of an earlier round) and fills every ' +
        'remaining slot with reagents from this file — produces derivatives of each BB across rounds.',
    });
    const REAGENTS_COL_TOOLTIP = 'Column in the reagents file that contains the reagent SMILES.';
    this.reagentsColInput = makeColInput('Reagent SMILES column', null,
      this.config.enumeration.reagent_smiles_column, isStringCol, REAGENTS_COL_TOOLTIP, true);
    // The platform marks this invalid because it has no table to pick a column from yet, not because
    // it's empty (nullable is true) — disable it and say why, instead of fighting the invalid style.
    const syncReagentsColEnabled = (): void => {
      const hasTable = this.reagentsInput.value != null;
      this.reagentsColInput.enabled = hasTable;
      this.reagentsColInput.setTooltip(hasTable ? REAGENTS_COL_TOOLTIP : 'Select a reagents file first.');
    };
    syncReagentsColEnabled();

    this.exclusionInput = ui.input.table('Exclusion SMARTS file (optional)', {
      value: exclusionDf ?? undefined, nullable: true,
      tooltipText: 'Optional table of SMARTS patterns. Any product matching one of these is rejected.',
    });
    this.exclusionColInput = makeColInput('Exclusion SMARTS column', exclusionDf,
      this.config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol,
      'Column in the exclusion substructures file that contains the SMARTS strings.', true);

    // Re-bind column inputs whenever the parent table changes.
    this.deps.view.subs.push(
      bindColToTable(this.smartsColInput, this.templatesInput, () => this.config.enumeration.smarts_col, isStringCol),
      bindColToTable(this.blockingColInput, this.templatesInput,
        () => this.config.enumeration.reactant_blocking_groups_per_template_column, isStringCol),
      bindColToTable(this.rxnNameColInput, this.templatesInput, () => this.config.enumeration.reaction_name_col, isStringCol),
      bindColToTable(this.bbColInput, this.bbsInput, () => this.config.enumeration.bb_smiles_column, isStringCol),
      bindColToTable(this.reagentsColInput, this.reagentsInput,
        () => this.config.enumeration.reagent_smiles_column, isStringCol),
      bindColToTable(this.exclusionColInput, this.exclusionInput,
        () => this.config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol),
      this.reagentsInput.onChanged.subscribe(syncReagentsColEnabled),
    );

    // nullable: true throughout — without it, clearing the box makes .value NaN rather than null,
    // and NaN slips past every `?? fallback` downstream uncaught.
    this.numRoundsInput = ui.input.int('Number of rounds',
      {value: this.config.enumeration.num_rounds, nullable: true, min: 1, max: MAX_ROUNDS, showPlusMinus: true});
    this.numRoundsInput.setTooltip(
      'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
      `products and (in depth-first mode) reacts each one with original BBs. Increase for deeper ` +
      `libraries (capped at ${MAX_ROUNDS} — a round tab is built for every round, and product counts ` +
      `grow combinatorially with each one).`);
    fixNullableIntStepper(this.numRoundsInput, 1);
    // `min`/`max` above only affect the tooltip/spinner, not validation — add that separately.
    this.numRoundsInput.addValidator((v) => {
      const n = Number(v);
      if (!Number.isFinite(n) || n < 1) return 'Must be at least 1.';
      if (n > MAX_ROUNDS) return `Must be at most ${MAX_ROUNDS} — the round strip shows the first ${MAX_ROUNDS} rounds only.`;
      return null;
    });

    this.depthFirstInput = ui.input.bool('Depth first', {value: this.config.enumeration.depth_first});
    this.depthFirstInput.setTooltip(
      'When checked, each round r > 1 must combine EXACTLY ONE round-(r-1) product with original ' +
      'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
      'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

    this.maxComponentsInput = ui.input.int('Max # components',
      {value: this.config.max_num_components, nullable: true, min: 1, showPlusMinus: true});
    this.maxComponentsInput.setTooltip('Max number of reactant components a template may have.');
    fixNullableIntStepper(this.maxComponentsInput, 1);
    // -1 is the config's "no cap" sentinel; blank means the same thing and reads better.
    this.maxRoutesInput = ui.input.int('Max routes per compound', {
      value: this.config.max_num_routes_per_compound < 0 ? undefined : this.config.max_num_routes_per_compound,
      nullable: true, min: 1, showPlusMinus: true,
    });
    this.maxRoutesInput.setTooltip('Cap on the number of routes saved per product. Leave blank for no cap.');
    // A cap of 0 would blank every product's route/template/reaction_name with no error, so floor
    // the stepper at 1; validate() rejects a typed or loaded 0.
    fixNullableIntStepper(this.maxRoutesInput, 1);

    // Bound to live factories so hovering always reflects current state.
    const mkIcon = (): HTMLElement => {
      const i = ui.iconFA('info-circle');
      i.style.marginLeft = '8px';
      i.style.color = 'var(--blue-2)';
      i.style.cursor = 'help';
      return i;
    };
    this.appInfoIcon = mkIcon();
    this.configInfoIcon = mkIcon();
    ui.tooltip.bind(this.appInfoIcon, () => this.buildAppHelp());
    ui.tooltip.bind(this.configInfoIcon, () => this.buildConfigCard());

    this.syncQuickInputsToConfig();

    // Re-validate on every input change so the Run button stays accurate.
    [this.smartsColInput, this.blockingColInput, this.rxnNameColInput, this.bbColInput, this.reagentsColInput,
      this.exclusionInput, this.exclusionColInput, this.numRoundsInput, this.depthFirstInput,
      this.maxComponentsInput, this.maxRoutesInput,
    ].forEach((inp) => this.deps.view.subs.push(this.wireValidationOne(inp)));
    [...this.combinationLimitFields.inputs, ...this.productFilterFields.inputs].forEach((inp) =>
      this.limitFieldSubs.push(this.wireValidationOne(inp)));
    // Covers a view that closes without any YAML load in between.
    this.deps.viewerHost.onClose(() => this.limitFieldSubs.forEach((s) => s.unsubscribe()));

    this.loadYamlBtn = ui.iconFA('folder-open', async () => {
      const f = await pickFile('.yaml,.yml');
      if (!f) return;
      try {
        const text = await f.text();
        this.config = configFromYaml(text);
        this.deps.viewerHost.withPreservedFilters(
          [this.templatesInput.value, this.bbsInput.value, this.reagentsInput.value].filter((d): d is DG.DataFrame => d != null),
          () => this.syncConfigToQuickInputs());
        this.deps.refreshValidation();
        grok.shell.info(`Loaded config from ${f.name}.`);
      } catch (e) {
        grok.shell.error(`Could not load YAML: ${e instanceof Error ? e.message : String(e)}`);
      }
    }, 'Load a YAML config file from disk and apply it to the form.');

    this.saveYamlBtn = ui.iconFA('arrow-to-bottom', () => {
      this.syncQuickInputsToConfig();
      DG.Utils.download('enumerator-config.yaml', configToYaml(this.config), 'text/yaml');
      let msg = 'Saved settings only — the template/building-block/reagent files aren\'t included; ' +
        'select those again yourself when loading this config.';
      if (this.deps.hasAnyPerRoundOverride())
        msg += ' Per-round subsets ("Subset by selection") aren\'t included either, for the same reason.';

      grok.shell.info(msg);
    }, 'Download the current config as a YAML file.');

    this.nav = new EnumeratorNav({
      view: this.deps.view,
      templatesInput: this.templatesInput, smartsColInput: this.smartsColInput,
      blockingColInput: this.blockingColInput, rxnNameColInput: this.rxnNameColInput,
      bbsInput: this.bbsInput, bbColInput: this.bbColInput,
      reagentsInput: this.reagentsInput, reagentsColInput: this.reagentsColInput,
      exclusionInput: this.exclusionInput, exclusionColInput: this.exclusionColInput,
      numRoundsInput: this.numRoundsInput, maxComponentsInput: this.maxComponentsInput,
      maxRoutesInput: this.maxRoutesInput, depthFirstInput: this.depthFirstInput,
      configInfoIcon: this.configInfoIcon,
      getCombinationLimitInputs: () => this.combinationLimitFields.inputs,
      getProductFilterInputs: () => this.productFilterFields.inputs,
      setAndFire: <T>(input: DG.InputBase<T>, v: T): void => this.setAndFire(input, v),
      openAccPaneAndSyncTab: this.deps.openAccPaneAndSyncTab,
      getPreviewRecapCard: this.deps.getPreviewRecapCard,
      getPreviewEnumerateBtnWrap: this.deps.getPreviewEnumerateBtnWrap,
    });
    this.accReactionsPane = this.nav.accReactionsPane;
    this.accBbsPane = this.nav.accBbsPane;
    this.accExtrasPane = this.nav.accExtrasPane;
    this.accCombinePane = this.nav.accCombinePane;
    this.accPreviewPane = this.nav.accPreviewPane;
    this.accPanes = this.nav.accPanes;
    this.leftPane = this.nav.leftPane;
    this.chipReactionsC = this.nav.chipReactionsC;
    this.chipBbsC = this.nav.chipBbsC;
    this.chipExtrasC = this.nav.chipExtrasC;
    this.chipCombineC = this.nav.chipCombineC;
    this.cfgEstEl = this.nav.cfgEstEl;
  }

  static async create(deps: EnumeratorConfigFormDeps): Promise<EnumeratorConfigForm> {
    const [templatesDf, bbsDf, exclusionDf] = await Promise.all([
      loadBundledCsv(BUNDLED_TEMPLATES), loadBundledCsv(BUNDLED_BBS), loadBundledCsv(BUNDLED_EXCLUSION),
    ]);
    return new EnumeratorConfigForm(deps, templatesDf, bbsDf, exclusionDf);
  }

  getConfig = (): EnumeratorConfig => this.config;

  currentMode = (): Mode => {
    return this.reagentsInput.value != null ? 'reagents' : (this.depthFirstInput.value ? 'depth' : 'breadth');
  };

  /** The raw round count as displayed, not DataPanel's clamped one — shared by the ribbon chip,
   * Strategy summary and Preview recap so they can't drift apart. */
  currentRounds = (): number => {
    return this.numRoundsInput.value ?? this.config.enumeration.num_rounds;
  };

  syncQuickInputsToConfig = (): void => {
    // No-op mid-push: each setAndFire() fires onChanged, which reaches here via validate() before
    // every input has been updated, and would write not-yet-pushed values back into config.
    if (this.pushingConfigToInputs) return;
    const config = this.config;
    config.enumeration.num_rounds = this.numRoundsInput.value ?? config.enumeration.num_rounds;
    config.enumeration.depth_first = !!this.depthFirstInput.value;
    config.max_num_components = this.maxComponentsInput.value ?? config.max_num_components;
    config.max_num_routes_per_compound = this.maxRoutesInput.value ?? -1; // blank == no cap
    // Column inputs hold a Column; persist its name, keeping the previous value if unselected.
    config.enumeration.smarts_col = this.smartsColInput.value?.name ?? config.enumeration.smarts_col;
    config.enumeration.reactant_blocking_groups_per_template_column =
      this.blockingColInput.value?.name ?? config.enumeration.reactant_blocking_groups_per_template_column;
    config.enumeration.reaction_name_col =
      this.rxnNameColInput.value?.name ?? config.enumeration.reaction_name_col;
    config.enumeration.bb_smiles_column = this.bbColInput.value?.name ?? config.enumeration.bb_smiles_column;
    config.enumeration.reagent_smiles_column =
      this.reagentsColInput.value?.name ?? config.enumeration.reagent_smiles_column;
    config.products_specs.exclusion_smarts_products_file_smarts_col =
      this.exclusionColInput.value?.name ?? config.products_specs.exclusion_smarts_products_file_smarts_col;
    this.combinationLimitFields.syncToConfig(config);
    this.productFilterFields.syncToConfig(config);
  };

  // The Dart widget doesn't always re-render its visible text when .value is set via the API
  // rather than typed, so push the value into the DOM element too.
  private setAndFire<T>(input: DG.InputBase<T>, v: T): void {
    input.value = v;
    try {
      const el = input.input as HTMLInputElement | undefined;
      if (el?.tagName === 'INPUT' && el.type !== 'checkbox') {
        const desired = v == null ? '' : String(v);
        if (el.value !== desired) el.value = desired;
      }
    } catch {/* ignore — non-textual inputs (column/table/bool/etc.) */}
    try {input.fireChanged();} catch {/* ignore — older API versions */}
  }

  private syncConfigToQuickInputs(): void {
    this.pushingConfigToInputs = true;
    try {
      const config = this.config;
      // Clamp on load too — a hand-edited/older YAML could carry num_rounds above the UI's max.
      if (config.enumeration.num_rounds > MAX_ROUNDS) config.enumeration.num_rounds = MAX_ROUNDS;
      this.setAndFire(this.numRoundsInput, config.enumeration.num_rounds);
      this.setAndFire(this.depthFirstInput, config.enumeration.depth_first);
      this.setAndFire(this.maxComponentsInput, config.max_num_components);
      this.setAndFire(this.maxRoutesInput, config.max_num_routes_per_compound < 0 ? null : config.max_num_routes_per_compound);
      const tDf = this.templatesInput.value;
      if (tDf) {
        const sc = tDf.col(config.enumeration.smarts_col);
        if (sc) this.setAndFire(this.smartsColInput, sc);
        const bc = tDf.col(config.enumeration.reactant_blocking_groups_per_template_column);
        if (bc) this.setAndFire(this.blockingColInput, bc);
        const rc = tDf.col(config.enumeration.reaction_name_col);
        if (rc) this.setAndFire(this.rxnNameColInput, rc);
      }
      const bDf = this.bbsInput.value;
      if (bDf) {
        const c = bDf.col(config.enumeration.bb_smiles_column);
        if (c) this.setAndFire(this.bbColInput, c);
      }
      const rDf = this.reagentsInput.value;
      if (rDf) {
        const c = rDf.col(config.enumeration.reagent_smiles_column);
        if (c) this.setAndFire(this.reagentsColInput, c);
      }
      const xDf = this.exclusionInput.value;
      if (xDf) {
        const c = xDf.col(config.products_specs.exclusion_smarts_products_file_smarts_col);
        if (c) this.setAndFire(this.exclusionColInput, c);
      }
      // Neither field group has a "set value" hook, so rebuild from the loaded config; the fresh
      // inputs need their own revalidation wiring.
      this.limitFieldSubs.forEach((s) => s.unsubscribe());
      this.limitFieldSubs = [];
      this.combinationLimitFields = buildCombinationLimitFields(config);
      this.productFilterFields = buildProductFilterFields(config);
      this.nav.invalidateLimitForms();
      this.combinationLimitFields.inputs.forEach((inp) => this.limitFieldSubs.push(this.wireValidationOne(inp)));
      this.productFilterFields.inputs.forEach((inp) => this.limitFieldSubs.push(this.wireValidationOne(inp)));
    } finally {
      this.pushingConfigToInputs = false;
    }
  }

  private wireValidationOne(input: DG.InputBase<unknown>): Subscription {
    return input.onChanged.subscribe(() => this.deps.refreshValidation());
  }

  validate = (): string | null => {
    this.syncQuickInputsToConfig();
    // A degenerate but valid DataFrame (empty upload) is truthy, so the row-count checks matter:
    // without them Enumerate stays clickable and silently does nothing.
    const tDf = this.templatesInput.value;
    if (!tDf) return 'Select a reaction templates file.';
    if (tDf.rowCount === 0) return 'Reaction templates file has no rows.';
    if (!this.smartsColInput.value) return 'Select a reaction SMARTS column.';
    if (!cachedSampleValid(this.smartsColInput.value, isValidReactionSmarts))
      return 'Selected column does not contain valid reaction SMARTS — pick a different column.';

    const bDf = this.bbsInput.value;
    if (!bDf) return 'Select a building blocks file.';
    if (bDf.rowCount === 0) return 'Building blocks file has no rows.';
    if (!this.bbColInput.value)
      return 'Building blocks are missing, or the wrong column is selected.';
    if (!cachedSampleValid(this.bbColInput.value, isValidSmiles))
      return 'Selected column does not contain valid SMILES — pick a different column.';

    const rDf = this.reagentsInput.value;
    if (rDf && rDf.rowCount === 0) return 'Reagents file has no rows — clear it or pick a non-empty one.';
    if (rDf && !this.reagentsColInput.value)
      return 'Select a reagent SMILES column or clear the reagents file.';

    const xDf = this.exclusionInput.value;
    if (xDf && xDf.rowCount === 0)
      return 'Exclusion substructures file has no rows — clear it or pick a non-empty one.';
    if (xDf && !this.exclusionColInput.value)
      return 'Select an exclusion substructures column or clear the exclusion file.';

    const rounds = this.numRoundsInput.value ?? 0;
    if (rounds < 1) return 'Number of rounds must be at least 1.';
    if (rounds > MAX_ROUNDS) return `Number of rounds must be at most ${MAX_ROUNDS}.`;

    // Reads the input, not this.config: the sync above falls back to the previous config value on a
    // cleared box, so config alone can never observe "the user just cleared this".
    if ((this.maxComponentsInput.value ?? 0) < 1) return 'Max # components must be at least 1.';
    if (this.config.max_num_routes_per_compound === 0)
      return 'Max routes per compound must be at least 1, or blank for no cap.';
    if (this.config.max_num_combinations_per_template === 0)
      return 'Max combinations per template must be at least 1, or blank for no cap.';

    // The only min+max pair in products_specs. Both active and inverted means every product fails
    // both checks: 0 rows, with nothing indicating the two thresholds contradict.
    const ps = this.config.products_specs;
    if (ps.min_num_carbon_atoms >= 0 && ps.max_num_carbon_atoms >= 0 &&
      ps.min_num_carbon_atoms > ps.max_num_carbon_atoms)
      return 'Min carbon atoms must not exceed Max carbon atoms.';

    return null;
  };

  refreshStrategyCards = (): void => {
    this.nav.refreshStrategyCards(this.currentMode());
  };

  private buildAppHelp(): HTMLElement {
    const card = ui.div([], {style: {fontSize: '12px', maxWidth: '520px', lineHeight: '1.5', padding: '4px 2px'}});
    card.innerHTML = `
      <div style="font-weight: bold; font-size: 13px; margin-bottom: 4px;">Chemical library enumeration</div>
      <p style="margin: 0 0 6px 0;">Generate a product library from reaction SMARTS templates and a set of starting materials. Each round can take products from the previous round and grow them further.</p>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Inputs</div>
      <ul style="margin: 0 0 6px 16px; padding: 0;">
        <li><b>Reaction templates file</b> — reaction SMARTS (LHS&gt;&gt;RHS). Optional blocking SMARTS exclude unwanted matches per template; optional reaction names surface in the output.</li>
        <li><b>Building blocks file</b> — SMILES of starting materials.</li>
        <li><b>Reagents file</b> (optional) — switches to <i>reagents mode</i>: every round uses exactly one BB or earlier-round product, with reagents in the remaining slots. Yields derivatives of each BB across rounds.</li>
        <li><b>Exclusion substructures</b> (optional) — SMARTS patterns; any product matching one is rejected.</li>
      </ul>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Enumeration modes</div>
      <ul style="margin: 0 0 6px 16px; padding: 0;">
        <li><b>Depth-first</b> — round 1 combines BBs; round R extends each round-(R-1) product with original BBs (linear chains, no merging of two complex products).</li>
        <li><b>Breadth-first</b> — each round may combine any products from earlier rounds with BBs (convergent routes possible).</li>
        <li><b>Reagents</b> — active whenever a reagents file is selected; overrides depth/breadth-first.</li>
      </ul>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Tips</div>
      <ul style="margin: 0 0 0 16px; padding: 0;">
        <li>Select rows on the right-pane grids and click <i>Subset by selection</i> to enumerate only a subset.</li>
        <li>Open the <i>Preview</i> tab to sample products at a reduced budget before kicking off a full run.</li>
        <li>The <i>i</i> icon next to <i>Enumeration options</i> shows the full current config.</li>
      </ul>`;
    return card;
  }

  private buildConfigCard(): HTMLElement {
    const config = this.config;
    const en = config.enumeration;
    const ps = config.products_specs;
    const mode = this.currentMode();
    const MODE_DESC = {
      depth: '(linear chain extension)',
      breadth: '(convergent allowed)',
      reagents: '(with reagents in other slots)',
    } as const;
    const modeLabel = `${MODE_LABEL[mode]} ${MODE_DESC[mode]}`;
    const fmtNum = (n: number, hint = 'unlimited') => n < 0 ? hint : String(n);
    const yn = (b: boolean) => b ? 'Yes' : 'No';

    const card = ui.div([], {style: {fontSize: '12px', maxHeight: '500px', maxWidth: '440px', overflow: 'auto', padding: '4px 2px', lineHeight: '1.5'}});

    const sectionTitle = (text: string) =>
      ui.divText(text, {style: {fontWeight: 'bold', marginTop: '8px', marginBottom: '2px', paddingBottom: '2px', borderBottom: '1px solid var(--grey-3)'}});
    const row = (label: string, value: string) => ui.divH([
      ui.divText(label, {style: {color: 'var(--grey-6)'}}),
      ui.divText(value, {style: {color: 'var(--text-color)', textAlign: 'right'}}),
    ], {style: {justifyContent: 'space-between', padding: '1px 0', gap: '12px'}});

    card.appendChild(sectionTitle('Enumeration'));
    card.appendChild(row('Rounds', String(en.num_rounds)));
    card.appendChild(row('Mode', modeLabel));
    card.appendChild(row('Max components', String(config.max_num_components)));
    card.appendChild(row('Max combinations / template', fmtNum(config.max_num_combinations_per_template)));
    card.appendChild(row('Max routes / compound', fmtNum(config.max_num_routes_per_compound)));
    card.appendChild(row('Keep BBs in final output', yn(config.keep_building_blocks_in_final_output)));

    card.appendChild(sectionTitle('Columns'));
    card.appendChild(row('Reaction SMARTS', en.smarts_col));
    card.appendChild(row('Blocking SMARTS', en.reactant_blocking_groups_per_template_column));
    card.appendChild(row('Reaction name', en.reaction_name_col));
    card.appendChild(row('Building block SMILES', en.bb_smiles_column));
    card.appendChild(row('Reagent SMILES', en.reagent_smiles_column));
    card.appendChild(row('Exclusion SMARTS', ps.exclusion_smarts_products_file_smarts_col));

    card.appendChild(sectionTitle('Product filters'));
    card.appendChild(row('Carbons (min – max)',
      `${fmtNum(ps.min_num_carbon_atoms, 'any')} – ${fmtNum(ps.max_num_carbon_atoms, 'any')}`));
    card.appendChild(row('Heavy atoms (max)', fmtNum(ps.max_num_heavy_atoms, 'any')));
    card.appendChild(row('Hetero atoms (max)', fmtNum(ps.max_num_hetero_atoms, 'any')));
    card.appendChild(row('Nitrogen (max)', fmtNum(ps.max_num_nitrogen, 'any')));
    card.appendChild(row('Sulfur (max)', fmtNum(ps.max_num_sulfur, 'any')));
    card.appendChild(row('Oxygen (max)', fmtNum(ps.max_num_oxygen, 'any')));
    card.appendChild(row('Metals (max)', fmtNum(ps.max_num_metals, 'any')));
    card.appendChild(row('Halogens (max)', fmtNum(ps.max_num_halogens, 'any')));
    card.appendChild(row('Aromatic atoms (max)', fmtNum(ps.max_num_aromatic_atoms, 'any')));
    card.appendChild(row('Unsat. non-arom. bonds (max)', fmtNum(ps.max_num_unsaturated_nonaromatic_bonds, 'any')));
    card.appendChild(row('Allowed atoms',
      ps.only_these_atoms_allowed.length ? ps.only_these_atoms_allowed.join(', ') : 'any'));
    card.appendChild(row('Reject radicals', yn(ps.remove_radicals)));
    card.appendChild(row('Strip isotopes', yn(ps.remove_isotope_information)));
    card.appendChild(row('Reject charged species', yn(ps.remove_charged_species)));

    return card;
  }

  /** Computes the chip/subtitle state from this form's inputs and config and hands it to `nav` as
   * plain data. `overrides` comes from the orchestrator, which owns the DataPanel instances. */
  refreshRibbonChips(overrides: RibbonOverrideFlags): void {
    const config = this.config;
    const tDf = this.templatesInput.value; const bDf = this.bbsInput.value; const rDf = this.reagentsInput.value;
    const mode = this.currentMode();
    const MODE_ABBR = {depth: 'DF', breadth: 'BF', reagents: 'RM'} as const;
    const roundsText = roundsLabel(this.currentRounds());
    const n = estimateProductCount(tDf, bDf);
    const combChanged = combinationLimitsChanged(config);
    const prodChangedCount = productFiltersChangedCount(config);
    this.nav.applyRibbonState({
      reactionsText: tDf ? `${tDf.rowCount} reactions` : 'No reaction table',
      reactionsOverride: overrides.templatesOverride,
      reactionsSubtitle: tDf ? `${tDf.rowCount} reactions` : 'No table selected',
      bbsText: bDf ? `${bDf.rowCount} BBs` : 'No BBs table',
      bbsOverride: overrides.bbsOverride,
      bbsSubtitle: bDf ? `${bDf.rowCount} building blocks` : 'No table selected',
      extrasText: rDf ? `${rDf.rowCount} reagents` : 'Reagents (None)',
      extrasOverride: overrides.reagentsOverride,
      extrasSubtitle: rDf ? `${rDf.rowCount} reagents` : 'Optional',
      combineChipText: `Strategy: ${MODE_ABBR[mode]} · ${roundsText}`,
      combineSubtitle: `${MODE_LABEL[mode]} · ${roundsText}`, // pane header already says "How to combine"
      estimateText: n > 0 ? `≈ ${n.toLocaleString('en-US')} products` : '',
      limitsChanged: combChanged, filtersChanged: prodChangedCount > 0,
    });
  }
}
