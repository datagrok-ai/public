import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {cloneConfig, EnumeratorConfig} from './config';
import {enumerate, OutputRow, PerRoundOverride} from './enumerate';
import {getRdKitModule} from '../chem-common-rdkit';
import {MountedViewerRegistry} from './viewer-mount';
import {
  BuiltInputs, buildInputs, buildResultDataFrame, DataKey, Mode, MODE_LABEL, panelHeader, roundsLabel,
  tabPanel,
} from './enumerator-app';

// Small enough to compute fast, large enough to show a representative mixed sample.
const PREVIEW_TARGET_ROWS = 20;
const PREVIEW_MAX_COMBOS_PER_TEMPLATE = 3;
const PREVIEW_MAX_ROUNDS = 2;

// Shuffles arr in place (Fisher-Yates) and returns it, for pickPreviewSamples' random sampling.
function shuffleInPlace<T>(arr: T[]): T[] {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
  return arr;
}

// Biases the sample toward multi-step routes (70/30 split) since those are the more interesting
// case to eyeball in a quick preview, while still showing some single-step ones.
function pickPreviewSamples(rows: OutputRow[], n: number): OutputRow[] {
  if (rows.length <= n) return shuffleInPlace(rows.slice());
  const stepCount = (r: OutputRow) => r.route ? Math.max(0, r.route.split('>>').length - 1) : 0;
  const multi = shuffleInPlace(rows.filter((r) => stepCount(r) > 1));
  const single = shuffleInPlace(rows.filter((r) => stepCount(r) <= 1));
  const targetMulti = Math.min(multi.length, Math.ceil(n * 0.7));
  const combined = [...multi.slice(0, targetMulti), ...single.slice(0, n - targetMulti)];
  if (combined.length < n) combined.push(...multi.slice(targetMulti), ...single.slice(n - targetMulti));
  return shuffleInPlace(combined.slice(0, n));
}

export interface PreviewPanelDeps {
  getConfig: () => EnumeratorConfig;
  currentMode: () => Mode;
  currentRounds: () => number;
  templatesInput: DG.InputBase<DG.DataFrame | null>;
  bbsInput: DG.InputBase<DG.DataFrame | null>;
  reagentsInput: DG.InputBase<DG.DataFrame | null>;
  exclusionInput: DG.InputBase<DG.DataFrame | null>;
  viewerHost: MountedViewerRegistry;
  buildPerRoundOverrides: (cfg: EnumeratorConfig) => PerRoundOverride[] | undefined;
  overrideCountFor: (overrides: PerRoundOverride[] | undefined, mode: Mode, r: number, key: DataKey) => number | null;
  validate: () => string | null;
}

/** Right-pane "Preview" tab (samples a small, budgeted product set) plus the compact recap card
 * embedded in the left-nav's own Preview breadcrumb pane. Re-renders whenever the user switches to
 * the tab; a stale in-flight run is short-circuited via an incrementing run id (isCancelled). */
export class PreviewPanel {
  readonly panel: HTMLElement;
  private readonly host: HTMLElement;
  private readonly status: HTMLElement;
  private readonly recapHost: HTMLElement;
  private runId = 0;

  constructor(private readonly deps: PreviewPanelDeps) {
    this.recapHost = ui.div([], {style: {fontSize: '12px'}});
    this.host = ui.div([]);
    this.status = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
    const header = panelHeader('Quick preview of a small product sample.', this.status);
    this.panel = tabPanel(header, this.host);
  }

  /** The left-nav recap card host — built once here, embedded by EnumeratorConfigForm's own
   * Preview accordion pane factory alongside its own nav row. */
  buildRecapCard(): HTMLElement {
    return this.recapHost;
  }

  renderRecap(): void {
    this.recapHost.innerHTML = '';
    const tDf = this.deps.templatesInput.value;
    const bDf = this.deps.bbsInput.value;
    const mode = this.deps.currentMode();
    const rounds = this.deps.currentRounds();

    const addRow = (label: string, value: string, indent = false): void => {
      this.recapHost.appendChild(ui.divH([
        ui.divText(label, {style: {color: 'var(--grey-5)', width: indent ? '80px' : '96px',
          flex: '0 0 auto', marginLeft: indent ? '16px' : '0'}}),
        ui.divText(value, {style: {fontWeight: indent ? '400' : '600'}}),
      ], {style: {gap: '8px', padding: '2px 0'}}));
    };

    addRow('Strategy', `${MODE_LABEL[mode]} · ${roundsLabel(rounds)}`);

    // Per-round breakdown only for rounds that actually have a custom subset — every round shows
    // "all N" identically otherwise, which just repeats the total for no benefit.
    const overrides = tDf && bDf ? this.deps.buildPerRoundOverrides(this.deps.getConfig()) : undefined;
    const addComponentRows = (
      label: string, df: DG.DataFrame | null, key: 'templates' | 'buildingBlocks',
    ): void => {
      if (!df) {addRow(label, 'Not selected'); return;}
      addRow(label, `${df.rowCount}`);
      for (let r = 1; r <= rounds; r++) {
        const oc = this.deps.overrideCountFor(overrides, mode, r, key);
        if (oc != null) addRow(`Round ${r}`, `${oc} of ${df.rowCount} (custom subset)`, true);
      }
    };
    addComponentRows('Reactions', tDf, 'templates');
    addComponentRows('Building blocks', bDf, 'buildingBlocks');
  }

  /** Bumps the run id without starting new work — called on tab-away, so an in-flight preview
   * (checked only at its own isCancelled points) stops updating an unwatched tab. */
  cancelPendingRun(): void {
    this.runId++;
  }

  private showInPreview(content: HTMLElement | null): void {
    this.deps.viewerHost.close(this.host);
    this.host.innerHTML = '';
    if (content) this.host.appendChild(content);
  }

  private showMessage(text: string, color = 'var(--grey-5)'): void {
    this.showInPreview(ui.divText(text, {style: {color, padding: '20px', textAlign: 'center'}}));
  }

  async refresh(): Promise<void> {
    const myRunId = ++this.runId;
    const err = this.deps.validate();
    if (err) {
      this.status.textContent = '';
      this.showMessage(`Fix the validation error first: ${err}`);
      return;
    }
    this.status.textContent = 'running preview…';
    this.showMessage('Computing preview…');

    const config = this.deps.getConfig();
    const tDf = this.deps.templatesInput.value!;
    const bDf = this.deps.bbsInput.value!;
    const xDf = this.deps.exclusionInput.value;
    const rDf = this.deps.reagentsInput.value;

    let inputs: BuiltInputs;
    try {
      inputs = buildInputs(config, tDf, bDf, xDf, rDf);
    } catch (e) {
      this.status.textContent = '';
      this.showMessage(e instanceof Error ? e.message : String(e), 'var(--red-3)');
      return;
    }

    let rdkit: ReturnType<typeof getRdKitModule>;
    try {rdkit = await getRdKitModule();} catch (e) {
      this.status.textContent = '';
      this.showMessage(`Could not load RDKit: ${e instanceof Error ? e.message : String(e)}`, 'var(--red-3)');
      return;
    }
    if (myRunId !== this.runId) return;

    const previewConfig = cloneConfig(config);
    previewConfig.enumeration.num_rounds = Math.min(previewConfig.enumeration.num_rounds, PREVIEW_MAX_ROUNDS);
    previewConfig.max_num_combinations_per_template =
      previewConfig.max_num_combinations_per_template < 0 ?
        PREVIEW_MAX_COMBOS_PER_TEMPLATE :
        Math.min(previewConfig.max_num_combinations_per_template, PREVIEW_MAX_COMBOS_PER_TEMPLATE);
    previewConfig.max_num_routes_per_compound = 1;
    previewConfig.keep_building_blocks_in_final_output = false;

    const perRoundOverrides = this.deps.buildPerRoundOverrides(config);
    let rows: OutputRow[] = [];
    try {
      const result = await enumerate({
        rdkit, config: previewConfig, ...inputs, perRoundOverrides,
        isCancelled: () => myRunId !== this.runId,
      });
      rows = result.rows;
    } catch (e) {
      if (myRunId !== this.runId) return;
      this.status.textContent = '';
      this.showMessage(`Preview failed: ${e instanceof Error ? e.message : String(e)}`, 'var(--red-3)');
      return;
    }
    if (myRunId !== this.runId) return;

    if (rows.length === 0) {
      this.status.textContent = '';
      this.showMessage('No products produced within the preview budget. Try relaxing filters or ' +
        'verifying that templates and building blocks are compatible.');
      return;
    }

    const samples = pickPreviewSamples(rows, PREVIEW_TARGET_ROWS);
    const df = buildResultDataFrame(samples, 'Preview');
    // rowHeight 110 (vs the data tabs' 75) fits the extra route-step lines; route isn't the last
    // column, so skip extendLastColumn.
    this.deps.viewerHost.mountDf(this.host, df, false, {rowHeight: 110, extendLastColumn: false});
    this.status.textContent =
      `${samples.length} samples of ${rows.length} preview rows (≤ ${previewConfig.enumeration.num_rounds} ` +
      `rounds, ≤ ${PREVIEW_MAX_COMBOS_PER_TEMPLATE} combos / template)`;
  }
}
