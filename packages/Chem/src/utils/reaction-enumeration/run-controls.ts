import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {enumerate, EnumerationProgress, PerRoundOverride} from './enumerate';
import {getRdKitModule} from '../chem-common-rdkit';
import {buildInputs, buildResultDataFrame} from './enumerator-app';

export interface RunControlsDeps {
  getConfig: () => EnumeratorConfig;
  templatesInput: DG.InputBase<DG.DataFrame | null>;
  bbsInput: DG.InputBase<DG.DataFrame | null>;
  reagentsInput: DG.InputBase<DG.DataFrame | null>;
  exclusionInput: DG.InputBase<DG.DataFrame | null>;
  validate: () => string | null;
  syncQuickInputsToConfig: () => void;
  buildPerRoundOverrides: (cfg: EnumeratorConfig) => PerRoundOverride[] | undefined;
  // Called once a run finishes (success, error, or cancel) — restores both buttons' disabled state
  // + tooltip via the orchestrator's own mediator, which also refreshes chips/cards/validationDiv.
  refreshValidation: () => void;
}

/** Ribbon Enumerate/Cancel/progress chrome, mirrored by a second Enumerate button at the end of the
 * Preview pane's nav chain. Both buttons share one run: `running` blocks a second concurrent run,
 * `cancelled` is checked by enumerate()'s own isCancelled callback. */
export class RunControls {
  readonly runBtn: HTMLButtonElement;
  readonly cancelBtn: HTMLButtonElement;
  readonly progressLabel: HTMLElement;
  readonly previewEnumerateBtnWrap: HTMLElement;

  private readonly previewEnumerateBtn: HTMLButtonElement;
  private cancelled = false;
  private running = false;
  private runBtnRibbonItem: HTMLElement | null = null;
  private lastValidationMsg: string | null = null;
  // WeakSet survives the platform replacing the ribbon-item node mid-session (e.g. after
  // "Subset by selection") — re-attaches on a new node, no-ops on an unchanged one.
  private readonly ancestorsWithClickListener = new WeakSet<HTMLElement>();

  private static readonly RUN_TOOLTIP_DEFAULT =
    'Run enumeration with the current config and add the result to the workspace.';
  private static readonly RUN_TOOLTIP_RUNNING = 'An enumeration is already running.';

  constructor(private readonly deps: RunControlsDeps) {
    this.progressLabel = ui.divText('', {style: {fontSize: '12px', color: 'var(--grey-5)'}});
    this.runBtn = ui.bigButton('Enumerate', () => this.runWithUi(() => this.runEnumeration()));
    // Mirrors the ribbon's Enumerate button — Preview is the end of the Next-button chain, so it
    // gets its own run action too.
    this.previewEnumerateBtn = ui.button('Enumerate', () => this.runWithUi(() => this.runEnumeration()));
    this.previewEnumerateBtn.classList.add('ui-btn-ok');
    // Disabled buttons get pointer-events:none, so hover/click never reaches them — bind tooltips to a
    // live ancestor instead. runBtn has the ribbon's own '.d4-ribbon-item'; previewEnumerateBtn gets
    // its own wrapper div for the same purpose.
    this.previewEnumerateBtnWrap = ui.div([this.previewEnumerateBtn]);
    this.armDisabledTooltip(this.previewEnumerateBtn, this.previewEnumerateBtnWrap);

    this.cancelBtn = ui.button('Cancel', () => {this.cancelled = true;});
    this.cancelBtn.style.display = 'none';
  }

  setValidation(err: string | null): void {
    this.runBtn.disabled = err != null;
    this.previewEnumerateBtn.disabled = err != null;
    this.bindRunTooltip(err);
  }

  private armDisabledTooltip(btn: HTMLButtonElement, ancestor: HTMLElement): void {
    if (this.ancestorsWithClickListener.has(ancestor)) return;
    this.ancestorsWithClickListener.add(ancestor);
    ancestor.addEventListener('click', (e) => {
      if (btn.disabled && this.lastValidationMsg)
        ui.tooltip.show(this.lastValidationMsg, (e as MouseEvent).clientX, (e as MouseEvent).clientY);
    });
  }

  private bindRunTooltip(msg: string | null): void {
    this.lastValidationMsg = msg;
    const el = this.runBtn.closest<HTMLElement>('.d4-ribbon-item');
    if (el) {
      this.runBtnRibbonItem = el;
      this.armDisabledTooltip(this.runBtn, el);
    }
    ui.tooltip.bind(this.runBtnRibbonItem ?? this.runBtn, msg ?? RunControls.RUN_TOOLTIP_DEFAULT);
    ui.tooltip.bind(this.previewEnumerateBtnWrap, msg ?? RunControls.RUN_TOOLTIP_DEFAULT);
  }

  // Shared run chrome: validate, disable both Run buttons, show Cancel + progress, restore on finish.
  // `running` blocks a second concurrent run — previewEnumerateBtn isn't otherwise disabled while one is in progress.
  private async runWithUi(fn: () => Promise<void>): Promise<void> {
    if (this.running || this.deps.validate() != null) return;
    this.deps.syncQuickInputsToConfig();
    this.running = true;
    this.cancelled = false;
    this.runBtn.disabled = true;
    this.previewEnumerateBtn.disabled = true;
    this.bindRunTooltip(RunControls.RUN_TOOLTIP_RUNNING);
    this.cancelBtn.style.display = '';
    this.progressLabel.textContent = 'Initializing…';
    try {
      await fn();
    } catch (e) {
      grok.shell.error(`Enumeration failed: ${e instanceof Error ? e.message : String(e)}`);
      console.error(e);
    } finally {
      this.running = false;
      this.cancelBtn.style.display = 'none';
      this.progressLabel.textContent = '';
      this.deps.refreshValidation(); // restores both buttons' disabled state + tooltip
    }
  }

  private async runEnumeration(): Promise<void> {
    this.progressLabel.textContent = 'Loading RDKit…';
    const rdkit = await getRdKitModule();
    const config = this.deps.getConfig();
    const tDf = this.deps.templatesInput.value!;
    const bDf = this.deps.bbsInput.value!;
    const xDf = this.deps.exclusionInput.value;
    const rDf = this.deps.reagentsInput.value;
    const inputs = buildInputs(config, tDf, bDf, xDf, rDf);

    const reagentsPart = inputs.reagents.length > 0 ? ` × ${inputs.reagents.length} reagents` : '';
    this.progressLabel.textContent =
      `Running: ${inputs.templates.length} templates × ${inputs.buildingBlocks.length} BBs${reagentsPart} × ` +
      `${config.enumeration.num_rounds} round(s)`;
    const onProgress = (p: EnumerationProgress): void => {
      const combo = p.combosTotal && p.combosTotal > 0 ?
        `, combos ${p.combosDone}/${p.combosTotal}` : '';
      this.progressLabel.textContent =
        `Round ${p.round}/${p.numRounds}, template ${p.templateIndex + 1}/${p.numTemplates}${combo}, ` +
        `products: ${p.productsSoFar}`;
    };

    const perRoundOverrides = this.deps.buildPerRoundOverrides(config);
    if (perRoundOverrides) this.progressLabel.textContent += ' · per-round subsets active';

    const start = performance.now();
    const {rows, warnings} = await enumerate({
      rdkit, config, ...inputs, perRoundOverrides, onProgress, isCancelled: () => this.cancelled,
    });
    const elapsed = ((performance.now() - start) / 1000).toFixed(1);

    if (this.cancelled)
      grok.shell.warning(`Enumeration cancelled. Partial results: ${rows.length} rows.`);
    else
      grok.shell.info(`Enumeration done in ${elapsed}s — ${rows.length} rows.`);

    if (warnings.length > 0) {
      console.warn('Enumeration warnings:', warnings);
      // Surface the actual warning TEXT, not just a count — e.g. a per-step override silently not
      // applying is only visible this way, not as a bare count.
      const preview = warnings.slice(0, 3).join(' | ');
      const more = warnings.length > 3 ? ` (+${warnings.length - 3} more; see console)` : '';
      grok.shell.warning(`${preview}${more}`);
    }
    if (rows.length > 0) grok.shell.addTableView(buildResultDataFrame(rows));
  }
}
