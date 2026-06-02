/*
 * Wizard shell for the Consensus Pharmacophore pipeline.
 *
 * Builds a DG.ViewBase with this layout:
 *   .cp-wizard-shell (flex column)
 *     header   — family legend + caption (32 px)
 *     ribbon   — Run step / Run pipeline / Save / Export / Reset (30 px)
 *     [stale banner — visible when downstream steps are stale]
 *     rail     — 5 numbered steps with status icons (48 px)
 *     content  — 3-column flex row (fills remaining height)
 *       left   — step-specific panel (460 px)
 *       center — persistent Mol* viewer (flex:1) — MOUNTED ONCE, NEVER REMOVED
 *       right  — PDB detail panel (320 px, visible only on Step 4)
 *     footer   — Prev / status / Next (48 px)
 *
 * Mol* persistence is the key design constraint: the Datagrok Biostructure
 * viewer must be attached to a stable DOM node so that step navigation does
 * not tear it down. `contentCenter` is that node; it is never emptied during
 * step navigation. Only the orchestrator's `refreshDockedViewer` clears it
 * (during full-rebuild for race-free state, per the comment in
 * `orchestrator.ts:788-795`).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ConsensusPharmacophoreApp} from './orchestrator';
import {DEFAULT_OPTIONS, PipelineOptions} from './orchestrator-types';
import {WizardState, StepId, StepStatus, STEP_COUNT, STEP_LABELS} from './wizard-state';
import {
  StepContext, buildStep1, buildStep2, buildStep3, buildStep4, buildStep5,
} from './wizard-steps';
import {FAMILY_CODES, FAMILY_MAP} from './family-map';

// ---------------------------------------------------------------------------
// Inline SVG icons used in the rail step circles
// ---------------------------------------------------------------------------

const CLOCK_SVG =
  '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" ' +
  'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">' +
  '<circle cx="12" cy="12" r="9"></circle><polyline points="12 7 12 12 16 14"></polyline></svg>';

const CHECK_SVG =
  '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" ' +
  'stroke-width="3" stroke-linecap="round" stroke-linejoin="round">' +
  '<polyline points="20 6 9 17 4 12"></polyline></svg>';

const XMARK_SVG =
  '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" ' +
  'stroke-width="3" stroke-linecap="round" stroke-linejoin="round">' +
  '<line x1="18" y1="6" x2="6" y2="18"></line>' +
  '<line x1="6" y1="6" x2="18" y2="18"></line></svg>';

// ---------------------------------------------------------------------------
// WizardShell
// ---------------------------------------------------------------------------

export class WizardShell {
  readonly view: DG.ViewBase;
  readonly state = new WizardState();
  private readonly app: ConsensusPharmacophoreApp;

  private activeStep: StepId = StepId.Fetch;
  /** True while a step (or the full pipeline) is mid-run. A run OWNS the Mol*
   *  viewer until it completes — concurrent navigation / isolation would kick
   *  off a SECOND render pipeline on the same viewer and corrupt it (stray
   *  overlay structures, mis-framed camera). So navigation + isolation are
   *  blocked while this is true. The footer Back/Next are already greyed via
   *  setRunButtonsEnabled, but the RAIL clicks and row-click isolation bypass
   *  that — hence this explicit guard in goToStep + onIsolatePdb. */
  private isStepRunning = false;
  private currentPdbIds: string[] = [];
  private currentOptions: PipelineOptions = {...DEFAULT_OPTIONS};
  private demoLoaded = false;
  /** Snapshot of PDB IDs that produced the last `done` results — for the
   *  "Revert" link in the stale banner. */
  private lastCommittedPdbIds: string[] | null = null;
  private staleBannerVisible = false;
  /** PDB ID currently isolated in the Mol* viewer (Step 2 row click).
   *  null = all PDBs visible (default). Reset on step change. */
  private isolatedPdb: string | null = null;
  /** Which step's visualization the Mol* viewer currently shows. Drives
   *  `syncViewerToStep` so navigating between steps re-renders the viewer to
   *  match the active step (Step 2 = aligned, Step 3 = pocket dots, etc.),
   *  instead of leaving a stale render from a previously-visited step.
   *  null = placeholder / nothing rendered yet. */
  private lastRenderedStep: StepId | null = null;
  /** Serializes viewer re-renders triggered by navigation. Rapid Back/Next
   *  clicks would otherwise overlap async renders — one navigation's
   *  `clearViewer()` could nil `this.viewer` mid-render of another, throwing
   *  "Cannot read properties of undefined (reading 'viewer')". Chaining the
   *  renders guarantees they run one-at-a-time. */
  private viewerSyncChain: Promise<void> = Promise.resolve();

  // Persistent DOM regions
  private readonly contentLeft: HTMLElement;
  private readonly contentCenter: HTMLElement;
  private readonly contentRight: HTMLElement;
  private readonly molstarPlaceholder: HTMLElement;
  /** The primary text line inside the placeholder, updated by
   *  showMolstarPlaceholder so it's context-appropriate (Step 1 before vs
   *  after run). Assigned in buildMolstarPlaceholder. */
  private molstarPlaceholderText!: HTMLElement;
  private readonly railEl: HTMLElement;
  private readonly railSteps: HTMLElement[] = [];
  private readonly railCircles: HTMLElement[] = [];
  private readonly railSublabels: HTMLElement[] = [];
  private readonly railConnectors: HTMLElement[] = [];
  private staleBannerEl?: HTMLElement;
  private readonly staleBannerSlot: HTMLElement;
  private readonly captionEl: HTMLElement;

  // Footer controls (kept as refs so we can enable/disable)
  private readonly footerPrev: HTMLButtonElement;
  private readonly footerNext: HTMLButtonElement;
  private readonly footerStatus: HTMLElement;
  private readonly footerRunStep: HTMLElement;
  private readonly footerRunPipeline: HTMLElement;
  private readonly footerReset: HTMLElement;
  private readonly footerStepIndicator: HTMLElement;

  constructor(seedPdbIds?: string[]) {
    this.app = new ConsensusPharmacophoreApp();
    if (seedPdbIds) {
      this.currentPdbIds = [...seedPdbIds];
      this.demoLoaded = true;
    }

    // Build DOM tree (sub-builders below).
    const header = this.buildHeader();
    this.staleBannerSlot = ui.div([]);
    this.railEl = this.buildRail();
    this.captionEl = (header.querySelector('.cp-wizard-header-caption') as HTMLElement) ??
      ui.divText('');

    // Content area — 3 columns
    this.contentLeft = ui.div([], 'cp-wizard-content-left');
    this.contentCenter = ui.div([], 'cp-wizard-content-center');
    this.contentRight = ui.div([], 'cp-wizard-content-right');
    this.molstarPlaceholder = this.buildMolstarPlaceholder();
    this.contentCenter.append(this.molstarPlaceholder);

    const contentArea = ui.divH(
      [this.contentLeft, this.contentCenter, this.contentRight],
      'cp-wizard-content-area');

    // Footer (hosts both navigation + all action buttons)
    const footerEls = this.buildFooter();
    this.footerPrev = footerEls.prev;
    this.footerNext = footerEls.next;
    this.footerStatus = footerEls.status;
    this.footerRunStep = footerEls.runStep;
    this.footerRunPipeline = footerEls.runPipeline;
    this.footerReset = footerEls.reset;
    this.footerStepIndicator = footerEls.stepIndicator;

    const root = ui.divV(
      [header, this.staleBannerSlot, this.railEl, contentArea, footerEls.root],
      'cp-wizard-shell');
    root.style.height = '100%';

    // Create the DG view from this root.
    this.view = DG.View.fromRoot(root);
    this.view.name = 'Consensus Pharmacophore';

    // Wire app's Mol*/detail mount points into our DOM regions.
    this.app.attachMolstarHost(this.contentCenter);
    this.app.attachDetailHost(this.contentRight);
    // Wire the Step 4 detail-panel per-interaction checkboxes to our option
    // state. Toggling one updates consensusExcludedInteractions (no stale —
    // re-running Step 5 is what applies it) and refreshes the left panel's
    // live count. The detail-panel checkbox itself is self-visual.
    this.app.attachConsensusSelectionHooks({
      isInteractionExcluded: (key) =>
        (this.currentOptions.consensusExcludedInteractions ?? []).includes(key),
      toggleInteraction: (key) => {
        const set = new Set(this.currentOptions.consensusExcludedInteractions ?? []);
        if (set.has(key)) set.delete(key); else set.add(key);
        this.currentOptions = {
          ...this.currentOptions, consensusExcludedInteractions: [...set]};
        this.showStepContent(this.activeStep);
      },
      isPdbExcluded: (pdbId) =>
        (this.currentOptions.consensusExcludedPdbs ?? [])
          .some((p) => p.toUpperCase() === pdbId.toUpperCase()),
      isInputPdbExcluded: (pdbId) =>
        (this.currentOptions.excludedInputPdbs ?? [])
          .some((p) => p.toUpperCase() === pdbId.toUpperCase()),
    });

    // Subscribe to state changes so the rail repaints when steps transition.
    // If the changed step happens to be the active one, also re-evaluate the
    // Next button enable state (it depends on the active step's status).
    this.state.changes.subscribe(({id, status}) => {
      this.applyStatus(id, status);
      if (id === this.activeStep) this.applyNextButtonGuard();
    });

    // Initial paint.
    this.goToStep(StepId.Fetch);
    // Paint all rail nodes once (state.changes emissions during init might
    // happen before the rail elements were built — repaint defensively).
    for (let s = 0; s < STEP_COUNT; s++)
      this.applyStatus(s as StepId, this.state.get(s as StepId));
    // Apply the PDB-IDs guard so Run step + Run full pipeline start disabled
    // until the user pastes IDs or clicks the demo button.
    this.applyPdbIdsGuard();
  }

  /** Public accessor — called from package.ts to add the view to the shell. */
  getView(): DG.ViewBase {
    return this.view;
  }

  // -------------------------------------------------------------------------
  // Builders — header, ribbon, rail, footer, placeholder
  // -------------------------------------------------------------------------

  private buildHeader(): HTMLElement {
    const title = ui.divText('Consensus Pharmacophore', 'cp-wizard-header-title');
    const caption = ui.divText(
      'Click Run step (or Run full pipeline) to render the consensus pharmacophore.',
      'cp-wizard-header-caption');

    // Family legend chips
    const chipRow = ui.divH([], 'cp-legend-row');
    chipRow.style.flexWrap = 'wrap';
    chipRow.style.gap = '8px';
    chipRow.style.padding = '0';
    chipRow.style.alignItems = 'center';
    for (const code of FAMILY_CODES) {
      const fam = FAMILY_MAP[code];
      const item = ui.divH([], 'cp-legend-item');
      item.style.gap = '4px';
      item.style.alignItems = 'center';
      item.style.fontSize = '11px';
      const dot = ui.div();
      dot.style.width = '10px';
      dot.style.height = '10px';
      dot.style.borderRadius = '50%';
      dot.style.background = fam.hexColor;
      dot.style.border = '1px solid rgba(0,0,0,0.2)';
      const label = ui.divText(`${code} · ${fam.name}`);
      label.style.color = '#444';
      item.append(dot, label);
      item.title = `Family ${code} (${fam.name})`;
      chipRow.append(item);
    }

    const titleRow = ui.divH([title, caption]);
    titleRow.style.alignItems = 'baseline';
    titleRow.style.flexWrap = 'wrap';

    return ui.divV([titleRow, chipRow], 'cp-wizard-header');
  }

  private buildRail(): HTMLElement {
    const rail = ui.divH([], 'cp-wizard-rail');

    for (let i = 0; i < STEP_COUNT; i++) {
      const circleWrap = ui.div([], 'cp-wizard-step-circle-wrap');
      const circle = ui.div([], 'cp-wizard-step-circle');
      circleWrap.append(circle);

      const labelEl = ui.divText(`${i + 1}. ${STEP_LABELS[i]}`, 'cp-wizard-step-label');
      const sublabel = ui.divText('Pending', 'cp-wizard-step-sublabel');
      const textWrap = ui.divV([labelEl, sublabel], 'cp-wizard-step-text');

      const step = ui.divH([circleWrap, textWrap], 'cp-wizard-step cp-step-pending');
      step.dataset['step'] = String(i);
      step.addEventListener('click', () => {
        const s = this.state.get(i as StepId);
        // Clicking is only allowed on done/stale steps (you can't skip ahead).
        if (s === 'done' || s === 'stale' || s === 'error') this.goToStep(i as StepId);
        else if (i === this.activeStep) {/* clicking active is a no-op */}
        else grok.shell.info('Complete previous steps first.');
      });
      rail.append(step);

      this.railSteps.push(step);
      this.railCircles.push(circle);
      this.railSublabels.push(sublabel);

      if (i < STEP_COUNT - 1) {
        const connector = ui.div([], 'cp-wizard-step-connector');
        rail.append(connector);
        this.railConnectors.push(connector);
      }
    }
    return rail;
  }

  private buildFooter(): {
    root: HTMLElement;
    prev: HTMLButtonElement;
    next: HTMLButtonElement;
    status: HTMLElement;
    runStep: HTMLElement;
    runPipeline: HTMLElement;
    reset: HTMLElement;
    stepIndicator: HTMLElement;
  } {
    // Click-guard: if the button has the disabled class (set by
    // applyPdbIdsGuard or setRunButtonsEnabled), short-circuit the handler.
    // We can't use `pointer-events: none` in the CSS because that also
    // suppresses the native `title` tooltip on hover.
    const guarded = (handler: () => void) => (el: HTMLElement) => () => {
      if (el.classList.contains('cp-wizard-footer-btn-disabled')) return;
      handler();
    };
    const mkSecondary = (label: string, title: string, handler: () => void): HTMLElement => {
      const btn = ui.div([], 'cp-wizard-footer-btn');
      btn.textContent = label;
      btn.title = title;
      btn.addEventListener('click', guarded(handler)(btn));
      return btn;
    };

    // --- left side: navigation + secondary actions ---
    const prev = ui.button('← Back', () => {}) as HTMLButtonElement;
    prev.classList.add('cp-wizard-footer-btn');
    prev.addEventListener('click', guarded(() => this.prevStep())(prev));

    const reset = mkSecondary('🔄 Reset',
      'Reset all step state and clear the Mol* viewer',
      () => this.resetAll());
    const runPipeline = mkSecondary('▶▶ Run full pipeline',
      'Run all stages end-to-end (opens cluster picker after Align if needed)',
      () => this.runFullPipeline());

    const leftGroup = ui.divH([prev, reset, runPipeline]);
    leftGroup.style.gap = '8px';
    leftGroup.style.alignItems = 'center';

    // --- middle: status text + step indicator ---
    const status = ui.divText('', 'cp-wizard-footer-status');
    const stepIndicator = ui.divText('Step 1 of 5', 'cp-wizard-footer-step-indicator');
    const middleGroup = ui.divH([status, stepIndicator]);
    middleGroup.style.gap = '12px';
    middleGroup.style.alignItems = 'center';

    // --- right side: primary per-step actions ---
    const runStep = ui.div([], 'cp-wizard-footer-btn cp-wizard-footer-btn-primary');
    runStep.textContent = '▶ Run step';
    runStep.title = 'Run the active step';
    runStep.addEventListener('click', guarded(() => this.runCurrentStep())(runStep));

    const next = ui.button('Next →', () => {}) as HTMLButtonElement;
    next.classList.add('cp-wizard-footer-btn');
    next.addEventListener('click', guarded(() => this.nextStep())(next));

    const rightGroup = ui.divH([runStep, next]);
    rightGroup.style.gap = '8px';
    rightGroup.style.alignItems = 'center';

    const root = ui.divH([leftGroup, middleGroup, rightGroup], 'cp-wizard-footer');
    return {root, prev, next, status, runStep, runPipeline, reset, stepIndicator};
  }

  private buildMolstarPlaceholder(): HTMLElement {
    const icon = ui.divText('⚛', 'cp-wizard-molstar-placeholder-icon');
    this.molstarPlaceholderText = ui.divText('Run Step 1 to populate this viewer.');
    const t2 = ui.divText('This is ONE viewer — it persists across all 5 steps.');
    t2.style.fontSize = '11px';
    t2.style.opacity = '0.6';
    t2.style.marginTop = '6px';
    return ui.divV([icon, this.molstarPlaceholderText, t2], 'cp-wizard-molstar-placeholder');
  }

  // -------------------------------------------------------------------------
  // Step navigation
  // -------------------------------------------------------------------------

  private goToStep(step: StepId): void {
    // Block navigation while a step is running — the run owns the Mol* viewer
    // (it rebuilds it + adds overlays asynchronously), and navigating away
    // would start a second, concurrent render that races with the in-flight
    // one and corrupts the scene. This catches RAIL clicks (the footer Back/
    // Next are already disabled via setRunButtonsEnabled during a run).
    if (this.isStepRunning) {
      grok.shell.info('Please wait for the current step to finish before navigating.');
      return;
    }
    // Leaving Step 2 (Align) — release any single-PDB isolation so the
    // viewer shows all structures by default when the user returns later.
    if (this.activeStep === StepId.Align && step !== StepId.Align && this.isolatedPdb)
      this.isolatedPdb = null;

    // The "Active" marker must track exactly the step the user is on. If we're
    // leaving a step that was only `active` (visited via Next, but never run —
    // no result yet), revert it to `pending`. Otherwise it stays marked
    // "Active" in the rail after the user navigates back to an earlier step
    // (e.g. Step 3 → Step 2), misleadingly suggesting Step 3 is the current
    // step. `done`/`stale`/`error` are real results and are preserved.
    if (step !== this.activeStep && this.state.get(this.activeStep) === 'active')
      this.state.set(this.activeStep, 'pending');

    this.activeStep = step;
    // Only flip pending -> active; preserve done/stale/error.
    if (this.state.get(step) === 'pending') this.state.set(step, 'active');

    this.showStepContent(step);
    this.updateFooter(step);
    this.updateCaption(step);
    this.footerStepIndicator.textContent = `Step ${step + 1} of ${STEP_COUNT}`;
    this.contentRight.style.display = (step === StepId.Features) ? '' : 'none';
    // The Next button enable state depends on the new active step's status.
    this.applyNextButtonGuard();
    // Re-render the Mol* viewer to match the step the user just navigated to,
    // so the viewer never lingers on a previously-visited step's render.
    this.syncViewerToStep(step);
  }

  /**
   * Queue a viewer re-render so it matches `step`. Serialized through
   * `viewerSyncChain` — see that field's comment for why. The actual work is
   * in `doSyncViewerToStep`; this just chains it so concurrent navigations
   * never tear down a render in flight.
   */
  private syncViewerToStep(step: StepId): void {
    this.viewerSyncChain = this.viewerSyncChain
      .catch(() => { /* a prior failure must not break the chain */ })
      .then(() => this.doSyncViewerToStep(step));
  }

  /**
   * Make the Mol* viewer reflect `step`'s visualization. Runs one-at-a-time
   * (see `syncViewerToStep`).
   *
   * - Bails if the user has navigated past `step` since this was queued
   *   (`step !== activeStep`) or the viewer already shows it.
   * - Step 1 (Fetch): re-renders the fetched structures in their native
   *   frames if Step 1 is done; placeholder otherwise.
   * - Steps 2-5: re-render their visualization, but ONLY when `done` (a
   *   completed result is cached). For stale/pending/running the viewer is
   *   left as-is — re-running would be wrong (stale) or premature (pending).
   *
   * The per-step preview methods hit the per-fingerprint cache, so this is a
   * Mol*-only re-render with no Python round-trip when the result is cached.
   */
  private async doSyncViewerToStep(step: StepId): Promise<void> {
    // Stale-navigation guard: the user may have clicked Back/Next again while
    // this was queued behind a slower render. Only render the CURRENT step.
    if (step !== this.activeStep) return;
    if (step === this.lastRenderedStep) return;
    if (this.currentPdbIds.length === 0) return;

    if (step === StepId.Fetch) {
      try {
        if (this.state.get(StepId.Fetch) === 'done') {
          // If alignment has also run, show the aligned MULTI-STRUCTURE view
          // (same as Step 2) so the Step 1 Accepted PDBs row-click can isolate
          // a PDB and the "Use" checkbox can hide it — both via the proven
          // toggleVisibility fast path against tracked `_alignedStructureRefs`.
          // The legacy native-frame concat view doesn't track per-PDB refs, so
          // isolation falls into the slow rebuild path and races with leftover
          // state from earlier steps' renders (dual display, etc.).
          // No alignment yet → fall back to the legacy native-frame render so
          // the user still sees what they fetched.
          if (this.state.get(StepId.Align) === 'done')
            await this.app.previewAlign(this.currentPdbIds, this.currentOptions);
          else
            await this.app.previewFetch(this.currentPdbIds, this.currentOptions);
          // Re-check: the user may have navigated away during the async
          // render. Don't claim the viewer if Step 1 is no longer active.
          if (step !== this.activeStep) return;
          this.isolatedPdb = null;
        } else {
          this.showMolstarPlaceholder();
        }
        this.lastRenderedStep = StepId.Fetch;
      } catch (e) {
        console.error('[wizard] doSyncViewerToStep(Fetch) failed:', e);
      }
      return;
    }

    if (this.state.get(step) !== 'done') return;

    try {
      switch (step) {
      case StepId.Align:
        await this.app.previewAlign(this.currentPdbIds, this.currentOptions); break;
      case StepId.Pocket:
        await this.app.previewPocket(this.currentPdbIds, this.currentOptions); break;
      case StepId.Features:
        await this.app.previewFeatures(this.currentPdbIds, this.currentOptions); break;
      case StepId.Consensus:
        await this.app.previewConsensus(this.currentPdbIds, this.currentOptions); break;
      }
      // Re-check after the async render: if the user navigated to a different
      // step while previewX was running, don't apply this step's panel/state
      // — a later queued doSync (or its goToStep) owns the viewer now. Without
      // this, a slow render for the step we just LEFT would stomp the panel.
      if (step !== this.activeStep) return;
      this.lastRenderedStep = step;
      // A fresh render shows every structure, so any prior isolation is gone.
      this.isolatedPdb = null;
      this.showStepContent(step);
    } catch (e) {
      console.error('[wizard] doSyncViewerToStep failed:', e);
    }
  }

  /** Close the BSV viewer and show the empty-state placeholder in the center
   *  region. Used for Step 1 (no 3D) and after Reset. The primary text line
   *  is set to match context: before any run vs. after Step 1's QC. */
  private showMolstarPlaceholder(): void {
    this.app.clearViewer();
    // clearViewer() empties the host (contentCenter); re-add the placeholder.
    this.contentCenter.append(this.molstarPlaceholder);
    if (this.molstarPlaceholderText) {
      this.molstarPlaceholderText.textContent =
        this.state.get(StepId.Fetch) === 'done'
          ? 'Structures fetched & QC-filtered. Go to Align (Step 2) to superimpose them in 3D.'
          : 'Run Step 1 to fetch + QC the PDBs.';
    }
  }

  private prevStep(): void {
    if (this.activeStep > 0) this.goToStep((this.activeStep - 1) as StepId);
  }

  private nextStep(): void {
    if (this.activeStep < STEP_COUNT - 1)
      this.goToStep((this.activeStep + 1) as StepId);
  }

  private showStepContent(step: StepId): void {
    ui.empty(this.contentLeft);
    const ctx = this.makeStepContext(step);
    const builders = [buildStep1, buildStep2, buildStep3, buildStep4, buildStep5];
    const panel = builders[step](ctx);
    this.contentLeft.append(panel);
  }

  private updateFooter(step: StepId): void {
    if (step === 0) this.footerPrev.classList.add('cp-wizard-footer-btn-hidden');
    else this.footerPrev.classList.remove('cp-wizard-footer-btn-hidden');

    if (step === STEP_COUNT - 1) {
      this.footerNext.textContent = '✓ Finish';
      this.footerNext.title = 'Wizard complete — Mol* viewer shows the consensus pharmacophore.';
    } else {
      this.footerNext.textContent = `Next: ${STEP_LABELS[step + 1]} →`;
      this.footerNext.title = `Move to step ${step + 2}.`;
    }

    // Footer "Inputs changed — re-run to refresh" status removed: too
    // noisy on top of the step rail's yellow stale dot and the in-panel
    // "Pending reference: …" hint. The footerStatus element is kept blank
    // so the layout doesn't shift if it's referenced elsewhere.
    this.footerStatus.textContent = '';
    this.footerStatus.classList.remove('cp-wizard-footer-status-stale');
  }

  private updateCaption(step: StepId): void {
    const CAPTIONS = [
      'Step 1: Fetch + QC. Structures show in their native crystal frames — Align (Step 2) superimposes them.',
      'Step 2: Mol* shows all PDBs superimposed in a shared frame.',
      'Step 3: Pocket residues within cutoff of each drug ligand.',
      'Step 4: Chain F = ProLIF interaction sites, colored by family.',
      'Step 5: Chain P = consensus pharmacophore points (B-factor = frequency).',
    ];
    this.captionEl.textContent = CAPTIONS[step];
  }

  // -------------------------------------------------------------------------
  // Per-step context construction
  // -------------------------------------------------------------------------

  private makeStepContext(step: StepId): StepContext {
    return {
      pdbIds: this.currentPdbIds,
      options: this.currentOptions,
      isStale: this.state.get(step) === 'stale',
      pdbQcData: this.app.previewCacheView()?.pdbQc ?? null,
      pdbQcDropped: this.app.previewDroppedPdbs(),
      // Step 2 is "Align" = Stage 2a (global Cα Kabsch). Stage 2b is pocket
      // refinement and conceptually belongs to Step 3 (Pocket). Prefer the
      // Stage 2a output.
      alignedData: this.app.previewCacheView()?.aligned ??
                   this.app.previewCacheView()?.alignedV2 ?? null,
      pocketAtoms: this.app.previewCacheView()?.pocketAtoms ?? null,
      summaryData: this.app.previewSummaryDataFrame(),
      consensusData: this.app.previewCacheView()?.consensusModel ?? null,
      // Live count of what the current PDB/interaction exclusions leave for the
      // consensus — drives the Step 4 + Step 5 "consensus will use N / M" notes.
      consensusUsed: this.app.consensusUsedCount(this.currentOptions),
      seedPdbIds: step === StepId.Fetch ? this.currentPdbIds : undefined,
      demoLoaded: this.demoLoaded,
      onPdbInputChanged: (ids, demoLoaded) => {
        this.currentPdbIds = ids;
        if (demoLoaded != null) this.demoLoaded = demoLoaded;
        // Repaint the Run buttons — disable if textarea is empty, re-enable
        // as soon as a valid ID is entered. Editing the PDB list is a pending
        // change (no stale cascade): re-running Step 1 is what applies it and
        // stales the downstream steps.
        this.applyPdbIdsGuard();
      },
      isolatedPdb: this.isolatedPdb,
      onIsolatePdb: async (pdbId) => {
        // Ignore isolation while a step is running — isolatePdb mutates the
        // Mol* viewer (hide/show structures), which races with the run's
        // in-flight viewer rebuild and corrupts the scene.
        if (this.isStepRunning) return;
        this.isolatedPdb = pdbId;
        try { await this.app.isolatePdb(pdbId); } catch (e) {
          console.error('[wizard] isolate failed:', e);
        }
        // Re-render the step content so the table reflects the new selection.
        this.showStepContent(this.activeStep);
      },
      onTogglePdbVisibility: (pdbId, hidden) => {
        // The Step 4 "Use" checkbox both excludes a PDB from the consensus
        // (via onOptionsChanged) AND hides it in Mol*. Guard against a
        // concurrent run rebuilding the viewer (same reason as onIsolatePdb).
        if (this.isStepRunning) return;
        void this.app.setPdbHidden(pdbId, hidden);
      },
      onOptionsChanged: (next) => {
        // Selecting/editing an option (reference pick, pocket radius, …) is a
        // PENDING change — it does NOT mark anything stale. The displayed
        // results are still the last completed run's, and the per-step hint
        // ("Pending reference: X — click Run step") tells the user the choice
        // isn't applied yet. Staleness cascades only when a step is actually
        // RE-RUN (see runCurrentStep), so "click a row to check it" never
        // throws premature yellow stale flags on the rail.
        this.currentOptions = next;
      },
      // Re-render only the active step's content panel. Used by Step 2's
      // row-click-sets-reference flow to refresh the dropdown selection,
      // table row highlight, and hint text in one go. Deliberately scoped:
      // text/float inputs in other steps must NOT trigger a re-render on
      // every keystroke (it would steal focus from the input). The step
      // builders decide when to call this.
      refreshStep: () => this.showStepContent(this.activeStep),
    };
  }

  // -------------------------------------------------------------------------
  // Rail rendering
  // -------------------------------------------------------------------------

  private applyStatus(step: StepId, status: StepStatus): void {
    const stepEl = this.railSteps[step];
    const circle = this.railCircles[step];
    const sublabel = this.railSublabels[step];
    if (!stepEl || !circle || !sublabel) return;

    stepEl.classList.remove(
      'cp-step-pending', 'cp-step-active', 'cp-step-running',
      'cp-step-done', 'cp-step-stale', 'cp-step-error',
    );
    stepEl.classList.add(`cp-step-${status}`);

    // Remove any existing stale dot.
    stepEl.querySelectorAll('.cp-step-stale-dot').forEach((el) => el.remove());

    const INNER: Record<StepStatus, string> = {
      pending: CLOCK_SVG,
      active:  String(step + 1),
      running: String(step + 1),
      done:    CHECK_SVG,
      stale:   CHECK_SVG,
      error:   XMARK_SVG,
    };
    const SUB: Record<StepStatus, string> = {
      pending: 'Pending',
      active:  'Active',
      running: 'Running…',
      done:    'Done ✓',
      stale:   'Stale ⚠',
      error:   'Error',
    };
    circle.innerHTML = INNER[status];
    sublabel.textContent = SUB[status];

    if (status === 'stale') {
      const dot = document.createElement('div');
      dot.className = 'cp-step-stale-dot';
      const wrap = circle.parentElement;
      if (wrap) wrap.append(dot);
    }

    // Connector — colored for the connector before this step (left edge).
    if (step > 0) {
      const connector = this.railConnectors[step - 1];
      connector.classList.remove(
        'cp-wizard-step-connector-done', 'cp-wizard-step-connector-stale');
      const prev = this.state.get((step - 1) as StepId);
      if (prev === 'done') connector.classList.add('cp-wizard-step-connector-done');
      else if (prev === 'stale') connector.classList.add('cp-wizard-step-connector-stale');
    }
  }

  // -------------------------------------------------------------------------
  // Ribbon action handlers
  // -------------------------------------------------------------------------

  private async runCurrentStep(): Promise<void> {
    if (this.currentPdbIds.length === 0) {
      grok.shell.warning('Paste at least one PDB ID, or click "Load EGFR demo".');
      return;
    }
    const step = this.activeStep;
    this.setRunButtonsEnabled(false);
    this.isStepRunning = true;
    this.state.markRunning(step);
    this.hideMolstarStaleOverlay();
    try {
      switch (step) {
      case StepId.Fetch:
        // Step 1 shows the fetched structures in their native crystal frames
        // (the caption explains they're superimposed at Align). Earlier this
        // skipped the 3D entirely, but that left the viewer empty/confusing —
        // users expect to SEE their loaded proteins here.
        await this.app.previewFetch(this.currentPdbIds, this.currentOptions);
        break;
      case StepId.Align:
        await this.app.previewAlign(this.currentPdbIds, this.currentOptions);
        break;
      case StepId.Pocket:
        await this.app.previewPocket(this.currentPdbIds, this.currentOptions);
        break;
      case StepId.Features:
        await this.app.previewFeatures(this.currentPdbIds, this.currentOptions);
        break;
      case StepId.Consensus:
        await this.app.previewConsensus(this.currentPdbIds, this.currentOptions);
        break;
      }
      this.state.markDone(step);
      // Re-running a step commits its (possibly changed) inputs and makes its
      // result authoritative — so every DOWNSTREAM step that was already
      // `done` is now built on outdated upstream data and must be re-run. Flip
      // those to `stale`. `markStaleStartingFrom` only touches `done` steps, so
      // the normal forward flow (downstream still `pending`) is unaffected —
      // this fires only when the user went BACK and re-ran an earlier step.
      // This is now the ONLY place staleness is introduced: selecting/editing
      // options never stales, an actual re-run does.
      this.state.markStaleStartingFrom((step + 1) as StepId);
      // Snapshot inputs that produced these results.
      this.lastCommittedPdbIds = [...this.currentPdbIds];
      this.hideStaleBanner();
      // A successful re-run rebuilds Mol* showing every structure, so any
      // pre-run row-click isolation is no longer in effect visually. Clear
      // the wizard's `isolatedPdb` state to match.
      this.isolatedPdb = null;
      // The viewer now reflects this step's result.
      this.lastRenderedStep = step;
      // Re-render the step panel so the new data DataFrame renders.
      this.showStepContent(step);
    } catch (e) {
      console.error(`[wizard] step ${step} failed:`, e);
      this.state.markError(step);
      grok.shell.error(`Step ${step + 1} failed: ${(e as Error)?.message ?? e}`);
    } finally {
      this.isStepRunning = false;
      this.setRunButtonsEnabled(true);
    }
  }

  /** Public — exposed for the demo entry point in package.ts. */
  async runFullPipeline(): Promise<void> {
    if (this.currentPdbIds.length === 0) {
      grok.shell.warning('Paste at least one PDB ID, or click "Load EGFR demo".');
      return;
    }
    this.setRunButtonsEnabled(false);
    this.isStepRunning = true;
    this.hideMolstarStaleOverlay();
    // Mark all steps as running progressively. Orchestrator owns the actual
    // execution + the mid-pipeline cluster picker dialog.
    for (let s = 0; s < STEP_COUNT; s++) this.state.set(s as StepId, 'pending');
    this.state.markRunning(StepId.Fetch);
    try {
      await this.app.runPipeline(this.currentPdbIds, this.currentOptions);
      // Mark all steps done. runPipeline is monolithic — it doesn't emit
      // per-step events. The end state is "all done" if no exception.
      for (let s = 0; s < STEP_COUNT; s++) this.state.markDone(s as StepId);
      this.lastCommittedPdbIds = [...this.currentPdbIds];
      this.hideStaleBanner();
      // runPipeline already rendered the consensus (Stage 5b) into the viewer,
      // so mark it as the rendered step BEFORE navigating — otherwise
      // goToStep → syncViewerToStep would redundantly rebuild it.
      this.lastRenderedStep = StepId.Consensus;
      // Clear the run lock BEFORE navigating — goToStep is guarded against
      // isStepRunning, so the final hop to Consensus would otherwise no-op.
      this.isStepRunning = false;
      // Navigate to Step 5 (Consensus) so the user sees the final result.
      this.goToStep(StepId.Consensus);
    } catch (e) {
      console.error('[wizard] full pipeline failed:', e);
      this.state.markError(this.activeStep);
      grok.shell.error(`Pipeline failed: ${(e as Error)?.message ?? e}`);
    } finally {
      this.isStepRunning = false;
      this.setRunButtonsEnabled(true);
    }
  }

  private resetAll(): void {
    // 1. Orchestrator-side cleanup: close BSV viewer instance, drop cache,
    //    drop the cached summary df, unsubscribe + clear detail panel.
    this.app.resetWizardState();

    // 2. Step-state machine back to "Step 1 active, others pending".
    this.state.reset();

    // 3. Wizard-local state.
    this.activeStep = StepId.Fetch;
    this.currentPdbIds = [];
    this.currentOptions = {...DEFAULT_OPTIONS};
    this.demoLoaded = false;
    this.lastCommittedPdbIds = null;
    this.isolatedPdb = null;
    this.lastRenderedStep = null;

    // 4. Visual: hide stale UI, restore the Mol* placeholder, clear the
    //    detail panel and hide it (only visible on Step 4 anyway).
    this.hideStaleBanner();
    this.hideMolstarStaleOverlay();
    ui.empty(this.contentCenter);
    this.contentCenter.append(this.molstarPlaceholder);
    ui.empty(this.contentRight);
    this.contentRight.style.display = 'none';

    // 5. Re-render Step 1 with a fresh (empty) context.
    this.goToStep(StepId.Fetch);

    // 6. Re-apply the PDB-IDs guard (currentPdbIds is now empty, so this
    //    disables Run step + Run full pipeline).
    this.applyPdbIdsGuard();

    grok.shell.info('Wizard reset — Mol* viewer cleared, PDB list emptied. ' +
      'Load the demo or paste PDB IDs to start over.');
  }

  // -------------------------------------------------------------------------
  // Stale banner + Mol* overlay
  // -------------------------------------------------------------------------

  /**
   * Single source of truth for "is the current form state stale?" — called
   * whenever the PDB list or any pipeline option changes.
   *
   * Logic:
   *  - If no run has happened yet → nothing to be stale about; no-op.
   *  - If the current `(pdbIds, options)` fingerprint MATCHES the fingerprint
   *    that produced the cached results → un-stale every step that still has
   *    cached data (revert from 'stale' to 'done'), hide the banner + Mol*
   *    overlay. This is the "user edited an option and edited it back" case.
   *  - Otherwise (fingerprint diverges from the last run) → mark the pipeline
   *    stale from Fetch downward; show the banner + Mol* overlay; keep the
   *    last run's DataFrames in the cache so the step panels stay populated.
   *
   * NOTE: we deliberately do NOT call `invalidatePreviewCache()`. The
   * fingerprint logic in `orchestrator.ensurePreviewCache` already handles
   * cache reuse vs recompute; keeping the cache around lets the user revert
   * back to a previous state without losing the displayed tables.
   */
  private reconcileStaleState(startStep: StepId = StepId.Fetch): void {
    const lastFp = this.app.lastRunFingerprint();
    if (lastFp == null) {
      // No prior run — nothing to reconcile. Make sure the banner isn't
      // somehow still up (e.g. after a Reset that left it open).
      this.hideStaleBanner();
      this.hideMolstarStaleOverlay();
      return;
    }
    const currentFp = this.app.previewFingerprint(this.currentPdbIds, this.currentOptions);
    if (currentFp === lastFp) {
      // Form state matches the last completed run — clear stale on every
      // step that still has cached data.
      const cache = this.app.previewCacheView();
      const hasCached: boolean[] = [
        !!cache?.pdbQc,
        !!cache?.aligned,
        !!cache?.pocketAtoms,
        !!cache?.ligandFeatures,
        !!cache?.consensusModel,
      ];
      this.state.clearStaleIfCached(hasCached);
      this.hideStaleBanner();
      this.hideMolstarStaleOverlay();
      // Repaint the active step's content so its hint text / row highlights
      // reflect the no-longer-pending state (e.g. Step 2 reverts from
      // "Pending reference: X" to "Reference: X").
      this.showStepContent(this.activeStep);
      return;
    }
    // Fingerprint diverges. Mark stale starting from `startStep` — the
    // EARLIEST step whose result is invalidated by the change. This is
    // pivotal for UX: changing a Pocket option (Step 3 input) must NOT
    // mark Step 2 (Align) stale, because the alignment is still the
    // result that's currently displayed in the viewer and is unchanged
    // by a Pocket-options edit.
    this.state.markStaleStartingFrom(startStep);
  }

  /**
   * Diff two PipelineOptions and return the earliest step whose result is
   * invalidated. Drives the per-option stale cascade:
   *
   *   Step 1 (Fetch)     ← QC filter changes (resolution, methods, ligand MW)
   *   Step 2 (Align)     ← refPdbId
   *   Step 3 (Pocket)    ← pocketMethod, pocketRadius, pocketRep
   *   Step 5 (Consensus) ← kq, minClusterSizeFraction, topClusterNumber
   *
   * Returns StepId.Fetch as a defensive default if the change doesn't fall
   * into any known bucket — better to over-invalidate than to silently
   * present stale data as fresh.
   */
  private earliestAffectedStep(prev: PipelineOptions, next: PipelineOptions): StepId | null {
    if (prev.maxResolution !== next.maxResolution ||
        prev.allowXray !== next.allowXray ||
        prev.allowNmr !== next.allowNmr ||
        prev.allowCryoEm !== next.allowCryoEm ||
        prev.allowAlphaFold !== next.allowAlphaFold ||
        prev.minLigandMw !== next.minLigandMw ||
        // The legacy requireXray flag is still in PipelineOptions for type
        // compat; include it for safety even though the UI doesn't expose it.
        prev.requireXray !== next.requireXray) return StepId.Fetch;
    if ((prev.refPdbId ?? '').toUpperCase() !== (next.refPdbId ?? '').toUpperCase())
      return StepId.Align;
    if (prev.pocketMethod !== next.pocketMethod ||
        prev.pocketRadius !== next.pocketRadius ||
        prev.pocketRep !== next.pocketRep) return StepId.Pocket;
    if (prev.kq !== next.kq ||
        prev.minClusterSizeFraction !== next.minClusterSizeFraction ||
        prev.topClusterNumber !== next.topClusterNumber) return StepId.Consensus;
    // No tracked field changed. Return null so the caller can skip the
    // stale cascade entirely — typical cause is the belt-and-braces raw
    // <select> 'change' listener firing the same handler twice (the second
    // call sees prev === next because currentOptions was already updated
    // by the first). Returning a default StepId.Fetch would blanket-stale
    // every step in this no-op case, which was the symptom in the QA pass.
    return null;
  }

  // Stale banner removed: noisy on top of the step rail's yellow stale dot
  // and the in-panel "Pending reference: …" hint. The functions are kept as
  // no-ops so the existing call sites compile; staleBannerSlot is left in
  // the DOM as an empty placeholder (cheap, keeps layout stable).
  // eslint-disable-next-line @typescript-eslint/no-empty-function
  private showStaleBanner(): void {}

  private hideStaleBanner(): void {
    // Belt-and-braces — clean up if an older build left a banner element
    // in the slot. Safe no-op when nothing's there.
    this.staleBannerVisible = false;
    ui.empty(this.staleBannerSlot);
    this.staleBannerEl = undefined;
  }

  // Mol* stale overlay removed: it was visually distracting and redundant
  // with the stale banner above and the footer status. The functions are
  // kept as no-ops so the existing call sites compile; they can be deleted
  // entirely in a later cleanup pass. The defensive `hideMolstarStaleOverlay`
  // call also cleans up any stale overlay element that survived a hot reload.
  // eslint-disable-next-line @typescript-eslint/no-empty-function
  private showMolstarStaleOverlay(): void {}

  private hideMolstarStaleOverlay(): void {
    // Belt-and-braces cleanup in case an older build left an overlay element
    // in the DOM. Safe no-op when nothing matches.
    this.contentCenter.querySelector('.cp-wizard-molstar-stale-overlay')?.remove();
  }

  // -------------------------------------------------------------------------
  // Helpers
  // -------------------------------------------------------------------------

  private setRunButtonsEnabled(enabled: boolean): void {
    const toggle = (el: HTMLElement, on: boolean) => {
      if (on) el.classList.remove('cp-wizard-footer-btn-disabled');
      else el.classList.add('cp-wizard-footer-btn-disabled');
    };
    toggle(this.footerRunStep, enabled);
    toggle(this.footerRunPipeline, enabled);
    toggle(this.footerReset, enabled);
    toggle(this.footerNext, enabled);
    toggle(this.footerPrev, enabled);
    // When re-enabling, fold in the PDB-IDs guard so we don't accidentally
    // enable Run buttons while the textarea is empty.
    if (enabled) this.applyPdbIdsGuard();
  }

  /** Disable Run step + Run full pipeline when no PDB IDs are entered.
   *  Adds a tooltip explaining the requirement. Called whenever the PDB
   *  list changes (textarea edit, demo load, target lookup, reset) and on
   *  step navigation. */
  private applyPdbIdsGuard(): void {
    const RUN_STEP_DEFAULT_TITLE = 'Run the active step';
    const RUN_PIPELINE_DEFAULT_TITLE =
      'Run all stages end-to-end (opens cluster picker after Align if needed)';
    const EMPTY_PDB_TITLE = 'Paste at least one PDB ID, or click "Load EGFR demo (6 PDBs)".';

    const empty = this.currentPdbIds.length === 0;
    const toggleDisabled = (el: HTMLElement, disable: boolean) => {
      if (disable) el.classList.add('cp-wizard-footer-btn-disabled');
      else el.classList.remove('cp-wizard-footer-btn-disabled');
    };

    toggleDisabled(this.footerRunStep, empty);
    toggleDisabled(this.footerRunPipeline, empty);
    this.footerRunStep.title = empty ? EMPTY_PDB_TITLE : RUN_STEP_DEFAULT_TITLE;
    this.footerRunPipeline.title = empty ? EMPTY_PDB_TITLE : RUN_PIPELINE_DEFAULT_TITLE;

    // Next/Finish also depends on the current step being done at least once.
    this.applyNextButtonGuard();
  }

  /** Disable Next (or Finish on the last step) until the current step has
   *  been completed at least once. `done` and `stale` both count as
   *  completed — stale means the user can still navigate forward to inspect
   *  old results. `pending`, `active`, `running`, and `error` all block. */
  private applyNextButtonGuard(): void {
    const status = this.state.get(this.activeStep);
    const isCompleted = status === 'done' || status === 'stale';
    const isLast = this.activeStep === STEP_COUNT - 1;
    const labelDefault = isLast ?
      'Finish — return to the consensus pharmacophore.' :
      `Move to step ${this.activeStep + 2}.`;
    const labelBlocked = isLast ?
      'Run step 5 first to compute the consensus pharmacophore.' :
      `Run step ${this.activeStep + 1} first to unlock the next step.`;

    if (isCompleted) {
      this.footerNext.classList.remove('cp-wizard-footer-btn-disabled');
      this.footerNext.title = labelDefault;
    } else {
      this.footerNext.classList.add('cp-wizard-footer-btn-disabled');
      this.footerNext.title = labelBlocked;
    }
  }
}
