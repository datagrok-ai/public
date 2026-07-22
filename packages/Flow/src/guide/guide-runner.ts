/** Runs a `Guide`: highlights each step's target and shows an instruction popup
 *  that auto-positions to stay fully on screen, then waits for the user to
 *  perform the step's action before advancing. Exit/Skip are always available.
 *  One guide at a time.
 *
 *  The popup is our own element (not `ui.hints.addHint`) so we control a single
 *  close affordance — the platform popup injects its own ✕ which collided with
 *  our Exit link — and the placement, which the platform's helper gets wrong for
 *  targets near a viewport edge (it never flips to the side with room). */

import * as ui from 'datagrok-api/ui';
import {
  Guide, GuideStep, GuideContext, GuideHost, isAborted, Side, computePlacement, openDialogEl,
  isScrolledIntoView, prefillSearch,
} from './guide-model';
import {setTid} from '../utils/test-ids';

type StepOutcome = 'next' | 'exit';

export class GuideRunner {
  private controller: AbortController | null = null;

  get isRunning(): boolean {
    return this.controller !== null && !this.controller.signal.aborted;
  }

  /** Abort any running guide (e.g. the user launched another, or left the view). */
  stop(): void {
    const wasRunning = this.isRunning;
    this.controller?.abort();
    this.controller = null;
    GuideRunner.clearAllHighlights();
    if (wasRunning) prefillSearch('');
  }

  async run(guide: Guide, host: GuideHost): Promise<void> {
    this.stop();
    const controller = new AbortController();
    this.controller = controller;
    const ctx: GuideContext = {host, signal: controller.signal};
    // The start overlay covers the canvas every early step points at.
    try {
      host.hideStartPanel?.();
    } catch { /* optional */ }

    // Steps skipped via skipIf must not leave holes in the "Step i of n"
    // counter ("Step 9" followed by "Step 11" reads as a missed instruction).
    // Show a running count instead, with the total estimated over the steps
    // whose skipIf isn't already satisfied — it self-corrects as state changes.
    let shown = 0;
    const remainingEstimate = (from: number): number => {
      let count = 0;
      for (let k = from; k < guide.steps.length; k++) {
        try {
          if (!guide.steps[k].skipIf?.(ctx)) count++;
        } catch {
          count++;
        }
      }
      return count;
    };

    for (let i = 0; i < guide.steps.length; i++) {
      if (controller.signal.aborted) return;
      const step = guide.steps[i];
      try {
        if (step.skipIf?.(ctx)) continue;
      } catch {/* a throwing predicate counts as not-satisfied */}
      shown++;
      const total = shown + remainingEstimate(i + 1);
      const outcome = await this.runStep(step, shown - 1, total, ctx, () => controller.abort());
      if (outcome === 'exit' || controller.signal.aborted) {
        this.controller = null;
        return;
      }
    }
    if (!controller.signal.aborted) this.showCompletion(guide, host);
    this.controller = null;
  }

  private async runStep(
    step: GuideStep, i: number, n: number, ctx: GuideContext, abort: () => void,
  ): Promise<StepOutcome> {
    // Nuke any stragglers from a previous step before we add ours — covers the
    // case where a step resolved instantly and a queued frame re-applied a
    // highlight after cleanup (see clearAllHighlights). Also drop any platform
    // tooltip left hanging by the previous step's hover — a stale "Run the
    // flow" tooltip beside a fresh highlight misleads.
    GuideRunner.clearAllHighlights();
    try {
      ui.tooltip.hide();
    } catch { /* tooltip host not ready */ }
    // Prerequisite steps declare `skipIf`: when already satisfied, skip silently
    // (no card, no setup) — that's what makes "ensure X exists" steps invisible
    // when X is already there.
    try {
      if (step.skipIf?.(ctx)) return 'next';
    } catch {/* a throwing predicate counts as not-satisfied */}
    try {
      await step.setup?.(ctx);
    } catch { /* setup is best-effort */ }
    if (ctx.signal.aborted) return 'exit';

    // What the popup anchors to; what gets highlighted (defaults to the anchor,
    // but a step may highlight several things — e.g. both pins to connect).
    const anchorOf = (): HTMLElement | null => step.target?.(ctx) ?? null;
    const highlightsOf = (): HTMLElement[] => {
      const list = step.highlights ? step.highlights(ctx) : [step.target?.(ctx) ?? null];
      const out: HTMLElement[] = [];
      for (const e of list) if (e) out.push(e);
      return out;
    };

    // One pulsing dot per highlighted element, pinned to its top-left corner.
    // Body-level (not a box-shadow on the element, which the canvas clips and
    // node themes override) and kept stable across ticks so the pulse doesn't
    // restart. The map is the source of truth for cleanup.
    const blobByEl = new Map<HTMLElement, HTMLElement>();
    const syncHighlights = (): void => {
      const want = highlightsOf();
      for (const [el, blob] of blobByEl) {
        if (!want.includes(el) || !document.body.contains(el)) {
          el.classList.remove('ff-guide-target', 'ff-guide-target-large');
          blob.remove();
          blobByEl.delete(el);
        }
      }
      for (const el of want) {
        if (!blobByEl.has(el)) {
          el.classList.add('ff-guide-target');
          const blob = document.createElement('div');
          blob.className = 'ff-guide-blob';
          document.body.appendChild(blob);
          blobByEl.set(el, blob);
        }
        // Big containers (a whole toolbox pane, the canvas) get an outline-only
        // highlight — the orange fill washes over every row inside and reads
        // as "everything is broken" rather than "look here".
        const r = el.getBoundingClientRect();
        el.classList.toggle('ff-guide-target-large', r.width * r.height > 30000);
      }
      for (const [el, blob] of blobByEl) {
        const r = el.getBoundingClientRect();
        if (r.width > 0 && r.height > 0) {
          blob.style.display = 'block';
          blob.style.left = `${Math.round(r.left)}px`;
          blob.style.top = `${Math.round(r.top)}px`;
        } else {
          blob.style.display = 'none';
        }
      }
    };
    const clearHighlights = (): void => {
      for (const [el, blob] of blobByEl) {
        el.classList.remove('ff-guide-target', 'ff-guide-target-large');
        blob.remove();
      }
      blobByEl.clear();
    };

    let resolveControl!: (o: StepOutcome) => void;
    const control = new Promise<StepOutcome>((res) => {
      resolveControl = res;
    });
    const onExit = (): void => {
      abort();
      resolveControl('exit');
    };
    const card = this.buildCard(step, i, n, !step.until, () => resolveControl('next'), onExit);
    document.body.appendChild(card);

    // If the anchor exists but sits scrolled out of view inside a scrollable
    // pane (a toolbox item below the fold, a context-panel row), bring it in —
    // otherwise the highlight lands on a clipped, invisible element. Canvas
    // nodes aren't inside a scroll container, so this is a no-op for them.
    // The one step that TEACHES scrolling targets the pane itself until its
    // file is visible, so it is not short-circuited by this.
    const anchor0 = anchorOf();
    if (anchor0 && !isScrolledIntoView(anchor0)) {
      try {
        anchor0.scrollIntoView({block: 'nearest'});
      } catch { /* detached mid-step — the reanchor timer recovers */ }
    }

    // Re-anchor on a timer: nodes re-render (replacing their DOM element), the
    // context panel opening shifts layout, and the user may drag the target.
    // Re-resolving each tick keeps highlights + popup glued to the live
    // elements and reflows the popup onto whichever side currently has room.
    let finished = false;
    const reanchor = (): void => {
      if (finished) return; // a queued frame must never re-highlight after cleanup
      syncHighlights();
      // A target-less step while a dialog is open must anchor BESIDE the
      // dialog — centered, it would sit underneath it (the card drops below
      // dialog z-index so the dialog stays interactive).
      const dialog = openDialogEl();
      const extraAvoid = (step.avoid?.(ctx) ?? []).filter((e): e is HTMLElement => !!e);
      this.place(card, anchorOf() ?? dialog, step.position,
        [...blobByEl.keys(), ...extraAvoid, ...(dialog ? [dialog] : [])]);
      card.style.zIndex = dialog ? '2900' : '5000';
    };
    reanchor();
    const timer = window.setInterval(reanchor, 250);
    const onResize = (): void => reanchor();
    window.addEventListener('resize', onResize);
    // One reflow after layout settles (offsetWidth is 0 until painted).
    const raf = requestAnimationFrame(reanchor);

    let outcome: StepOutcome = 'next';
    try {
      if (step.until) {
        const done = step.until(ctx)
          .then<StepOutcome>(() => 'next')
          .catch<StepOutcome>((e) => (isAborted(e) ? 'exit' : Promise.reject(e)));
        outcome = await Promise.race([done, control]);
      } else {
        outcome = await control;
      }
    } catch {
      outcome = 'exit';
    } finally {
      finished = true;
      window.clearInterval(timer);
      cancelAnimationFrame(raf);
      window.removeEventListener('resize', onResize);
      clearHighlights();
      GuideRunner.clearAllHighlights();
      card.remove();
    }
    return outcome;
  }

  /** Remove every guide highlight from the DOM — the target tint/outline and all
   *  pulsing dots — regardless of which step created them. Called on every step
   *  boundary and on finish/exit so a highlight can never linger (e.g. when a
   *  step resolves instantly and its per-step cleanup races a queued frame). */
  static clearAllHighlights(): void {
    document.querySelectorAll('.ff-guide-target')
      .forEach((el) => el.classList.remove('ff-guide-target', 'ff-guide-target-large'));
    document.querySelectorAll('.ff-guide-blob').forEach((el) => el.remove());
  }

  /** Place `popup` next to `target`, choosing the side with room and clamping it
   *  fully into the viewport. `preferred` is honored when it fits, else we flip
   *  to its opposite, else fall back to right/left/bottom/top — so a target low
   *  on screen gets a popup above it, one far left gets it to the right, etc.
   *  Sides whose popup would sit on top of a highlighted element (e.g. the
   *  other socket of a connect step) lose to sides that keep every highlight
   *  visible. */
  private place(popup: HTMLElement, target: HTMLElement | null, preferred?: Side,
    avoidEls: HTMLElement[] = []): void {
    const pw = popup.offsetWidth || 300;
    const ph = popup.offsetHeight || 170;
    const rect = target && document.body.contains(target) ? target.getBoundingClientRect() : null;
    // A detached/hidden target has a zero-size rect → fall back to centered.
    const usable = rect && rect.width > 0 && rect.height > 0 ? rect : null;
    const avoid = avoidEls
      .filter((el) => document.body.contains(el))
      .map((el) => el.getBoundingClientRect())
      .filter((r) => r.width > 0 && r.height > 0);
    const p = computePlacement(usable, pw, ph, window.innerWidth, window.innerHeight, preferred,
      undefined, undefined, avoid);
    popup.style.left = `${p.x}px`;
    popup.style.top = `${p.y}px`;
    popup.dataset.side = p.side;
  }

  /** The instruction popup body. A single ✕ (top-right) exits the guide; the
   *  footer offers Next/Finish (manual steps) or Skip (action steps). */
  private buildCard(
    step: GuideStep, i: number, n: number, manual: boolean,
    onNext: () => void, onExit: () => void,
  ): HTMLElement {
    const progress = ui.divText(`Step ${i + 1} of ${n}`, 'ff-guide-progress');
    const close = ui.iconFA('times', onExit, 'Exit the tour');
    close.classList.add('ff-guide-close');
    setTid(close, 'guide-exit');
    const header = ui.divH([progress, close], 'ff-guide-header');

    const title = ui.divText(step.title, 'ff-guide-title');
    const text = ui.divText(step.text, 'ff-guide-text');

    const footer = ui.divH([], 'ff-guide-footer');
    if (manual) {
      const next = ui.button(i + 1 < n ? 'Next' : 'Finish', onNext, 'Continue');
      setTid(next, 'guide-next');
      footer.appendChild(next);
    } else {
      const skip = ui.link('Skip this step', onNext, 'Move on without doing it', 'ff-guide-skip');
      setTid(skip, 'guide-skip');
      footer.appendChild(skip);
    }

    const card = setTid(ui.divV([header, title, text, footer], 'ff-guide-card'), 'guide-card');
    card.style.position = 'fixed';
    card.style.zIndex = '5000';
    return card;
  }

  private showCompletion(guide: Guide, host: GuideHost): void {
    GuideRunner.clearAllHighlights();
    // Steps prefill the toolbox search to spotlight one function — don't leave
    // the catalog filtered once the guide is over.
    prefillSearch('');
    let timer = 0;
    const close = (): void => {
      window.clearTimeout(timer);
      card.remove();
    };
    const done = ui.button('Done', close, 'Close');
    setTid(done, 'guide-done');
    const card = setTid(ui.divV([
      ui.divText('All done! 🎉', 'ff-guide-title'),
      ui.divText(`You finished “${guide.title}”. Open the help button anytime for more.`, 'ff-guide-text'),
      ui.divH([done], 'ff-guide-footer'),
    ], 'ff-guide-card ff-guide-card-done'), 'guide-card');
    card.style.position = 'fixed';
    card.style.zIndex = '5000';
    document.body.appendChild(card);
    this.place(card, host.anchorEl, 'top');
    // Auto-dismiss so the "all done" note doesn't linger at the bottom.
    timer = window.setTimeout(close, 5000);
  }
}
