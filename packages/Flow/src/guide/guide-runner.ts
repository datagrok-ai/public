/** Runs a `Guide`: highlight each step's target, show an instruction popup, wait for the action.
 *  The popup is NOT `ui.hints.addHint` — that injects its own ✕ and never flips to the side with room. */

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

    // Skipped steps must not leave holes in the "Step i of n" counter — show a running
    // count, with the total re-estimated over the not-yet-satisfied steps.
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
    // Stragglers from a previous step: a queued frame may re-apply a highlight after
    // cleanup, and a stale platform tooltip beside a fresh highlight misleads.
    GuideRunner.clearAllHighlights();
    try {
      ui.tooltip.hide();
    } catch { /* tooltip host not ready */ }
    try {
      if (step.skipIf?.(ctx)) return 'next';
    } catch {/* a throwing predicate counts as not-satisfied */}
    try {
      await step.setup?.(ctx);
    } catch { /* setup is best-effort */ }
    if (ctx.signal.aborted) return 'exit';

    const anchorOf = (): HTMLElement | null => step.target?.(ctx) ?? null;
    const highlightsOf = (): HTMLElement[] => {
      const list = step.highlights ? step.highlights(ctx) : [step.target?.(ctx) ?? null];
      const out: HTMLElement[] = [];
      for (const e of list) if (e) out.push(e);
      return out;
    };

    // Body-level pulse dots (a box-shadow on the element is clipped by the canvas),
    // kept stable across ticks so the pulse doesn't restart.
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
        // Big containers get outline-only — the orange fill washes over every row inside.
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
    // stop() must settle this step even when its `until` ignores the signal —
    // otherwise the finally below never runs and the card outlives the view.
    const onAbort = (): void => resolveControl('exit');
    ctx.signal.addEventListener('abort', onAbort, {once: true});
    const onExit = (): void => {
      abort();
      resolveControl('exit');
    };
    const card = this.buildCard(step, i, n, !step.until, () => resolveControl('next'), onExit);
    document.body.appendChild(card);

    // Bring a clipped anchor into view; the step that TEACHES scrolling targets
    // the pane itself, so it is not short-circuited by this.
    const anchor0 = anchorOf();
    if (anchor0 && !isScrolledIntoView(anchor0)) {
      try {
        anchor0.scrollIntoView({block: 'nearest'});
      } catch { /* detached mid-step — the reanchor timer recovers */ }
    }

    // Re-anchor on a timer: nodes re-render, layout shifts, the user drags the target.
    let finished = false;
    const reanchor = (): void => {
      if (finished) return; // a queued frame must never re-highlight after cleanup
      syncHighlights();
      // A target-less step while a dialog is open must anchor BESIDE the dialog —
      // centered, it would sit underneath (the card drops below dialog z-index).
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
      ctx.signal.removeEventListener('abort', onAbort);
      window.clearInterval(timer);
      cancelAnimationFrame(raf);
      window.removeEventListener('resize', onResize);
      clearHighlights();
      GuideRunner.clearAllHighlights();
      card.remove();
    }
    return outcome;
  }

  /** Called on every step boundary and on finish/exit so a highlight can never linger. */
  static clearAllHighlights(): void {
    document.querySelectorAll('.ff-guide-target')
      .forEach((el) => el.classList.remove('ff-guide-target', 'ff-guide-target-large'));
    document.querySelectorAll('.ff-guide-blob').forEach((el) => el.remove());
  }

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
    // Don't leave the catalog filtered by a step's search prefill.
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
    timer = window.setTimeout(close, 5000);
  }
}
