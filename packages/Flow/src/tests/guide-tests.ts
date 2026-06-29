/** Tests for the in-app guide system: the pure condition helpers (poll, click,
 *  section-expanded) and content integrity of the tutorials/questions. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  poll, waitForClick, isAborted, untilSectionExpanded, computePlacement,
  byNodeFunc, byParam, paramFieldSelector, untilValueContains, untilValueMatches,
  GuideContext, GuideHost,
} from '../guide/guide-model';
import {TUTORIALS, QUESTIONS} from '../guide/guide-content';
import {GuideRunner} from '../guide/guide-runner';
import {Guide} from '../guide/guide-model';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {makeEditor, destroyEditor, until} from './test-utils';

/** Build a PlaceRect from x/y/w/h. */
function rect(x: number, y: number, w: number, h: number) {
  return {left: x, top: y, right: x + w, bottom: y + h, width: w, height: h};
}

function fakeHost(): GuideHost {
  const el = document.createElement('div');
  return {getFlow: () => undefined, showFunctionBrowser: () => {}, anchorEl: el};
}

category('Flow: guide', () => {
  test('poll resolves when the predicate becomes true', async () => {
    const ac = new AbortController();
    let flag = false;
    setTimeout(() => {flag = true;}, 40);
    await poll(() => flag, ac.signal, 15); // resolves or the test times out
    expect(flag, true);
  });

  test('poll rejects (AbortedError) when the signal aborts', async () => {
    const ac = new AbortController();
    const p = poll(() => false, ac.signal, 15);
    ac.abort();
    let aborted = false;
    try {
      await p;
    } catch (e) {
      aborted = isAborted(e);
    }
    expect(aborted, true);
  });

  test('waitForClick resolves on a click of the (late-appearing) target', async () => {
    const ac = new AbortController();
    const ctx: GuideContext = {host: fakeHost(), signal: ac.signal};
    const div = document.createElement('div');
    document.body.appendChild(div);
    try {
      const p = waitForClick(() => div, ctx);
      setTimeout(() => div.click(), 30);
      await p;
      expect(true, true);
    } finally {
      div.remove();
    }
  });

  test('untilValueMatches ignores case and whitespace', async () => {
    const ac = new AbortController();
    const ctx: GuideContext = {host: fakeHost(), signal: ac.signal};
    const input = document.createElement('input');
    input.setAttribute('data-testid', 'tmp-search');
    document.body.appendChild(input);
    try {
      const p = untilValueMatches('[data-testid="tmp-search"]', 'openfile')(ctx);
      // The user types the spaced, mixed-case form — still satisfies it.
      setTimeout(() => {input.value = 'Open File';}, 30);
      await p;
      expect(true, true, 'resolved on "Open File" for needle "openfile"');
    } finally {
      ac.abort();
      input.remove();
    }
  });

  test('a step that resolves instantly leaves no lingering highlight', async () => {
    // Regression: when `until` is already true, the step's cleanup raced a queued
    // frame that re-applied the highlight, leaving the orange tint/dot forever.
    const runner = new GuideRunner();
    const target = document.createElement('div');
    target.style.cssText = 'position:fixed;left:10px;top:10px;width:40px;height:20px;';
    document.body.appendChild(target);
    const guide: Guide = {
      id: 't', kind: 'tutorial', title: 'T', summary: 's',
      steps: [{title: 'a', text: 'a', target: () => target, until: () => Promise.resolve()}],
    };
    try {
      await runner.run(guide, fakeHost());
      // Let any queued frame / interval fire — the fix must keep them inert.
      await new Promise((r) => requestAnimationFrame(() => r(null)));
      await new Promise((r) => setTimeout(r, 30));
      expect(document.querySelectorAll('.ff-guide-target').length, 0, 'no lingering target highlight');
      expect(document.querySelectorAll('.ff-guide-blob').length, 0, 'no lingering pulse dots');
    } finally {
      runner.stop();
      target.remove();
      document.querySelectorAll('[data-testid="ff-guide-card"]').forEach((e) => e.remove());
    }
  });

  test('a prerequisite step whose skipIf is satisfied is skipped (no hang)', async () => {
    // The prereq's `until` never resolves; only skipIf lets the guide proceed.
    const runner = new GuideRunner();
    const guide: Guide = {
      id: 'sk', kind: 'question', title: 'Q', summary: 's',
      steps: [
        {title: 'prereq', text: 'p', skipIf: () => true, until: () => new Promise<void>(() => { /* never */ })},
        {title: 'answer', text: 'a'},
      ],
    };
    try {
      const p = runner.run(guide, fakeHost());
      await new Promise((r) => setTimeout(r, 60));
      const next = document.querySelector('[data-testid="ff-guide-next"]') as HTMLElement | null;
      expect(!!next, true, 'reached the manual answer step — prereq was skipped, not stuck');
      next?.click();
      await p;
    } finally {
      runner.stop();
      document.querySelectorAll('[data-testid="ff-guide-card"]').forEach((e) => e.remove());
    }
  });

  test('untilSectionExpanded reacts to the header collapsed class', async () => {
    const ac = new AbortController();
    const ctx: GuideContext = {host: fakeHost(), signal: ac.signal};
    const header = document.createElement('div');
    header.className = 'funcflow-section-header collapsed';
    header.dataset.section = 'Inputs';
    document.body.appendChild(header);
    try {
      const p = untilSectionExpanded('Inputs')(ctx);
      setTimeout(() => header.classList.remove('collapsed'), 30);
      await p;
      expect(header.classList.contains('collapsed'), false);
    } finally {
      header.remove();
    }
  });

  test('placement: no target → centered horizontally', async () => {
    const p = computePlacement(null, 300, 160, 1200, 800);
    expect(p.side, 'center');
    expect(p.x, Math.round((1200 - 300) / 2));
  });

  test('placement: target hugging the left edge → popup goes right', async () => {
    // A toolbar icon at the far left: 'left' can't fit, 'right' can.
    const p = computePlacement(rect(4, 60, 30, 30), 300, 160, 1200, 800, 'left');
    expect(p.side, 'right');
    expect(p.x >= 10, true);
  });

  test('placement: target low on screen → popup flips above, stays on screen', async () => {
    // A node near the bottom; 'bottom' would overflow, so flip to 'top'.
    const vh = 800;
    const p = computePlacement(rect(500, 760, 120, 30), 300, 200, 1200, vh, 'bottom');
    expect(p.side, 'top');
    expect(p.y >= 10, true, 'top edge on screen');
    expect(p.y + 200 <= vh - 10, true, 'bottom edge on screen');
  });

  test('placement: result is always fully inside the viewport', async () => {
    const vw = 1000; const vh = 700; const pw = 300; const ph = 180;
    const cases = [rect(0, 0, 20, 20), rect(980, 680, 20, 20), rect(500, 0, 40, 10),
      rect(500, 690, 40, 10), rect(-50, 350, 30, 30), rect(990, 350, 30, 30)];
    for (const r of cases) {
      for (const pref of ['top', 'bottom', 'left', 'right'] as const) {
        const p = computePlacement(r, pw, ph, vw, vh, pref);
        expect(p.x >= 10 && p.x + pw <= vw - 10, true, `x in bounds for ${JSON.stringify(r)}/${pref}`);
        expect(p.y >= 10 && p.y + ph <= vh - 10, true, `y in bounds for ${JSON.stringify(r)}/${pref}`);
      }
    }
  });

  test('content: 4 tutorials, 10+ questions, unique ids, every step has text', async () => {
    expect(TUTORIALS.length, 4);
    expect(QUESTIONS.length >= 10, true, `expected 10+ questions, got ${QUESTIONS.length}`);

    const all = [...TUTORIALS, ...QUESTIONS];
    const ids = all.map((g) => g.id);
    expect(new Set(ids).size, ids.length, 'guide ids are unique');

    for (const g of all) {
      expect(g.steps.length >= 1, true, `${g.id} has steps`);
      expect(!!g.title && !!g.summary, true, `${g.id} has title + summary`);
      for (const s of g.steps)
        expect(!!s.title && !!s.text, true, `${g.id} step "${s.title}" has title + text`);
    }
    // Tutorials are multi-step; questions are single-answer.
    for (const t of TUTORIALS) expect(t.steps.length >= 3, true, `${t.id} is multi-step`);
    // The Start-panel tour launches the first tutorial: keep it the hands-on one.
    expect(TUTORIALS[0].id, 'load-data-add-column');
  });
});

/** Drive the flagship tutorial's core against a live editor + property panel:
 *  the same DOM the guide resolvers query at runtime. */
category('Flow: guide playthrough', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('OpenFile node + fullPath field are targetable and the typed-value wait fires', async () => {
    const reg = getRegisteredFuncs().find((r) => r.func.name === 'OpenFile');
    if (!reg) {
      expect(true, true, 'OpenFile not on this stand — skipped');
      return;
    }
    const e = makeEditor();
    const ac = new AbortController();
    const host: GuideHost = {getFlow: () => e.flow, showFunctionBrowser: () => {}, anchorEl: e.container};
    const ctx: GuideContext = {host, signal: ac.signal};
    try {
      const node = createNode(reg.nodeTypeName)!;
      await e.flow.addNodeAt(node, 0, 0);
      // byNodeFunc('OpenFile') finds the rendered canvas node.
      const ok = await until(() => byNodeFunc('OpenFile')(ctx) !== null, 2000);
      expect(ok, true, 'OpenFile node renders and is targetable by data-func');
      const nodeEl = byNodeFunc('OpenFile')(ctx)!;
      expect(nodeEl.dataset.func, 'OpenFile');

      // Property panel: the fullPath row is targetable, and its editable field
      // is reachable via paramFieldSelector. Mount the panel root in the
      // document (in the app it's grok.shell.o) so the document-wide resolvers
      // can find it.
      const panel = new PropertyPanel(e.flow);
      panel.root.style.cssText = 'position:absolute;left:-10000px;';
      document.body.appendChild(panel.root);
      panel.showNode(node);
      const paramOk = await until(() => byParam('fullPath')(ctx) !== null, 1500);
      expect(paramOk, true, 'fullPath row targetable by data-param');
      const field = document.querySelector(paramFieldSelector('fullPath')) as HTMLTextAreaElement | null;
      expect(!!field, true, 'fullPath editable field found');

      // Simulate the paste and confirm the guide's "until value contains" resolves.
      field!.value = 'System:DemoFiles/demog.csv';
      field!.dispatchEvent(new Event('input', {bubbles: true}));
      await untilValueContains(paramFieldSelector('fullPath'), 'demog.csv')(ctx);
      expect(true, true, 'typed-value wait resolved');
      panel.root.remove();
    } finally {
      ac.abort();
      destroyEditor(e);
    }
  });
});
