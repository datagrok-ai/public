import { test, expect, Page, Locator } from '@playwright/test';
import { CONTEXT_PANEL, TREE_EXPAND_ARROW } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  expandTreeGroup,
  watchErrors,
  collectBalloonErrors,
  ErrorSink,
} from './helpers';

/**
 * Browse > Apps — open EVERY installed app in the Browse tree, walk into the apps that
 * expand into a data subtree, open the views their nodes produce, interact with each view,
 * and assert there are no errors anywhere while surfing:
 *   - no uncaught page errors (`pageerror`)
 *   - no `console.error`
 *   - no error balloons
 *   - no error raised while the Context Panel renders all its panes
 *
 * Scope: every top-level node under `Apps` EXCEPT `Demo` (covered by demo_apps.test.ts) and
 * `Compute > Model Hub` (covered by modelhub.test.ts).
 *
 * Coverage shape ("sample" — keeps the suite runnable under the CI single-worker timeout):
 *   - leaf app  -> open it, interact, check.
 *   - group app -> open ALL of its top-level children (the app's tools / instances), and for
 *     every child that itself expands, drill ONE representative chain (first child at each
 *     level, down to a leaf). This opens one representative data-instance deep per branch
 *     instead of the full combinatorial data explosion (e.g. Preclinical Case = 216 nodes).
 *
 * Portability: the matrix is captured on dev (full app set). On a minimal stack a plugin's
 * package isn't installed, so the app's tree node is simply absent -> `test.skip()` for that
 * whole app, and absent child nodes are skipped individually. The same file is therefore
 * correct on dev and on a reduced CI / customer stack with no hardcoded allowlist to drift.
 */

type Tag = 'heavy' | 'runOnly';

interface AppSpec {
  /** Top-level node label under `Apps` (exactly as shown in the tree). */
  app: string;
  tags?: Tag[];
  /** Child labels (direct or any descendant segment) to skip — e.g. `Model Hub`. */
  exclude?: string[];
}

const APPS: AppSpec[] = [
  { app: 'Bio', tags: ['heavy'] },
  { app: 'Chem', tags: ['heavy'] },
  { app: 'Clinical Case', tags: ['heavy'] },
  { app: 'Compute', tags: ['heavy'], exclude: ['Model Hub'] },
  { app: 'Flow' },
  { app: 'Misc' },
  { app: 'Peptides', tags: ['heavy'] },
  { app: 'Plates' },
  { app: 'PopPK', tags: ['heavy'] },
  { app: 'Preclinical Case', tags: ['heavy'] },
  { app: 'Tutorials', tags: ['runOnly'] },
];

/**
 * Errors that are NEVER a Browse/UI defect, dropped on EVERY node (they fire async and bleed
 * across nodes): connectivity drops, the platform explicitly reporting a docker container is not
 * up, and raw backend HTTP responses. Specific signatures so a real UI bug whose message merely
 * contains a number is not masked.
 */
const BLANKET_IGNORE: RegExp[] = [
  /Connection .*(refused|failed|timed out|reset)/i,
  /ECONNREFUSED|ETIMEDOUT|ENOTFOUND|ECONNRESET/i,
  /Container .*not started/i,
  /Failed to load resource:\s*the server responded with a status of (4\d\d|5\d\d)/i,
];

/**
 * Benign UI noise provoked by our own poking (clicking / hovering a React-based editor while it
 * is mounting), NOT a real navigation defect — confirmed by opening the same views patiently in
 * isolation (zero errors). The ResizeObserver "observe non-element" warning is a well-known React
 * mount/resize race; the "above error occurred in <…> component" line is React's error-boundary
 * echo of it. Dropped on every node.
 */
const NOISE_IGNORE: RegExp[] = [
  /Failed to execute 'observe' on 'ResizeObserver'/i,
  /ResizeObserver loop/i,
  /The above error occurred in the .* component/i,
];

/** A thrown error from a docker-backed app's API layer (e.g. "MolTrack API error: 400 ..."). */
const DOCKER_API_ERROR = /API error:\s*\d{3}/i;

/**
 * Docker-backed apps: node-label / error-text keyword -> docker container name (substring).
 * When such an error fires we look up the LIVE container status: `started` => the app genuinely
 * failed (real bug, keep it); anything else (`starting`/`error`/`stopped`/absent) => the backend
 * isn't available, which is infrastructure, not a Browse/UI defect => tolerate.
 */
const DOCKER_DEPS: { key: string; container: string }[] = [
  { key: 'MolTrack', container: 'moltrack' },
  { key: 'Admetica', container: 'admetica' },
  { key: 'Retrosynthesis', container: 'retrosynthesis' },
  { key: 'Docking', container: 'bio' },
  { key: 'Boltz', container: 'bio' },
  { key: 'Preclinical Case', container: 'preclinicalcase' },
];

/** Snapshot of docker container statuses (`name` lowercased -> status), taken once per app. */
type DockerStatuses = Record<string, string>;

async function fetchDockerStatuses(page: Page): Promise<DockerStatuses> {
  return page.evaluate(async () => {
    const out: Record<string, string> = {};
    try {
      const list = await (window as any).grok.dapi.docker.dockerContainers.list();
      for (const c of list) out[String(c.name ?? '').toLowerCase()] = String(c.status ?? '');
    } catch { /* leave empty -> unknown */ }
    return out;
  });
}

function containerStarted(statuses: DockerStatuses, containerSub: string): boolean {
  const sub = containerSub.toLowerCase();
  const name = Object.keys(statuses).find((n) => n.includes(sub));
  return !!name && statuses[name] === 'started';
}

/** Return the docker dependency of this node IF its container is NOT started (else null). */
function dockerDownDep(segments: string[], statuses: DockerStatuses): { key: string; container: string } | null {
  const dep = DOCKER_DEPS.find((d) => segments.some((s) => s.toLowerCase().includes(d.key.toLowerCase())));
  return dep && !containerStarted(statuses, dep.container) ? dep : null;
}

/**
 * External-system connectors (SaaS): the app's view depends entirely on a third-party service that
 * is configured only on dev and absent elsewhere. Its errors are environment-dependent, not a
 * Browse/UI defect, so all errors on these nodes are tolerated (user decision 2026-06-23).
 */
const CONNECTOR_KEYS = ['Benchling', 'CDD Vault', 'Revvity Signals', 'Chemspace', 'Alation'];

function isConnectorNode(segments: string[]): boolean {
  const norm = (s: string) => s.toLowerCase().replace(/[^a-z0-9]+/g, '');
  return segments.some((s) => CONNECTOR_KEYS.some((k) => norm(s) === norm(k)));
}

/**
 * Known, filed platform/plugin bugs surfaced by this suite. We DO NOT suppress them — the test
 * goes red honestly whenever the bug fires; we only attach the ticket to the report so reviewers
 * can attribute the failure. When the bug is fixed the error stops, the test goes green, and the
 * entry can be removed. (Same philosophy as KNOWN_BUGS.md / demo_apps.)
 */
interface KnownIssue { ref: string; note: string; match: (e: string) => boolean; }
const KNOWN_ISSUES: KnownIssue[] = [
  // GROK-20254 (Hit Triage / Hit Design recurring null.id) — fixed on dev, entry removed.
];

/** Safety cap on opened nodes per app — surfing must stay within the per-test timeout. */
const MAX_NODES_PER_APP = 60;

/** `['Chem','Reactions']` -> `tree-Apps---Chem---Reactions` (spaces -> dashes, other chars kept). */
function treeName(segments: string[]): string {
  return 'tree-Apps---' + segments.map((s) => s.replace(/ /g, '-')).join('---');
}

interface ChildNode { name: string; label: string; expandable: boolean; }

/** Direct children of a tree node (must be expanded first). Reads stable `name=` attributes. */
async function directChildren(page: Page, parentName: string): Promise<ChildNode[]> {
  return page.evaluate((parent) => {
    const els = Array.from(document.querySelectorAll(`[name^="${parent}---"]`));
    const seen = new Set<string>();
    const out: { name: string; label: string; expandable: boolean }[] = [];
    for (const el of els) {
      const name = el.getAttribute('name');
      if (!name || seen.has(name)) continue;
      const rest = name.slice(parent.length + 3); // strip "parent---"
      if (rest.includes('---')) continue; // not a direct child
      seen.add(name);
      const tri = el.querySelector('.d4-tree-view-tri');
      const label = (el.querySelector('.d4-tree-view-group-label, .d4-tree-view-item-label')?.textContent
        || el.textContent || '').trim().split('\n')[0].trim();
      out.push({ name, label, expandable: !!tri });
    }
    return out;
  }, parentName);
}

/** Locator for a tree node by its stable `name` attribute. */
function nodeByName(page: Page, name: string): Locator {
  return page.locator(`[name="${name}"]`).first();
}

/**
 * Expand a node and return its direct children, polling for lazy-loaded subtrees. Some apps
 * build their Browse subtree asynchronously (an `@appTreeBrowser` that hits the DB) so a single
 * triangle click + fixed wait misses the children. Strategy, escalating only as needed:
 *   1. expand the triangle, poll.
 *   2. collapse + re-expand (some trees populate on a second expand), poll.
 *   3. open the app node to force it to build its tree, then RE-OPEN Browse (the click may
 *      have swapped the panel for Toolbox), re-reveal the path, expand, poll.
 */
async function expandAndChildren(page: Page, segments: string[], thorough: boolean): Promise<ChildNode[]> {
  const name = treeName(segments);
  // Data apps (Clinical/Preclinical Case studies, etc.) STREAM their children in over several
  // seconds (and a static node like "Import study" appears first), so returning on the first
  // non-empty read under-covers. When we need the FULL set (top-level, `thorough`), poll the
  // whole window and keep the LARGEST set seen. When we only need a representative child (deep
  // drill), return as soon as anything appears.
  const poll = async (ms: number): Promise<ChildNode[]> => {
    const deadline = Date.now() + ms;
    const start = Date.now();
    let best: ChildNode[] = [];
    let stableReads = 0;
    while (Date.now() < deadline) {
      const cur = await directChildren(page, name);
      if (cur.length > best.length) { best = cur; stableReads = 0; }
      else stableReads++;
      if (!thorough && best.length > 0) return best;
      // Thorough: exit early once the count has stopped growing (children streamed in) instead of
      // burning the whole window. Require a small minimum elapsed so a static first node (e.g.
      // "Import study") doesn't look "stable" before the real instances start streaming.
      if (thorough && best.length > 0 && stableReads >= 3 && Date.now() - start >= 4_000) return best;
      await page.waitForTimeout(700);
    }
    return best;
  };
  await expandByName(page, name);
  let children = await poll(thorough ? 12_000 : 5_000);
  if (children.length) return children;

  // Collapse + re-expand toggle.
  for (let attempt = 0; attempt < 2 && !children.length; attempt++) {
    const tri = nodeByName(page, name).locator(TREE_EXPAND_ARROW).first();
    await tri.click().catch(() => undefined); // collapse
    await page.waitForTimeout(500);
    await tri.click().catch(() => undefined); // expand
    children = await poll(8_000);
  }
  if (children.length) return children;

  // Last resort: open the app to build its subtree, then re-open Browse and re-expand.
  await nodeByName(page, name).click().catch(() => undefined);
  await page.waitForTimeout(2_500);
  await ensureBrowsePanelOpen(page);
  await revealNode(page, segments);
  await expandByName(page, name);
  children = await poll(10_000);
  return children;
}

/** Expand a node (by name) if it is collapsible and not yet expanded. No-op for leaves. */
async function expandByName(page: Page, name: string, timeoutMs = 6_000): Promise<boolean> {
  const node = nodeByName(page, name);
  try {
    await node.waitFor({ state: 'visible', timeout: timeoutMs });
  } catch {
    return false;
  }
  await node.scrollIntoViewIfNeeded().catch(() => undefined);
  const tri = node.locator(TREE_EXPAND_ARROW).first();
  if (!(await tri.isVisible().catch(() => false))) return true; // leaf
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false);
  if (!expanded) {
    await tri.click().catch(() => undefined);
    await page.waitForTimeout(700);
  }
  return true;
}

/**
 * Make a node reachable: ensure Browse is open, expand `Apps` and every ancestor segment,
 * then return the node locator — or null if any node along the path is absent (plugin not
 * installed on this stack).
 */
async function revealNode(page: Page, segments: string[]): Promise<Locator | null> {
  await ensureBrowsePanelOpen(page);
  await expandTreeGroup(page, 'Apps').catch(() => undefined);
  for (let i = 1; i < segments.length; i++) // expand every ancestor (not the leaf itself)
    if (!await expandByName(page, treeName(segments.slice(0, i)))) return null;
  const leaf = nodeByName(page, treeName(segments));
  try {
    await leaf.waitFor({ state: 'visible', timeout: 6_000 });
  } catch {
    return null;
  }
  await leaf.scrollIntoViewIfNeeded().catch(() => undefined);
  return leaf;
}

/** Close all open views for a clean slate (the platform piles up restored views across opens). */
async function closeAllViews(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).grok?.shell?.closeAll?.()).catch(() => undefined);
  await page.waitForTimeout(300);
}

/** Best-effort wait for some non-Home view to be current after clicking a node. */
async function waitForOpened(page: Page, heavy: boolean): Promise<void> {
  await expect
    .poll(async () => page.evaluate(() => {
      const v = (window as any).grok?.shell?.v;
      return v ? String(v.name ?? v.type ?? 'view') : '';
    }).catch(() => ''), { timeout: heavy ? 30_000 : 12_000, intervals: [400, 800, 1500, 2500] })
    .not.toBe('')
    .catch(() => undefined);
}

/** "Poke" the opened view to exercise render / event handlers. Best-effort, asserts nothing. */
async function pokeView(page: Page): Promise<void> {
  const doc = page.locator('.d4-root .grok-view, .layout-workarea .d4-view, .document-manager .tab-content').first();
  if (await doc.isVisible().catch(() => false)) {
    const box = await doc.boundingBox().catch(() => null);
    if (box) await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2).catch(() => undefined);
  }
  const gridCanvas = page.locator('canvas.d4-grid, .d4-grid canvas').first();
  if (await gridCanvas.isVisible().catch(() => false)) {
    const gb = await gridCanvas.boundingBox().catch(() => null);
    if (gb) await page.mouse.click(gb.x + Math.min(60, gb.width / 2), gb.y + Math.min(80, gb.height / 2)).catch(() => undefined);
  }
  const viewer = page.locator('.d4-viewer:not(.d4-grid)').first();
  if (await viewer.isVisible().catch(() => false)) await viewer.hover().catch(() => undefined);
  await page.waitForTimeout(600);
}

interface Offender { node: string; errors: string[]; }

/** Total errors currently recorded across all three channels. */
function errSnapshot(sink: ErrorSink): string[] {
  return [...sink.pageErrors, ...sink.consoleErrors, ...sink.balloonErrors];
}

/**
 * Open one node, interact with the resulting view, render the Context Panel, and record any
 * NEW error that appeared while doing so (attributed to this node). Returns the running count
 * of nodes opened (for the per-app cap).
 */
async function openAndCheck(
  page: Page, segments: string[], label: string, heavy: boolean,
  sink: ErrorSink, offenders: Offender[], docker: DockerStatuses,
): Promise<boolean> {
  // NB: we do NOT closeAll between sibling/child nodes — only once per top-level app (in
  // beforeEach's goHome reload). Closing a view mid-walk while its debounced async search is
  // still in flight makes the late callback fire on a disposed view and throw — a test-induced
  // race, not a real defect. Keeping prior views open lets that callback land safely.
  const node = await revealNode(page, segments);
  if (!node) {
    console.log(`  · absent: ${label}`);
    return false;
  }
  await collectBalloonErrors(page, sink);
  const before = errSnapshot(sink);

  await node.click().catch(() => undefined);
  await waitForOpened(page, heavy);
  await pokeView(page);
  await ensureContextPanelOpen(page, true);
  await expect(page.locator(CONTEXT_PANEL), 'Context Panel should be visible')
    .toBeVisible({ timeout: 10_000 }).catch(() => undefined);
  // Brief settle so a view's debounced async render finishes and any error it raises is recorded
  // before we judge the node. We deliberately do NOT wait for `networkidle`: Datagrok polls
  // continuously (it rarely goes idle), so that wait usually burned its full timeout for nothing.
  // A short fixed settle is predictable and far faster; correctness is unaffected because the
  // per-app pass/fail is the SUM of errors, not their exact per-node attribution.
  await page.waitForTimeout(1_500);
  await collectBalloonErrors(page, sink);

  const after = errSnapshot(sink);
  let fresh = after.slice(before.length);

  // External-system connector node: its view depends on a third-party SaaS — tolerate all errors.
  if (isConnectorNode(segments)) {
    if (fresh.length) console.log(`  ⓘ ${label}: ${fresh.length} error(s) tolerated — external connector`);
    fresh = [];
  }

  // Docker-backed node whose container isn't `started`: the whole view is non-functional for an
  // infrastructure reason (not a Browse/UI defect) — tolerate all errors. (User rule: docker down
  // => skip; docker up => errors are real and stay.)
  const nodeDep = dockerDownDep(segments, docker);
  if (nodeDep) {
    if (fresh.length) console.log(`  ⓘ ${label}: ${fresh.length} error(s) tolerated — docker '${nodeDep.container}' not started`);
    fresh = [];
  }

  fresh = fresh.filter((e) => {
    // Pure infra (connectivity, container-not-started, raw HTTP) — never a Browse/UI defect.
    if (BLANKET_IGNORE.some((re) => re.test(e))) return false;
    // Benign poking-induced UI noise (ResizeObserver / React boundary echo).
    if (NOISE_IGNORE.some((re) => re.test(e))) return false;
    // Docker-backed app API error that bled onto a NON-docker node (matched by error text): honest
    // only if that app's container is up; otherwise it is the backend being unavailable.
    if (DOCKER_API_ERROR.test(e)) {
      const dep = DOCKER_DEPS.find((d) => e.toLowerCase().includes(d.key.toLowerCase()));
      if (dep && !containerStarted(docker, dep.container)) return false;
    }
    return true;
  });
  if (fresh.length) {
    offenders.push({ node: label, errors: fresh });
    console.log(`  ✘ ${label}: ${fresh.length} new error(s)`);
  } else
    console.log(`  ✓ ${label}`);
  return true;
}

/**
 * Walk a subtree. `openAllChildren=true` at the app's top level opens every child; deeper
 * levels follow only the first child (representative chain) so one data-instance is opened
 * deep without the full data explosion. `counter` enforces the per-app node cap.
 */
async function walk(
  page: Page, segments: string[], openAllChildren: boolean, spec: AppSpec, heavy: boolean,
  sink: ErrorSink, offenders: Offender[], counter: { n: number }, docker: DockerStatuses,
): Promise<void> {
  if (counter.n >= MAX_NODES_PER_APP) return;
  await revealNode(page, segments);
  let children = await expandAndChildren(page, segments, openAllChildren);
  console.log(`  [walk] ${segments.join('/')} -> ${children.length} child(ren): ${children.map((c) => c.label).join(', ')}`);
  if (spec.exclude)
    children = children.filter((c) => !spec.exclude!.some((e) => c.label === e || c.name.includes('---' + e.replace(/ /g, '-'))));
  if (children.length === 0) return;

  const toOpen = openAllChildren ? children : children.slice(0, 1);
  for (const child of toOpen) {
    if (counter.n >= MAX_NODES_PER_APP) {
      console.log(`  … node cap (${MAX_NODES_PER_APP}) reached — remaining ${segments.join('/')} children not opened`);
      return;
    }
    const childSegs = [...segments, child.label];
    const opened = await openAndCheck(page, childSegs, childSegs.slice(1).join(' / '), heavy, sink, offenders, docker);
    if (opened) counter.n++;
    if (opened && child.expandable) // drill one representative chain deeper
      await walk(page, childSegs, false, spec, heavy, sink, offenders, counter, docker);
    // A connector / docker-down sub-app keeps emitting error balloons in the background; those
    // bleed onto the next sibling node and red it falsely. Close its views to silence it before
    // moving on. (Targeted exception to "closeAll only between top-level apps".)
    if (opened && (isConnectorNode(childSegs) || dockerDownDep(childSegs, docker)))
      await closeAllViews(page);
  }
}

test.describe('Browse Apps matrix (Browse-AppsMatrix-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await closeAllViews(page);
    await ensureBrowsePanelOpen(page);
  });

  for (const spec of APPS) {
    const heavy = !!spec.tags?.includes('heavy');
    const runOnly = !!spec.tags?.includes('runOnly');

    test(`Browse-AppsMatrix — ${spec.app}`, async ({ page }) => {
      test.setTimeout(runOnly ? 90_000 : heavy ? 420_000 : 240_000);
      // Generic browser noise is filtered in helpers; infra/docker errors are filtered per-node
      // inside openAndCheck (BLANKET_IGNORE + live docker-status check), not globally here.
      const sink = watchErrors(page);
      const offenders: Offender[] = [];

      // App must be present on this stack; otherwise skip the whole app (plugin not installed).
      const appNode = await revealNode(page, [spec.app]);
      test.skip(!appNode, `Apps > ${spec.app} not deployed on this stack`);

      // Live docker container statuses: a docker-backed app whose container isn't `started` is an
      // infra condition (tolerated), not an app bug. Taken once per app.
      const docker = await fetchDockerStatuses(page);

      const appTri = nodeByName(page, treeName([spec.app])).locator(TREE_EXPAND_ARROW).first();
      const isGroup = await appTri.isVisible().catch(() => false);

      const counter = { n: 0 };
      if (runOnly || !isGroup) {
        // Leaf app (or run-only): open the app's own view and check.
        const opened = await openAndCheck(page, [spec.app], spec.app, heavy, sink, offenders, docker);
        if (opened) counter.n++;
      } else {
        // Group app: open all top-level children; drill one representative chain per branch.
        await walk(page, [spec.app], true, spec, heavy, sink, offenders, counter, docker);
      }

      console.log(`[${spec.app}] opened ${counter.n} node(s), ${offenders.length} with errors`);

      // Attribute any known, filed bug to the report (WITHOUT suppressing it — the test still
      // goes red). When the platform fix lands the error disappears and the test goes green.
      const flagged = new Set<string>();
      for (const o of offenders)
        for (const ki of KNOWN_ISSUES)
          if (!flagged.has(ki.ref) && o.errors.some(ki.match)) {
            test.info().annotations.push({ type: 'issue', description: `${ki.ref}: ${ki.note}` });
            console.log(`  ⚑ known issue ${ki.ref} reproduced (${o.node})`);
            flagged.add(ki.ref);
          }

      // Guard against a false-green: a group app whose subtree silently failed to load would
      // open 0 nodes and trivially "pass". The app node IS present (we passed the skip), so
      // an empty walk is a load failure, not a clean run.
      expect(counter.n, `Apps > ${spec.app} is present but no nodes were opened — subtree failed to load`)
        .toBeGreaterThan(0);
      const report = offenders.map((o) => `• ${o.node}\n    ${o.errors.join('\n    ')}`).join('\n');
      expect(offenders, `Errors while surfing Apps > ${spec.app}:\n${report}`).toEqual([]);
    });
  }
});
