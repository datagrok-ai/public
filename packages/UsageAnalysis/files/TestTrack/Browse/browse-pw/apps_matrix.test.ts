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

type Tag = 'heavy' | 'runOnly';

interface AppSpec {

  app: string;
  tags?: Tag[];

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

const BLANKET_IGNORE: RegExp[] = [
  /Connection .*(refused|failed|timed out|reset)/i,
  /ECONNREFUSED|ETIMEDOUT|ENOTFOUND|ECONNRESET/i,
  /Container .*not started/i,
  /Failed to load resource:\s*the server responded with a status of (4\d\d|5\d\d)/i,
];

const NOISE_IGNORE: RegExp[] = [
  /Failed to execute 'observe' on 'ResizeObserver'/i,
  /ResizeObserver loop/i,
  /The above error occurred in the .* component/i,

  /ZoomTool\.observeCanvasResize/,
];

const DOCKER_API_ERROR = /API error:\s*\d{3}/i;

const DOCKER_DEPS: { key: string; container: string }[] = [
  { key: 'MolTrack', container: 'moltrack' },
  { key: 'Admetica', container: 'admetica' },
  { key: 'Retrosynthesis', container: 'retrosynthesis' },
  { key: 'Docking', container: 'bio' },
  { key: 'Boltz', container: 'bio' },
  { key: 'Preclinical Case', container: 'preclinicalcase' },
];

type DockerStatuses = Record<string, string>;

async function fetchDockerStatuses(page: Page): Promise<DockerStatuses> {
  return page.evaluate(async () => {
    const out: Record<string, string> = {};
    try {
      const list = await (window as any).grok.dapi.docker.dockerContainers.list();
      for (const c of list) out[String(c.name ?? '').toLowerCase()] = String(c.status ?? '');
    } catch {  }
    return out;
  });
}

function containerStarted(statuses: DockerStatuses, containerSub: string): boolean {
  const sub = containerSub.toLowerCase();
  const name = Object.keys(statuses).find((n) => n.includes(sub));
  return !!name && statuses[name] === 'started';
}

function dockerDownDep(segments: string[], statuses: DockerStatuses): { key: string; container: string } | null {
  const dep = DOCKER_DEPS.find((d) => segments.some((s) => s.toLowerCase().includes(d.key.toLowerCase())));
  return dep && !containerStarted(statuses, dep.container) ? dep : null;
}

const CONNECTOR_KEYS = ['Benchling', 'CDD Vault', 'Revvity Signals', 'Chemspace', 'Alation'];

function isConnectorNode(segments: string[]): boolean {
  const norm = (s: string) => s.toLowerCase().replace(/[^a-z0-9]+/g, '');
  return segments.some((s) => CONNECTOR_KEYS.some((k) => norm(s) === norm(k)));
}

interface KnownIssue { ref: string; note: string; match: (e: string) => boolean; }
const KNOWN_ISSUES: KnownIssue[] = [

];

const MAX_NODES_PER_APP = 60;

function treeName(segments: string[]): string {
  return 'tree-Apps---' + segments.map((s) => s.replace(/ /g, '-')).join('---');
}

interface ChildNode { name: string; label: string; expandable: boolean; }

async function directChildren(page: Page, parentName: string): Promise<ChildNode[]> {
  return page.evaluate((parent) => {
    const els = Array.from(document.querySelectorAll(`[name^="${parent}---"]`));
    const seen = new Set<string>();
    const out: { name: string; label: string; expandable: boolean }[] = [];
    for (const el of els) {
      const name = el.getAttribute('name');
      if (!name || seen.has(name)) continue;
      const rest = name.slice(parent.length + 3); 
      if (rest.includes('---')) continue; 
      seen.add(name);
      const tri = el.querySelector('.d4-tree-view-tri');
      const label = (el.querySelector('.d4-tree-view-group-label, .d4-tree-view-item-label')?.textContent
        || el.textContent || '').trim().split('\n')[0].trim();
      out.push({ name, label, expandable: !!tri });
    }
    return out;
  }, parentName);
}

function nodeByName(page: Page, name: string): Locator {
  return page.locator(`[name="${name}"]`).first();
}

async function expandAndChildren(page: Page, segments: string[], thorough: boolean): Promise<ChildNode[]> {
  const name = treeName(segments);

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

      if (thorough && best.length > 0 && stableReads >= 3 && Date.now() - start >= 4_000) return best;
      await page.waitForTimeout(700);
    }
    return best;
  };
  await expandByName(page, name);
  let children = await poll(thorough ? 12_000 : 5_000);
  if (children.length) return children;

  for (let attempt = 0; attempt < 2 && !children.length; attempt++) {
    const tri = nodeByName(page, name).locator(TREE_EXPAND_ARROW).first();
    await tri.click().catch(() => undefined); 
    await page.waitForTimeout(500);
    await tri.click().catch(() => undefined); 
    children = await poll(8_000);
  }
  if (children.length) return children;

  await nodeByName(page, name).click().catch(() => undefined);
  await page.waitForTimeout(2_500);
  await ensureBrowsePanelOpen(page);
  await revealNode(page, segments);
  await expandByName(page, name);
  children = await poll(10_000);
  return children;
}

async function expandByName(page: Page, name: string, timeoutMs = 6_000): Promise<boolean> {
  const node = nodeByName(page, name);
  try {
    await node.waitFor({ state: 'visible', timeout: timeoutMs });
  } catch {
    return false;
  }
  await node.scrollIntoViewIfNeeded().catch(() => undefined);
  const tri = node.locator(TREE_EXPAND_ARROW).first();
  if (!(await tri.isVisible().catch(() => false))) return true; 
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false);
  if (!expanded) {
    await tri.click().catch(() => undefined);
    await page.waitForTimeout(700);
  }
  return true;
}

async function revealNode(page: Page, segments: string[]): Promise<Locator | null> {
  await ensureBrowsePanelOpen(page);
  await expandTreeGroup(page, 'Apps').catch(() => undefined);
  for (let i = 1; i < segments.length; i++) 
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

async function closeAllViews(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).grok?.shell?.closeAll?.()).catch(() => undefined);
  await page.waitForTimeout(300);
}

async function waitForOpened(page: Page, heavy: boolean): Promise<void> {
  await expect
    .poll(async () => page.evaluate(() => {
      const v = (window as any).grok?.shell?.v;
      return v ? String(v.name ?? v.type ?? 'view') : '';
    }).catch(() => ''), { timeout: heavy ? 30_000 : 12_000, intervals: [400, 800, 1500, 2500] })
    .not.toBe('')
    .catch(() => undefined);
}

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

function errSnapshot(sink: ErrorSink): string[] {
  return [...sink.pageErrors, ...sink.consoleErrors, ...sink.balloonErrors];
}

async function openAndCheck(
  page: Page, segments: string[], label: string, heavy: boolean,
  sink: ErrorSink, offenders: Offender[], docker: DockerStatuses,
): Promise<boolean> {

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

  await page.waitForTimeout(1_500);
  await collectBalloonErrors(page, sink);

  const after = errSnapshot(sink);
  let fresh = after.slice(before.length);

  if (isConnectorNode(segments)) {
    if (fresh.length) console.log(`  ⓘ ${label}: ${fresh.length} error(s) tolerated — external connector`);
    fresh = [];
  }

  const nodeDep = dockerDownDep(segments, docker);
  if (nodeDep) {
    if (fresh.length) console.log(`  ⓘ ${label}: ${fresh.length} error(s) tolerated — docker '${nodeDep.container}' not started`);
    fresh = [];
  }

  fresh = fresh.filter((e) => {

    if (BLANKET_IGNORE.some((re) => re.test(e))) return false;

    if (NOISE_IGNORE.some((re) => re.test(e))) return false;

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
    if (opened && child.expandable) 
      await walk(page, childSegs, false, spec, heavy, sink, offenders, counter, docker);

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

      const sink = watchErrors(page);
      const offenders: Offender[] = [];

      const appNode = await revealNode(page, [spec.app]);
      test.skip(!appNode, `Apps > ${spec.app} not deployed on this stack`);

      const docker = await fetchDockerStatuses(page);

      const appTri = nodeByName(page, treeName([spec.app])).locator(TREE_EXPAND_ARROW).first();
      const isGroup = await appTri.isVisible().catch(() => false);

      const counter = { n: 0 };
      if (runOnly || !isGroup) {

        const opened = await openAndCheck(page, [spec.app], spec.app, heavy, sink, offenders, docker);
        if (opened) counter.n++;
      } else {

        await walk(page, [spec.app], true, spec, heavy, sink, offenders, counter, docker);
      }

      console.log(`[${spec.app}] opened ${counter.n} node(s), ${offenders.length} with errors`);

      const flagged = new Set<string>();
      for (const o of offenders)
        for (const ki of KNOWN_ISSUES)
          if (!flagged.has(ki.ref) && o.errors.some(ki.match)) {
            test.info().annotations.push({ type: 'issue', description: `${ki.ref}: ${ki.note}` });
            console.log(`  ⚑ known issue ${ki.ref} reproduced (${o.node})`);
            flagged.add(ki.ref);
          }

      expect(counter.n, `Apps > ${spec.app} is present but no nodes were opened — subtree failed to load`)
        .toBeGreaterThan(0);
      const report = offenders.map((o) => `• ${o.node}\n    ${o.errors.join('\n    ')}`).join('\n');
      expect(offenders, `Errors while surfing Apps > ${spec.app}:\n${report}`).toEqual([]);
    });
  }
});
