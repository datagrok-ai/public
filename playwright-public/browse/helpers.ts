import { Page, expect, Locator } from '@playwright/test';
import {
  BROWSE_HEADER,
  BROWSE_HEADER_HOME,
  BROWSE_HEADER_COLLAPSE_ALL,
  BROWSE_PANEL_TAB_PANE,
  BALLOON_CONTAINER,
  CONTEXT_PANEL,
  CONTEXT_PANEL_EXPAND_ALL,
  RIBBON,
  SIDEBAR_BROWSE_ICON,
  TREE_EXPAND_ARROW,
  TREE_EXPAND_ARROW_EXPANDED,
  treeGroupByName,
} from './selectors';

export const BASE: string = process.env.DATAGROK_URL!;

/** Recording sink for errors observed during a test scenario. */
export interface ErrorSink {
  pageErrors: string[];
  consoleErrors: string[];
  balloonErrors: string[];
}

/**
 * Console messages we deliberately ignore for Browse tests — they are browser-level
 * noise unrelated to app behavior (missing resources, ad blockers, etc).
 */
const IGNORED_CONSOLE_PATTERNS: RegExp[] = [
  // Generic browser-level resource noise (missing icons, lazy-loaded assets, etc.).
  /Failed to load resource:\s*the server responded with a status of 404/i,
  /Failed to load resource:\s*net::ERR_/i,
  // Aborted fetches during fast page navigation: fires when an in-flight request is
  // cancelled because the user navigated away. Not a real failure.
  /TypeError:\s*Failed to fetch/i,
  /Error fetching .*:\s*TypeError:\s*Failed to fetch/i,
  // Known platform noise from the global search "similarity" subroutine: fires on every
  // keystroke, does not affect search results. Track separately if it becomes blocking.
  /Error during similarity search:/i,
  // External images embedded in demo descriptions / Markup (e.g. Wikipedia) are CORS-blocked
  // by the browser. This is demo-content noise, not a platform defect.
  /Access to image at .* has been blocked by CORS policy/i,
];

/**
 * Attach console.error / pageerror listeners and start watching for error balloons.
 * Call once per test; pass the returned sink to `expectNoErrors` at the end (or wherever needed).
 */
export function watchErrors(page: Page, extraIgnore: RegExp[] = []): ErrorSink {
  const sink: ErrorSink = { pageErrors: [], consoleErrors: [], balloonErrors: [] };
  const ignored = (text: string): boolean =>
    IGNORED_CONSOLE_PATTERNS.some((re) => re.test(text)) || extraIgnore.some((re) => re.test(text));
  page.on('pageerror', (err) => {
    if (!ignored(String(err))) sink.pageErrors.push(String(err));
  });
  page.on('console', (msg) => {
    if (msg.type() !== 'error') return;
    const text = msg.text();
    if (ignored(text)) return;
    sink.consoleErrors.push(text);
  });
  return sink;
}

/** Snapshot current visible error balloons and merge into the sink. */
export async function collectBalloonErrors(page: Page, sink: ErrorSink): Promise<void> {
  const texts = await page
    .locator(`${BALLOON_CONTAINER} .d4-balloon-error, ${BALLOON_CONTAINER} .grok-balloon-error`)
    .allTextContents()
    .catch(() => [] as string[]);
  for (const t of texts) {
    const trimmed = t.trim();
    if (trimmed && !sink.balloonErrors.includes(trimmed)) sink.balloonErrors.push(trimmed);
  }
}

/** Assert no errors were collected. Convenience for `assertNoErrors(page)` calls. */
export async function expectNoErrors(page: Page, sink: ErrorSink): Promise<void> {
  await collectBalloonErrors(page, sink);
  expect(
    sink.pageErrors,
    `pageerror occurred:\n${sink.pageErrors.join('\n')}`,
  ).toEqual([]);
  expect(
    sink.consoleErrors,
    `console.error occurred:\n${sink.consoleErrors.join('\n')}`,
  ).toEqual([]);
  expect(
    sink.balloonErrors,
    `error balloon shown:\n${sink.balloonErrors.join('\n')}`,
  ).toEqual([]);
}

/**
 * Open the app, wait until the ribbon is visible (= app loaded).
 * Uses the storageState from global-setup, so no login is performed here.
 */
export async function goHome(page: Page): Promise<void> {
  await page.goto(BASE);
  await page.waitForSelector(RIBBON, { timeout: 30_000 });
}

/** Ensure the Browse side panel is visible. Idempotent. */
export async function ensureBrowsePanelOpen(page: Page): Promise<void> {
  const header = page.locator(BROWSE_HEADER);
  if (!(await header.isVisible().catch(() => false))) {
    await page.locator(SIDEBAR_BROWSE_ICON).click();
    await header.waitFor({ state: 'visible', timeout: 5_000 });
  }
}

/** Ensure the Context Panel (right side) is open and (optionally) all panes expanded. */
export async function ensureContextPanelOpen(page: Page, expandAll = false): Promise<void> {
  const cp = page.locator(CONTEXT_PANEL);
  if (!(await cp.isVisible().catch(() => false))) {
    await page.keyboard.press('F4');
    await cp.waitFor({ state: 'visible', timeout: 5_000 });
  }
  if (expandAll) {
    const btn = page.locator(CONTEXT_PANEL_EXPAND_ALL);
    if (await btn.isVisible().catch(() => false)) await btn.click();
  }
}

/** Click the "Home" icon in the Browse header. */
export async function clickHome(page: Page): Promise<void> {
  await page.locator(BROWSE_HEADER_HOME).click();
}

/**
 * Re-open Browse + ensure the parent group is expanded, then click the named item.
 * Convenient for tests that need to interact with multiple tree items in sequence,
 * because opening any item replaces the side panel with Toolbox.
 */
export async function openTreeItem(
  page: Page,
  parent: string,
  itemName: string,
  options?: { dblclick?: boolean },
): Promise<void> {
  await ensureBrowsePanelOpen(page);
  await expandTreeGroup(page, parent);
  const { treeItemByName } = await import('./selectors');
  const item = treeItemByName(page, itemName);
  await item.waitFor({ state: 'visible', timeout: 10_000 });
  await item.scrollIntoViewIfNeeded();
  if (options?.dblclick) await item.dblclick();
  else await item.click();
  await page.waitForTimeout(1000);
}

/** Click the "Collapse all" icon in the Browse header. */
export async function clickCollapseAll(page: Page): Promise<void> {
  await page.locator(BROWSE_HEADER_COLLAPSE_ALL).click();
  await page.waitForTimeout(400);
}

/** Expand a tree node (group) by name if not already expanded. Returns the label locator. */
export async function expandTreeGroup(page: Page, name: string): Promise<Locator> {
  const label = treeGroupByName(page, name);
  await label.waitFor({ state: 'visible', timeout: 10_000 });
  await label.scrollIntoViewIfNeeded();
  // Click the arrow if collapsed; click the label as a fallback.
  const node = label.locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]');
  const tri = node.locator(TREE_EXPAND_ARROW).first();
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false);
  if (!expanded) {
    if (await tri.isVisible().catch(() => false)) {
      await tri.click();
    } else {
      await label.click();
    }
    // Wait for the arrow to flip to expanded; fall back to a short sleep if the wait misses.
    await tri.evaluate((el) => new Promise<void>((resolve) => {
      if (el.classList.contains('d4-tree-view-tri-expanded')) { resolve(); return; }
      const o = new MutationObserver(() => {
        if (el.classList.contains('d4-tree-view-tri-expanded')) { o.disconnect(); resolve(); }
      });
      o.observe(el, { attributes: true, attributeFilter: ['class'] });
      setTimeout(() => { o.disconnect(); resolve(); }, 800);
    })).catch(() => page.waitForTimeout(400));
  }
  return label;
}

/** Returns how many tree-view nodes are currently in the "expanded" state. */
export async function countExpandedNodes(page: Page): Promise<number> {
  return page.locator(TREE_EXPAND_ARROW_EXPANDED).count();
}

/** Returns the top-level Browse tree node names (depth 0 of the tree). */
export async function getTopLevelTreeNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'));
    return labels
      .map((el) => {
        const node = el.closest('.d4-tree-view-node');
        if (!node) return null;
        // top-level node: no ancestor d4-tree-view-node above it (parents are tree-view-root-like)
        let p: HTMLElement | null = node.parentElement;
        while (p) {
          if (p.classList?.contains('d4-tree-view-node')) return null;
          p = p.parentElement;
        }
        const text = (el.textContent || '').trim();
        return text || null;
      })
      .filter((s): s is string => !!s);
  });
}
