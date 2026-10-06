import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const window: any;

export const ROW_COUNT = 5850;
export const GLYPH = {unchecked: 0xF0C8, checked: 0xF14A, indeterminate: 0xF146};

export interface HierNodeState {
  found: boolean;
  glyph: number;
  expanded: boolean;
  childCaptions: string[];
  childHeaders: string[];
}

export const glyphName = (g: number): string => 'U+' + g.toString(16).toUpperCase().padStart(4, '0');

export async function hierNode(page: Page, path: string[],
  action: 'probe' | 'read' | 'expand' | 'toggle'): Promise<HierNodeState> {
  return await page.evaluate(({path, action}: {path: string[], action: string}) => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
    if (!card)
      throw new Error('no hierarchical filter card is painted in the Filter Panel — its caption is the '
        + 'only one carrying "/", and no card in the panel carries one');
    const hostOf = (n: HTMLElement) =>
      n.parentElement?.querySelector(':scope > .d4-tree-view-group-host') as HTMLElement | null;
    const directChildren = (host: HTMLElement) => [...host.children].flatMap((c) =>
      c.matches('.d4-tree-view-node') ? [c as HTMLElement]
        : [...c.querySelectorAll(':scope > .d4-tree-view-node')] as HTMLElement[]);
    const captionOf = (n: HTMLElement) =>
      (n.querySelector('.d4-hierarchical-filter-caption-value')?.textContent ?? '').trim();
    const headerLabelOf = (n: HTMLElement) =>
      (n.querySelector('.d4-tree-view-item-label')?.textContent ?? n.textContent ?? '').trim();

    let node: HTMLElement | null = null;
    for (let level = 0; level < path.length; level++) {
      const value = path[level];
      let candidates: HTMLElement[];
      if (level === 0) {
        const rootHost = card.querySelector('.d4-tree-view-root > .d4-tree-view-group-host') as HTMLElement | null;
        if (!rootHost)
          throw new Error('the hierarchical filter card paints no tree root '
            + '(.d4-tree-view-root > .d4-tree-view-group-host), so its level-0 nodes cannot be addressed');
        candidates = directChildren(rootHost).filter((c) => captionOf(c) === value);
      }
      else {
        const host = hostOf(node!);
        candidates = host ? directChildren(host).filter((c) => captionOf(c) === value) : [];
      }
      if (candidates.length === 0) {
        if (action === 'probe')
          return {found: false, glyph: 0, expanded: false, childCaptions: [], childHeaders: []};
        throw new Error(`hierarchical tree node "${value}" not found (path ${path.join(' / ')}) — `
          + `it was looked for ${level === 0 ? 'inside the hierarchical card' : 'among the direct children of "' + path[level - 1] + '"'}`);
      }
      if (candidates.length > 1)
        throw new Error(`hierarchical tree node "${value}" is AMBIGUOUS (path ${path.join(' / ')}): `
          + `${candidates.length} nodes ${level === 0 ? 'inside the hierarchical card' : 'among the direct children of "' + path[level - 1] + '"'} `
          + 'carry that caption, so addressing it by caption would silently pick one of them');
      node = candidates[0];
    }
    const tri = node!.querySelector(':scope > .d4-tree-view-tri') as HTMLElement | null;
    if (action === 'expand') {
      if (!tri) throw new Error(`hierarchical tree node "${path.join(' / ')}" carries no expander`);
      if (!tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
    }
    if (action === 'toggle') {
      const cb = node!.querySelector('input.d4-hierarchical-filter-checkbox') as HTMLElement | null;
      if (!cb) throw new Error(`hierarchical tree node "${path.join(' / ')}" carries no checkbox`);
      cb.click();
    }
    const sub = node!.querySelector('.d4-hierarchical-filter-checkbox-substitute');
    const host = hostOf(node!);
    const children = host ? directChildren(host) : [];
    return {
      found: true,
      glyph: (sub?.textContent ?? '').codePointAt(0) ?? 0,
      expanded: !!tri?.classList.contains('d4-tree-view-tri-expanded'),
      childCaptions: children.map(captionOf).filter((c) => c !== ''),
      childHeaders: children.filter((c) => captionOf(c) === '').map(headerLabelOf),
    };
  }, {path, action});
}

export const trueCountOf = async (page: Page): Promise<number> =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

export async function hierCaption(page: Page): Promise<string> {
  return page.evaluate(() => {
    for (const c of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
      const cn = c.querySelector('.d4-filter-column-name');
      if (cn && cn.textContent!.includes('/')) return cn.textContent!.trim();
    }
    return '';
  });
}

export async function applyHierarchyState(page: Page, state: Record<string, any>): Promise<void> {
  await page.evaluate((s: Record<string, any>) => {
    const fg = grok.shell.tv.getFiltersGroup();
    for (const f of fg.filters) {
      if (f.filterType === 'hierarchical') {
        window.grok_GridFilterBase_ApplyState(f.dart || f, s);
        grok.shell.tv.dataFrame.rows.requestFilter();
        break;
      }
    }
  }, state);
}

export async function applyBoolState(page: Page, state: Record<string, any>): Promise<number> {
  return page.evaluate((s: Record<string, any>) => {
    const fg = grok.shell.tv.getFiltersGroup();
    for (const f of fg.filters) {
      if (f.filterType === 'bool-columns') {
        window.grok_GridFilterBase_ApplyState(f.dart ?? f, s);
        grok.shell.tv.dataFrame.rows.requestFilter();
        break;
      }
    }
    return grok.shell.tv.dataFrame.filter.trueCount;
  }, state);
}

// Step 5: the hierarchical card added through the panel menu and pointed at SEX / RACE.
export async function addHierarchicalCard(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Hierarchical');
  await expect.poll(async () => page.evaluate(() =>
    grok.shell.tv.getFiltersGroup().filters.filter((f: any) => f.filterType === 'hierarchical').length),
  {timeout: 15_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(1);
  await applyHierarchyState(page, {type: 'hierarchical', active: true, colNames: ['SEX', 'RACE'], allEnabled: true});
  await expect.poll(async () => (await hierNode(page, ['F'], 'probe')).found,
    {message: 'the SEX root node "F" never rendered — the hierarchy was not applied to the card',
      timeout: 10_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(true);
}

// Step 12: a fresh demog view with a second boolean column beside CONTROL.
export async function openDemogWithSexBool(page: Page, datasetPath: string): Promise<{type: string; hasControl: boolean}> {
  return page.evaluate(async (path: string) => {
    const w = window as any;
    grok.shell.closeAll();
    const df = await w.__readCsv(path);
    grok.shell.addTableView(df);
    await w.__tableReady(1500);
    await df.columns.addNewCalculated('SEX_bool', '${SEX} == "F"');
    const col = df.col('SEX_bool');
    return {type: col.type, hasControl: !!df.col('CONTROL')};
  }, datasetPath);
}

export function closeFilterPanelInPage(): void {
  const fv = document.querySelector('[name="viewer-Filters"]');
  let el: any = fv;
  while (el && !el.classList.contains('panel-base')) el = el.parentElement;
  const closeBtn = el?.querySelector('[name="Close"]') as HTMLElement | null;
  if (!closeBtn) throw new Error('the Filter Panel titlebar carries no [name="Close"] control — the '
    + 'close gesture the step is measured against never happened');
  closeBtn.click();
}
