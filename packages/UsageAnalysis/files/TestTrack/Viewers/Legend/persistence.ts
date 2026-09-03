import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

/** What a round-trip reads back, and what its settle stamps: every field the caller asserts. */
export interface StateProbe {
  viewerType?: string;
  props?: string[];
  column?: string;
  tags?: string[];
  legendItems?: boolean;
}

export interface ProbeResult {
  viewerCount: number;
  props: Record<string, any>;
  tags: Record<string, string | null>;
  legendItems: number | null;
}

export interface LayoutRoundTrip extends ProbeResult { layoutId: string; }

export interface ProjectRoundTrip extends ProbeResult {
  ok: boolean;
  phase: 'save' | 'reopen' | 'verified';
  error?: string;
  projectId: string | null;
}

async function installProbe(page: Page): Promise<void> {
  await v.installEventWaits(page);
  await page.evaluate(() => {
    const w = window as any;
    if (w.__legendProbe) return;
    w.__legendProbe = (spec: any) => {
      const tv = w.grok.shell.tv;
      const viewers = Array.from(tv?.viewers ?? []) as any[];
      const same = viewers.filter((x) => x.type === spec.viewerType);
      const props: Record<string, any> = {};
      for (const p of spec.props ?? []) props[p] = same[0] ? same[0].props[p] ?? null : null;
      const col = spec.column ? tv?.dataFrame?.col(spec.column) : null;
      const tags: Record<string, string | null> = {};
      for (const t of spec.tags ?? []) tags[t] = col ? col.tags[t] ?? null : null;
      const legendItems = spec.legendItems && same.length > 0
        ? Math.min(...same.map((x) => x.root.querySelectorAll('[name="legend"] .d4-legend-item').length))
        : null;
      return {viewerCount: viewers.length, props, tags, legendItems};
    };
    // a 30s server call is a failure, not a wait (a dev stall once held a save for 383s)
    w.__legendTo = (pr: Promise<any>, label: string) => Promise.race([pr,
      new Promise((_, rej) => setTimeout(() => rej(new Error(label + ' timed out after 30s')), 30000))]);
    // the settle has to see the legend exist before it can see it stop changing
    w.__legendSettle = async (spec: any, capMs: number) => {
      const stamp = () => JSON.stringify(w.__legendProbe(spec));
      if (spec.legendItems)
        await w.__poll(() => w.__legendProbe(spec).legendItems, (n: number | null) => (n ?? 0) > 0, capMs);
      await w.__settledFor(stamp, 250, 1500, 25);
    };
  });
}

/** Saves the current view as a layout, re-applies it and reads the probe back once the rebuilt view settles. */
export async function layoutRoundTrip(page: Page, namePrefix: string, spec: StateProbe): Promise<LayoutRoundTrip> {
  await installProbe(page);
  return page.evaluate(async ({prefix, spec}) => {
    const w = window as any;
    const grok = w.grok;
    const tv = grok.shell.tv;
    const layout = tv.saveLayout();
    layout.name = prefix + '_' + Date.now();
    const saved = await w.__legendTo(grok.dapi.layouts.save(layout), 'layouts.save');
    const found = await w.__findSaved(() => w.__legendTo(grok.dapi.layouts.find(saved.id), 'layouts.find'));
    const gen = w.__viewerGen();
    tv.loadLayout(found);
    await w.__rebuilt(gen, () => JSON.stringify(w.__legendProbe(spec)), 4500);
    await w.__legendSettle(spec, 3500);
    return {layoutId: String(saved.id), ...w.__legendProbe(spec)};
  }, {prefix: namePrefix, spec});
}

/** Saves the current view as a project, closes everything, reopens it and reads the probe back. */
export async function projectRoundTrip(page: Page, namePrefix: string, spec: StateProbe): Promise<ProjectRoundTrip> {
  await installProbe(page);
  return page.evaluate(async ({prefix, spec}) => {
    const w = window as any;
    const grok = w.grok;
    const DG = w.DG;
    const empty = {viewerCount: 0, props: {}, tags: {}, legendItems: null};
    let projectId: string | null = null;
    try {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const proj = DG.Project.create();
      proj.name = prefix + '_' + Date.now();
      const tableInfo = df.getTableInfo();
      const viewInfo = tv.getInfo();
      proj.addChild(tableInfo);
      proj.addChild(viewInfo);
      // a relation must point at an entity already persisted server-side, or projects.save
      // throws a project_relations FK violation — upload/save the table and view first
      await w.__legendTo(grok.dapi.tables.uploadDataFrame(df), 'tables.uploadDataFrame');
      await w.__legendTo(grok.dapi.tables.save(tableInfo), 'tables.save');
      await w.__legendTo(grok.dapi.views.save(viewInfo), 'views.save');
      const saved = await w.__legendTo(grok.dapi.projects.save(proj), 'projects.save');
      projectId = String(saved.id);
    } catch (e: any) {
      return {ok: false, phase: 'save', error: String(e).slice(0, 200), projectId, ...empty};
    }
    grok.shell.closeAll();
    await w.__poll(() => Array.from(grok.shell.tableViews).length, (c: number) => c === 0, 1200);
    try {
      const reopened = await w.__legendTo(grok.dapi.projects.find(projectId), 'projects.find');
      await w.__legendTo(reopened.open(), 'project.open');
    } catch (e: any) {
      return {ok: false, phase: 'reopen', error: String(e).slice(0, 200), projectId, ...empty};
    }
    // a reopened project lands the view, the dataFrame and the restored look in that
    // order, so readiness is the table and the settle is on what the caller reads
    await w.__tableReady(3500);
    await w.__legendSettle(spec, 3500);
    if (!grok.shell.tv)
      return {ok: false, phase: 'reopen', error: 'no tv after reopen', projectId, ...empty};
    return {ok: true, phase: 'verified', projectId, ...w.__legendProbe(spec)};
  }, {prefix: namePrefix, spec});
}

/** Deletes what a spec created, all deletes in flight together, then closes the shell. */
export async function deleteEntities(
  page: Page, ids: {layoutIds?: (string | null)[]; projectId?: string | null},
): Promise<void> {
  await page.evaluate(async ({layoutIds, projectId}) => {
    const w = window as any;
    const grok = w.grok;
    const jobs: Promise<void>[] = [];
    for (const id of layoutIds ?? [])
      if (id) jobs.push(grok.dapi.layouts.find(id).then((l: any) => l && grok.dapi.layouts.delete(l)).catch(() => {}));
    if (projectId)
      jobs.push(grok.dapi.projects.find(projectId).then((p: any) => p && grok.dapi.projects.delete(p)).catch(() => {}));
    await Promise.all(jobs);
    grok.shell.closeAll();
    await w.__poll(() => Array.from(grok.shell.tableViews).length, (c: number) => c === 0, 500);
  }, {layoutIds: ids.layoutIds ?? [], projectId: ids.projectId ?? null});
}
