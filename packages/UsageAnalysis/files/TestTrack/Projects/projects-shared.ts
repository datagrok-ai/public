import {Page, expect} from '@playwright/test';
import {baseUrl, getSecondUserLogin, resolveSecondUserToken} from '../spec-login';
import {PROVENANCE_PATTERNS, SourceClass} from '../helpers/openers';
import {ShareSecondUserResult, SavedAllTables, SavedProject} from '../helpers/projects';

/**
 * The second user, on a page of its own, booted once per worker.
 *
 * `helpers/projects.ts shareWithSecondUserAndVerify` re-injects a token into the shared owner
 * page and then calls `loginToDatagrok` to switch back — but that cannot switch back: its
 * `alreadyUp` early return (spec-login.ts) sees a booted shell and returns without injecting
 * the owner token, so the shared page stays authenticated as the recipient for every later
 * step and spec. Holding the recipient on a second context keeps both identities live and
 * costs one boot for the whole section instead of one per share step.
 */
let recipient: Page | null = null;

export async function recipientPage(owner: Page): Promise<Page> {
  if (recipient && !recipient.isClosed()) return recipient;
  const browser = owner.context().browser();
  if (!browser)
    throw new Error('recipientPage: the owner page has no browser to open a second context on');
  const token = await resolveSecondUserToken();
  const context = await browser.newContext({viewport: {width: 1920, height: 1080}});
  context.setDefaultTimeout(15_000);
  context.setDefaultNavigationTimeout(120_000);
  const p = await context.newPage();
  await p.goto(baseUrl + '/oauth/');
  await context.addCookies([{name: 'auth', value: token, domain: new URL(baseUrl).hostname, path: '/'}]);
  await p.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await p.goto(baseUrl);
  await p.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await p.locator('[name="Browse"]').waitFor({timeout: 60_000});
  recipient = p;
  return p;
}

/**
 * Grants the project to the second user and verifies the recipient can reach it — the same
 * contract as `helpers/projects.ts shareWithSecondUserAndVerify`, read on the recipient's own
 * page instead of re-authenticating the owner's.
 */
export async function shareWithSecondUserAndVerify(
  page: Page,
  project: {id?: string; name: string},
  options?: {full?: boolean},
): Promise<ShareSecondUserResult> {
  const secondLogin = await getSecondUserLogin();

  const grant = await page.evaluate(async (args: {id: string | null; name: string; wanted: string; full: boolean}) => {
    const grok = (window as any).grok;
    try {
      const target = await grok.dapi.users.filter(`login = "${args.wanted}"`).first();
      if (!target || !target.group || !target.group.id)
        return {shared: false, reason: `second user "${args.wanted}" / group not resolvable`, login: null};
      const p = args.id
        ? await grok.dapi.projects.find(args.id)
        : await grok.dapi.projects.filter(`name = "${args.name}"`).first();
      if (!p) return {shared: false, reason: 'project not found for share', login: null};
      await grok.dapi.permissions.grant(p, target.group, false);
      if (args.full) await grok.dapi.permissions.grant(p, target.group, true);

      let confirmed = true;
      try {
        const perms = await grok.dapi.permissions.get(p);
        const groups = [...(perms.view || []), ...(perms.edit || [])];
        confirmed = groups.some((g: any) => g && g.id === target.group.id);
      } catch (_) {  }
      return {
        shared: confirmed,
        reason: confirmed ? undefined : 'recipient group not in permissions.get',
        login: target.login,
      };
    } catch (e) {
      return {shared: false, reason: String(e).slice(0, 200), login: null};
    }
  }, {id: project.id ?? null, name: project.name, wanted: secondLogin, full: options?.full ?? false});

  const result: ShareSecondUserResult = {
    shared: grant.shared,
    recipientLogin: grant.login,
    recipientVisible: null,
    reason: grant.reason,
  };

  if (grant.shared && grant.login === secondLogin) {
    const rp = await recipientPage(page);
    let visible = false;
    await expect.poll(async () => {
      visible = await rp.evaluate(async (args: {id: string | null; name: string}) => {
        const grok = (window as any).grok;
        if (args.id) {
          try {
            if ((await grok.dapi.projects.find(args.id)) != null) return true;
          } catch (_) {  }
        }
        return (await grok.dapi.projects.filter(`name = "${args.name}"`).first()) != null;
      }, {id: project.id ?? null, name: project.name});
      return visible;
    }, {timeout: 30_000, intervals: [500, 1000, 2000, 5000]}).toBe(true).catch(() => {});
    result.recipientVisible = visible;
  }
  return result;
}

/**
 * In-page preamble shared by the save and reopen paths below: a bounded race so a dev stall
 * fails by name instead of holding a spec to its whole budget (one Case 5 upload sat until the
 * 420s test timeout and took the page down with it), and a state poll to replace the fixed
 * settles the shared helpers sleep through. The uploads stay serial — running them together
 * bought nothing measurable and cost one Dart-side failure.
 */
const PAGE_UTIL = `
  const bounded = (p, what, ms) => Promise.race([p,
    new Promise((_, rej) => setTimeout(() => rej(new Error(what + ' timed out after ' + ms + 'ms')), ms))]);
  const until = async (cond, cap) => {
    const t0 = Date.now();
    while (Date.now() - t0 < cap) {
      try { if (await cond()) return true; } catch (_) {  }
      await new Promise((r) => setTimeout(r, 100));
    }
    return false;
  };
`;

/** Saves every open table into one project, uploading the tables concurrently. */
export async function saveAllTablesWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedAllTables> {
  const run = () => page.evaluate(`(async () => {
    ${PAGE_UTIL}
    const n = ${JSON.stringify(projectName)};
    const tables = Array.from(grok.shell.tables);
    if (tables.length === 0) throw new Error('saveAllTablesWithProvenance: no tables in shell');
    const project = DG.Project.create();
    project.name = n;
    const infos = tables.map((df) => df.getTableInfo());
    for (const ti of infos) project.addChild(ti);
    for (let i = 0; i < tables.length; i++) {
      await bounded(grok.dapi.tables.uploadDataFrame(tables[i]), 'tables.uploadDataFrame', 60000);
      await bounded(grok.dapi.tables.save(infos[i]), 'tables.save', 45000);
    }
    const tv = grok.shell.tv;
    const layout = tv && tv.saveLayout ? tv.saveLayout() : null;
    if (layout) {
      project.addChild(layout);
      await bounded(grok.dapi.layouts.save(layout), 'layouts.save', 45000);
    }
    await bounded(grok.dapi.projects.save(project), 'projects.save', 45000);
    return {
      projectId: project.id,
      tableInfoIds: infos.map((t) => t.id),
      primaryTableInfoId: infos[0].id,
      layoutId: layout ? layout.id : null,
      resolvedName: project.name,
    };
  })()`) as Promise<SavedAllTables>;
  try { return await run(); }
  catch (e: any) {
    if (!/timed out after/.test(String(e?.message ?? e))) throw e;
    return await run();
  }
}

/** Saves the active TableView into one project. */
export async function saveProjectWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedProject> {
  const run = () => page.evaluate(`(async () => {
    ${PAGE_UTIL}
    const n = ${JSON.stringify(projectName)};
    const tv = grok.shell.tv;
    if (!tv || !tv.dataFrame) throw new Error('saveProjectWithProvenance: no active TableView');
    const df = tv.dataFrame;
    const project = DG.Project.create();
    project.name = n;
    const tableInfo = df.getTableInfo();
    project.addChild(tableInfo);
    const layout = tv.saveLayout ? tv.saveLayout() : null;
    if (layout) project.addChild(layout);
    await bounded(grok.dapi.tables.uploadDataFrame(df), 'tables.uploadDataFrame', 60000);
    await bounded(grok.dapi.tables.save(tableInfo), 'tables.save', 45000);
    if (layout) await bounded(grok.dapi.layouts.save(layout), 'layouts.save', 45000);
    await bounded(grok.dapi.projects.save(project), 'projects.save', 45000);
    return {
      projectId: project.id,
      tableInfoId: tableInfo.id,
      layoutId: layout ? layout.id : null,
      resolvedName: project.name,
    };
  })()`) as Promise<SavedProject>;
  try { return await run(); }
  catch (e: any) {
    if (!/timed out after/.test(String(e?.message ?? e))) throw e;
    return await run();
  }
}

export interface ReopenResult {
  reopenMs: number;
  tablesAfter: number;
  reopenedName: string;
  reopenedRowCount: number;
  reopenedScript: string;
}

/**
 * closeAll → reopen by id → read the re-materialized table, waiting on shell state rather
 * than the 1s + 2s the shared helper sleeps.
 */
export async function reopenAndAssertProvenance(
  page: Page,
  projectId: string,
  expectedScriptPattern?: RegExp,
): Promise<ReopenResult> {
  const result = await page.evaluate(`(async () => {
    ${PAGE_UTIL}
    const id = ${JSON.stringify(projectId)};
    grok.shell.closeAll();
    await until(() => Array.from(grok.shell.tableViews).length === 0, 3000);
    const t0 = performance.now();
    const proj = await grok.dapi.projects.find(id);
    if (!proj) throw new Error('Project not found by id: ' + id);
    await proj.open();
    const reopenMs = Math.round(performance.now() - t0);
    await until(() => {
      const d = grok.shell.tv && grok.shell.tv.dataFrame;
      return !!d && d.rowCount > 0;
    }, 20000);
    // a project can carry more than one table and they land a beat apart: hold until the
    // count stops growing, capped at the 2s settle this replaces
    let last = -1;
    for (let i = 0; i < 8; i++) {
      const n = Array.from(grok.shell.tables).length;
      if (n === last && n > 0) break;
      last = n;
      await new Promise((r) => setTimeout(r, 250));
    }
    const tv = grok.shell.tv;
    const df = tv && tv.dataFrame;
    if (!df) throw new Error('Reopen of ' + id + ': no df after open');
    return {
      reopenMs,
      tablesAfter: Array.from(grok.shell.tables).length,
      reopenedName: df.name,
      reopenedRowCount: df.rowCount,
      reopenedScript: (df.tags && df.tags.get ? df.tags.get('.script') : '') || '',
    };
  })()`) as ReopenResult;

  if (expectedScriptPattern && !expectedScriptPattern.test(result.reopenedScript)) {
    throw new Error(
      `reopenAndAssertProvenance: post-reopen .script = "${result.reopenedScript.slice(0, 200)}" ` +
      `does not match ${expectedScriptPattern}. The project saved successfully but provenance did ` +
      `not survive — likely cause: prior save path skipped tables.uploadDataFrame / tables.save / ` +
      `project.addChild(tableInfo).`);
  }
  return result;
}

export const PROVENANCE = PROVENANCE_PATTERNS;
export type {SourceClass};

/**
 * Queues a server-side delete that runs when the worker fixture drains
 * `window.__pendingDeletes`, not now.
 *
 * `Promise.all` subscribes to a thenable only when it runs, so the delete of a fixture that
 * has to outlive the test that created it is expressed as one — pushing a real Promise would
 * start the delete immediately and take the fixture out from under the next test.
 */
export async function deleteAtWorkerClose(
  page: Page,
  what: 'projects' | 'spaces' | 'queries' | 'scripts' | 'tables',
  id: string,
): Promise<void> {
  await page.evaluate(({ds, entityId}) => {
    const w = window as any;
    w.__pendingDeletes = w.__pendingDeletes ?? [];
    w.__pendingDeletes.push({
      then(res: any) {
        (async () => {
          try {
            const e = await w.grok.dapi[ds].find(entityId);
            if (e) await w.grok.dapi[ds].delete(e);
          } catch (_) {  }
          res();
        })();
      },
    });
  }, {ds: what, entityId: id}).catch(() => {});
}

