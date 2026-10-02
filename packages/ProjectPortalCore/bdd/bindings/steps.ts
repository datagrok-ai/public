/* The steps only this package can define: how to reach the app, a scratch program to take
   through a gate, and reading a board's stage columns in order. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

export const APP_PATH = '/apps/ProjectPortalCore';

export const openApp = Given('user opens the Project Portal app', async (page: Page) => {
  await page.goto(APP_PATH, {waitUntil: 'domcontentloaded'});
  await page.locator('[data-u2-name="portfolio-list"] .pp-portfolio').first().waitFor({timeout: 60000});
}, {tier: 'ui', enters: 'Project Portal app',
  description: 'needs this package published and ProjectPortalCore:loadPortfolios + syncAllPortfolios run on the stand'});

/** Removes every live program with [code] in [portfolio]; the sweep at the start and the cleanup at the end. */
async function removeScratch(page: Page, code: string, portfolio: string): Promise<number> {
  return page.evaluate(async ([c, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const programs = grok.dapi.domains.table('portal.program');
    const rows = await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`});
    for (const r of rows)
      await programs.delete(r.id);
    return (await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`})).length;
  }, [code, portfolio]);
}

export const scratchProject = Given('a scratch project {string} in the {string} portfolio',
  async (page: Page, code: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    await removeScratch(page, code, portfolio);
    atFeatureEnd(page, async () => {
      const left = await removeScratch(page, code, portfolio);
      if (left !== 0)
        throw new Error(`${left} scratch program(s) ${code} survived the cleanup`);
    });
    await page.evaluate(async ([c, p]) => {
      const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
      await grok.dapi.domains.table('portal.program').insert({portfolio_id: pf.id, code: c, name: c,
        stage_model_id: pf.stage_model_id, source: 'manual'});
    }, [code, portfolio]);
    // the board was drawn before the program existed: redraw it
    await page.evaluate(async () => {
      const view = grok.shell.v;
      await view?.portal?.load?.();
    });
  }, {tier: 'api', description: 'a manual program, swept first and deleted at the end of the feature'});

export const boardStages = Then('the board should show the stages {string}', async (page: Page, codes: string) => {
  const expected = codes.split(',').map((s) => s.trim());
  await expect.poll(() => page.locator('[data-u2-name="stage-board"] .pp-column-title').allTextContents()
    .then((titles) => titles.map((t) => t.split(' · ')[0]).filter((c) => c !== 'Unstaged')),
  {message: `stage columns in order`}).toEqual(expected);
}, {description: 'the stage codes of the selected portfolio\'s model, left to right, unstaged column excluded'});

export const sidePanelsHidden = Given('the side panels are hidden', async (page: Page) => {
  silent(page);   // framing for a guide, not something the viewer does
  const before = await page.evaluate(() => {
    const w = grok.shell.windows;
    const was = {browse: w.showBrowse, context: w.showContextPanel, help: w.showHelp};
    w.showBrowse = false;
    w.showContextPanel = false;
    w.showHelp = false;
    // a dev stand announces its debug-mode packages; not part of the story
    document.querySelectorAll('.d4-balloon').forEach((b) => b.remove());
    return was;
  });
  atFeatureEnd(page, async () => {
    await page.evaluate((was) => {
      grok.shell.windows.showBrowse = was.browse;
      grok.shell.windows.showContextPanel = was.context;
      grok.shell.windows.showHelp = was.help;
    }, before);
  });
}, {tier: 'ui', description: 'hides the browse tree, context panel and help for the rest of the feature'});

export const addressEndsWith = Then('the address should end with {string}', async (page: Page, tail: string) => {
  await expect.poll(() => new URL(page.url()).pathname.endsWith(tail),
    {message: `the address ends with ${tail}`}).toBe(true);
}, {description: 'the path of the browser address, so a board can be bookmarked and shared'});

export const currentSpace = Then('the open view should be the space {string}', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => grok.shell.v?.project?.friendlyName ?? grok.shell.v?.name ?? ''),
    {message: 'the open view'}).toBe(name);
}, {description: 'the Spaces view the Browse tree opens for a space, by the space\'s display name'});

/** The open portfolio table's viewers by title: those over the programs and those over its linked tables. */
export const tableViewers = Then('the portfolio table should show the viewers {string}', async (page: Page, titles: string) => {
  const expected = titles.split(',').map((s) => s.trim());
  await expect.poll(() => page.evaluate(() => {
    const tv = grok.shell.tv;
    const viewers = [...Array.from(tv?.viewers ?? []), ...(tv?.extraViewers ?? [])] as any[];
    return viewers.filter((v) => v.root.isConnected && v.root.offsetHeight > 0).map((v) => v.props.title);
  }).then((shown: string[]) => expected.filter((t) => !shown.includes(t))),
  {message: 'viewers of the portfolio table not shown'}).toEqual([]);
}, {description: 'every titled viewer is on screen in the current portfolio table view'});

export const matrixFollows = Then('the readiness matrix should show the programs that pass the filter', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const tv = grok.shell.tv;
    const t = tv?.portfolioTables;
    if (!t)
      return 'not a portfolio table';
    const ids = (df: any) => {
      const out: string[] = [];
      for (let i = df.filter.findNext(-1, true); i !== -1; i = df.filter.findNext(i, true))
        out.push(df.get('id', i));
      return out.sort().join();
    };
    return ids(t.readiness) === ids(tv.dataFrame) ? 'same' : `${ids(t.readiness)} ≠ ${ids(tv.dataFrame)}`;
  }), {message: 'the readiness matrix against the programs'}).toBe('same');
}, {description: 'the readiness matrix is a linked table: its rows passing the filter are those of the programs'});

/** Removes every live program with [code] in [portfolio], with what saving runs to it left: the
 * model_run rows, the saved calls and their uploaded tables, and its workspace folder. */
async function removeScratchWithRuns(page: Page, code: string, portfolio: string): Promise<number> {
  return page.evaluate(async ([c, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const programs = grok.dapi.domains.table('portal.program');
    const runs = grok.dapi.domains.table('portal.model_run');
    const removeSpace = async (id: string): Promise<void> => {
      const children = await grok.dapi.spaces.id(id).children.ofTypes('Project', false).list({pageSize: 1000}).catch(() => []);
      for (const child of children)
        await removeSpace(child.id);
      const space = await grok.dapi.spaces.find(id).catch(() => null);
      if (space)
        await grok.dapi.spaces.delete(space);
    };
    for (const r of await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`})) {
      for (const run of await runs.query({filter: `program_id = "${r.id}"`})) {
        const call = run.run_id ? await grok.dapi.functions.calls.find(run.run_id).catch(() => null) : null;
        for (const id of (call?.options?.portalTables ?? []) as string[]) {
          const t = await grok.dapi.tables.find(id).catch(() => null);
          if (t)
            await grok.dapi.tables.delete(t);
        }
        if (call)
          await grok.dapi.functions.calls.delete(call);
        await runs.delete(run.id);
      }
      if (r.space_id)
        await removeSpace(r.space_id);
      await programs.delete(r.id);
    }
    return (await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`})).length;
  }, [code, portfolio]);
}

export const scratchProjectAtStage = Given('a scratch project {string} at stage {string} in the {string} portfolio',
  async (page: Page, code: string, stage: string, portfolio: string) => {
    silent(page);
    await removeScratchWithRuns(page, code, portfolio);
    atFeatureEnd(page, async () => {
      const left = await removeScratchWithRuns(page, code, portfolio);
      if (left !== 0)
        throw new Error(`${left} scratch program(s) ${code} survived the cleanup`);
    });
    await page.evaluate(async ([c, s, p]) => {
      const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
      const st = await grok.dapi.domains.table('portal.stage').getByKey({stage_model_id: pf.stage_model_id, code: s});
      await grok.dapi.domains.table('portal.program').insert({portfolio_id: pf.id, code: c, name: c,
        stage_model_id: pf.stage_model_id, stage_id: st.id, source: 'manual', modality: 'mAb'});
    }, [code, stage, portfolio]);
  }, {tier: 'api', description: 'a manual program at a stage; it, its saved runs, their calls and tables, and its workspace folder are swept first and deleted at the end of the feature'});

export const exampleRunSaved = Given('a run of the example model saved to the scratch project {string} in the {string} portfolio',
  async (page: Page, code: string, portfolio: string) => {
    silent(page);
    await page.evaluate(async ([c, p]) => {
      const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
      const call = await f.prepare({feedRate: 3, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
      await grok.functions.call('ProjectPortalCore:saveModelRun', {call, portfolio: p, program: c,
        unitOperation: 'EX-PROD', runType: 'EX-SIM', purpose: 'characterization'});
    }, [code, portfolio]);
  }, {tier: 'api', description: 'runs Fed-batch titer (example) and saves the run to the program with saveModelRun'});

export const scratchWorkspace = Given('the workspace of the scratch project {string} in the {string} portfolio is built',
  async (page: Page, code: string, portfolio: string) => {
    silent(page);
    await page.evaluate(async ([c, p]) => grok.functions.call('ProjectPortalCore:buildWorkspace',
      {code: p, dryRun: false, program: c}), [code, portfolio]);
  }, {tier: 'api', description: 'builds the folders of that one program; the scratch cleanup removes them'});

export const runExampleModel = When('user runs the example model', async (page: Page) => {
  await page.evaluate(async () => {
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    await f.prepare({feedRate: 3, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
  });
}, {tier: 'api', description: 'runs Fed-batch titer (example) on its default inputs, as its run view would'});

export const titledRunSaved = Given('a run of the example model at feed rate {int} titled {string} saved to the scratch project {string} in the {string} portfolio',
  async (page: Page, feedRate: number, title: string, code: string, portfolio: string) => {
    silent(page);
    await page.evaluate(async ([rate, t, c, p]) => {
      const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
      const call = await f.prepare({feedRate: rate, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
      await grok.functions.call('ProjectPortalCore:saveModelRun', {call, portfolio: p, program: c,
        unitOperation: 'EX-PROD', runType: 'EX-SIM', purpose: 'characterization', title: t});
    }, [feedRate, title, code, portfolio] as const);
  }, {tier: 'api', description: 'runs Fed-batch titer (example) at that feed rate and saves it to the program under that title'});

/** The current table view's table: its row count and name. */
const openTable = (page: Page) => page.evaluate(() => {
  const df = grok.shell.tv?.dataFrame;
  return df ? {rows: df.rowCount as number, name: df.name as string} : null;
});

export const openTableStarting = Then('the open table should have {int} rows and a name starting with {string}',
  async (page: Page, rows: number, head: string) => {
    await expect.poll(() => openTable(page).then((t) => t && t.name.startsWith(head) ? t.rows : t?.name ?? 'no table view'),
      {message: `a table view of ${rows} rows named ${head}…`}).toBe(rows);
  }, {description: 'the table of the current table view'});

export const openTableEnding = Then('the open table should have {int} rows and a name ending with {string}',
  async (page: Page, rows: number, tail: string) => {
    await expect.poll(() => openTable(page).then((t) => t && t.name.endsWith(tail) ? t.rows : t?.name ?? 'no table view'),
      {message: `a table view of ${rows} rows named …${tail}`}).toBe(rows);
  }, {description: 'the table of the current table view'});

export const openTableViewers = Then('the open table should show the viewer types {string}', async (page: Page, types: string) => {
  const expected = types.split(',').map((s) => s.trim());
  await expect.poll(() => page.evaluate(() => Array.from(grok.shell.tv?.viewers ?? []).map((v: any) => v.type))
    .then((shown: string[]) => expected.filter((t) => !shown.includes(t))),
  {message: 'viewer types not in the current table view'}).toEqual([]);
}, {description: 'every viewer type is in the current table view, besides its grid'});

export const openTableCurrent = Then('the open table\'s current row should have {string} in {string}',
  async (page: Page, value: string, column: string) => {
    await expect.poll(() => page.evaluate((c) => {
      const df = grok.shell.tv?.dataFrame;
      return df && df.currentRowIdx >= 0 ? String(df.get(c, df.currentRowIdx)) : 'no current row';
    }, column), {message: `the current row's ${column}`}).toBe(value);
  }, {description: 'a cell of the current row of the current table view'});

export const portfolioSynced = Given('a sync of the {string} portfolio has run', async (page: Page, code: string) => {
  silent(page);
  await page.evaluate(async (c) => {
    await grok.functions.call('ProjectPortalCore:syncPortfolio', {code: c, dryRun: false});
  }, code);
}, {tier: 'api', description: 'syncs the portfolio by hand through syncPortfolio, which records the run'});

/** Turns off every follow the user has on [program] of [portfolio]: the program's and its work units'. */
async function unfollowAll(page: Page, program: string, portfolio: string): Promise<number> {
  return page.evaluate(async ([c, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const pr = (await grok.dapi.domains.table('portal.program').query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`}))[0];
    const follows = grok.dapi.domains.table('portal.follow');
    const mine = (await follows.query({filter: `program_id = "${pr.id}" and user_id = "${DG.User.current().id}"`}))
      .filter((f: any) => f.active !== false);
    for (const f of mine)
      await follows.update(f.id, {active: false});
    return mine.length;
  }, [program, portfolio]);
}

export const notFollowed = Given('the user does not follow {string} in the {string} portfolio',
  async (page: Page, program: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    await unfollowAll(page, program, portfolio);
    atFeatureEnd(page, async () => {
      await unfollowAll(page, program, portfolio);
    });
  }, {tier: 'api', description: 'turns off the user\'s follows of the program and its series, now and at the end of the feature'});

/** Forgets the user's mode and persona, so the portal is in Advanced with the persona from the portfolio. */
async function clearPreferences(page: Page): Promise<void> {
  await page.evaluate(async () => {
    for (const key of ['mode', 'persona']) {
      // deleting a key that is not there throws
      if (grok.userSettings.getValue('ProjectPortalCore', key) != null)
        grok.userSettings.delete('ProjectPortalCore', key);
    }
    await grok.userSettings.flush();
  });
}

export const preferencesCleared = Given('the user\'s portal mode and persona are not set', async (page: Page) => {
  silent(page);   // set-up a person does not do: kept out of guide videos
  await clearPreferences(page);
  atFeatureEnd(page, async () => {
    await clearPreferences(page);
  });
}, {tier: 'api', description: 'clears the remembered mode and persona from the user\'s settings, now and at the end of the feature'});
