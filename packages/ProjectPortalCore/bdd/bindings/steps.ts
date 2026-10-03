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

export const currentViewNamed = Then('the open view should be named {string}', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => grok.shell.v?.name ?? ''), {message: 'the open view'}).toBe(name);
}, {description: 'the current view by its name, e.g. a view another package opens'});

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
 * model_run rows, the saved calls and their uploaded tables, its snapshot files and its workspace folder. */
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
      // snapshot tables are files in the program's folder
      for (const sn of await grok.dapi.domains.table('portal.program_snapshot').query({filter: `program_id = "${r.id}"`})) {
        if (sn.table_file && sn.table_folder)
          await grok.dapi.spaces.id(sn.table_folder).files.delete(sn.table_file).catch(() => undefined);
      }
      // follows and pending approval requests would outlive the program: in My feed, in Waiting for me
      const follows = grok.dapi.domains.table('portal.follow');
      for (const f of await follows.query({filter: `program_id = "${r.id}"`, limit: 10000}))
        await follows.update(f.id, {active: false}).catch(() => undefined);
      const requests = grok.dapi.domains.table('portal.approval_request');
      for (const q of await requests.query({filter: `program_id = "${r.id}" and status = "pending"`, limit: 10000}))
        await requests.update(q.id, {status: 'withdrawn', result: 'withdrawn by the BDD cleanup'}).catch(() => undefined);
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
    if (!pr)
      return 0;
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

/** Removes [login]'s own role on [program] of [portfolio], if it has one. */
async function noRole(page: Page, login: string, program: string, portfolio: string): Promise<void> {
  await page.evaluate(async ([l, c, p]) => {
    await grok.functions.call('ProjectPortalCore:setProgramRole', {portfolio: p, program: c, principal: l, role: 'none'});
  }, [login, program, portfolio]);
}

export const noProgramRole = Given('{string} has no role on {string} in the {string} portfolio',
  async (page: Page, login: string, program: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    await noRole(page, login, program, portfolio);
    atFeatureEnd(page, async () => {
      await noRole(page, login, program, portfolio);
    });
  }, {tier: 'api', description: 'removes the user\'s program role through setProgramRole, now and at the end of the feature'});

export const scratchFolderFile = Given('a file {string} in the folder of the scratch project {string} in the {string} portfolio',
  async (page: Page, name: string, code: string, portfolio: string) => {
    silent(page);
    await page.evaluate(async ([n, c, p]) => {
      const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
      const pr = (await grok.dapi.domains.table('portal.program').query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`}))[0];
      if (!pr?.space_id)
        throw new Error(`${c} has no workspace folder`);
      await grok.dapi.spaces.id(pr.space_id).files.writeString(n, `# ${n}\n\nWritten by a BDD scenario; removed with the scratch project.\n`);
    }, [name, code, portfolio]);
  }, {tier: 'api', description: 'a small text file in the scratch program\'s workspace folder; the scratch cleanup removes the folder'});

/** Removes [model]'s predictions on [program] of [portfolio], with the scoring runs, calls and tables behind them. */
async function removePredictions(page: Page, model: string, program: string, portfolio: string): Promise<void> {
  await page.evaluate(async ([m, c, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const pr = (await grok.dapi.domains.table('portal.program').query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`}))[0];
    const predictions = grok.dapi.domains.table('portal.prediction');
    const runs = grok.dapi.domains.table('portal.model_run');
    const rows = await predictions.query({filter: `program_id = "${pr.id}" and model_code = "${m}"`});
    for (const r of rows)
      await predictions.delete(r.id).catch(() => undefined);
    for (const id of [...new Set(rows.map((r: any) => r.model_run_id).filter((x: any) => !!x))]) {
      const run = await runs.get(id).catch(() => null);
      const call = run?.run_id ? await grok.dapi.functions.calls.find(run.run_id).catch(() => null) : null;
      for (const t of call?.options?.['portalTables'] ?? [])
        await grok.dapi.tables.find(t).then((x: any) => grok.dapi.tables.delete(x)).catch(() => undefined);
      if (call)
        await grok.dapi.functions.calls.delete(call).catch(() => undefined);
      if (run)
        await runs.delete(id).catch(() => undefined);
    }
  }, [model, program, portfolio]);
}

export const notScored = Given('the example solubility model has not scored {string} in the {string} portfolio',
  async (page: Page, program: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    const model = 'ProjectPortalCore:SolubilityEstimate';
    await removePredictions(page, model, program, portfolio);
    atFeatureEnd(page, async () => {
      await removePredictions(page, model, program, portfolio);
    });
  }, {tier: 'api', description: 'removes the model\'s predictions on the program and the scoring runs behind them, now and at the end of the feature'});

/** Removes the folder [name] under [portfolio]'s workspace root, wherever parking put it: unparked
 * first through unparkFolder (which removes an empty _Parked), then deleted with its records. */
async function removeRootFolder(page: Page, name: string, portfolio: string): Promise<void> {
  await page.evaluate(async ([n, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const parked = grok.dapi.domains.table('portal.parked_folder');
    for (const r of await parked.query({filter: `portfolio_id = "${pf.id}" and name = "${n}"`})) {
      if (r.status === 'parked')
        await grok.functions.call('ProjectPortalCore:unparkFolder', {id: r.id}).catch(() => undefined);
      await parked.delete(r.id).catch(() => undefined);
    }
    const key = (s: string) => String(s ?? '').replace(/[^A-Za-z0-9]/g, '').toLowerCase();
    const root = (await grok.dapi.spaces.filter(`name = "${pf.workspace_root}"`).list({pageSize: 100}))
      .find((s: any) => key(s.friendlyName || s.name) === key(pf.workspace_root));
    if (!root)
      return;
    const kids = await grok.dapi.spaces.id(root.id).children.ofTypes('Project', false).list({pageSize: 1000});
    for (const k of kids.filter((s: any) => key(s.friendlyName || s.name) === key(n)))
      await grok.dapi.spaces.delete(k);
  }, [name, portfolio]);
}

export const rootFolder = Given('a folder {string} no program claims under the workspace root of the {string} portfolio',
  async (page: Page, name: string, portfolio: string) => {
    silent(page);
    await removeRootFolder(page, name, portfolio);
    atFeatureEnd(page, () => removeRootFolder(page, name, portfolio));
    await page.evaluate(async ([n, p]) => {
      const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
      const key = (s: string) => String(s ?? '').replace(/[^A-Za-z0-9]/g, '').toLowerCase();
      const root = (await grok.dapi.spaces.filter(`name = "${pf.workspace_root}"`).list({pageSize: 100}))
        .find((s: any) => key(s.friendlyName || s.name) === key(pf.workspace_root));
      if (!root)
        throw new Error(`the workspace of ${p} is not built`);
      await grok.dapi.spaces.id(root.id).addSubspace(n);
    }, [name, portfolio]);
  }, {tier: 'api', description: 'an empty folder at the top of the portfolio\'s workspace; unparked if parked and deleted at the end of the feature'});

/** Whether [program]'s folder is among the user's own favorites; pins it or not when [set] is given. */
async function pinnedOwn(page: Page, program: string, portfolio: string, set?: boolean): Promise<boolean> {
  return page.evaluate(async ([c, p, s]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const row = (await grok.dapi.domains.table('portal.program').query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`}))[0];
    if (!row?.space_id)
      throw new Error(`${c} has no workspace folder: build the ${p} workspace first`);
    const space = await grok.dapi.spaces.find(row.space_id);
    const on = (await grok.dapi.entities.getFavorites(DG.User.current().group)).some((e: any) => e.id === space.id);
    if (s !== null && s !== on)
      await (s ? DG.Favorites.add(space) : DG.Favorites.remove(space));
    return on;
  }, [program, portfolio, set ?? null] as [string, string, boolean | null]);
}

export const notPinned = Given('the user has not pinned {string} in the {string} portfolio',
  async (page: Page, program: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    const was = await pinnedOwn(page, program, portfolio, false);
    // the user's favorites end as they were found
    atFeatureEnd(page, async () => {
      await pinnedOwn(page, program, portfolio, was);
    });
  }, {tier: 'api', description: 'unpins the program\'s folder from the user\'s own favorites; at the end of the feature, puts back what was there'});

export const openHome = When('user opens Home', async (page: Page) => {
  await page.evaluate(() => {
    const home = Array.from(grok.shell.views as any[]).find((v) => v.name === 'Home');
    if (home)
      grok.shell.v = home;
    else
      grok.shell.route('/');
  });
  await page.locator('[data-u2-name="my-programs"]').waitFor({timeout: 60000});
}, {tier: 'ui', description: 'the platform\'s Home view with its widgets'});

export const modelPinned = Given('the user has pinned the {string} model',
  async (page: Page, model: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    const pin = (set: boolean | null) => page.evaluate(async ([m, s]) => {
      const [pkg, name] = m.split(':');
      const f = DG.Func.find({package: pkg, name})[0];
      const on = (await grok.dapi.entities.getFavorites(DG.User.current().group)).some((e: any) => e.id === f.id);
      if (s !== null && s !== on)
        await (s ? DG.Favorites.add(f) : DG.Favorites.remove(f));
      return on;
    }, [model, set] as [string, boolean | null]);
    const was = await pin(true);
    // the user's favorites end as they were found
    atFeatureEnd(page, async () => {
      await pin(was);
    });
  }, {tier: 'api', description: 'adds the model function to the user\'s own favorites; at the end of the feature, puts back what was there'});

export const selectInWorkspace = When('user selects the {string} model in the Workspace on Home',
  async (page: Page, model: string) => {
    const id = await page.evaluate((m) => {
      const [pkg, name] = m.split(':');
      return DG.Func.find({package: pkg, name})[0].id;
    }, model);
    // Home builds its Workspace once, so it is opened afresh to list the pin just made
    await page.goto('/', {waitUntil: 'domcontentloaded'});
    await page.locator('.d4-tab-header', {hasText: /^Workspace$/}).first().click({timeout: 120000});
    const row = page.locator(`.pp-workspace-item[data-entity-id="${id}"]`).first();
    await row.click({timeout: 60000});
  }, {tier: 'ui', description: 'Home, its Spotlight widget\'s Workspace tab, the model\'s row (PowerPack addresses the rows by entity id)'});

export const exampleDesignSets = Given('the example design sets are loaded', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => grok.functions.call('ProjectPortalCore:loadExampleDesignSets', {}));
}, {tier: 'api', description: 'the synthetic KIN-A design sets and the K1 hit list on EXS-KIN1; idempotent'});
