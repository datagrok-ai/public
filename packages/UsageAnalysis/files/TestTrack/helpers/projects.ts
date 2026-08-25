import {Page, expect} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser, getSecondUserLogin} from '../spec-login';

export interface ShareSecondUserResult {

  shared: boolean;

  recipientLogin: string | null;

  recipientVisible: boolean | null;

  reason?: string;
}

export async function shareWithSecondUserAndVerify(
  page: Page,
  project: {id?: string; name: string},
  options?: {full?: boolean},
): Promise<ShareSecondUserResult> {

  const secondLogin = await getSecondUserLogin();

  const grant = await page.evaluate(async (args: {id: string | null; name: string; wanted: string | null; full: boolean}) => {
    const grok = (window as any).grok;
    try {
      const target = args.wanted ? await grok.dapi.users.filter(`login = "${args.wanted}"`).first() : null;
      if (!target || !target.group || !target.group.id)
        return {shared: false, reason: `second user "${args.wanted}" / group not resolvable`, login: null, groupId: null};
      const p = args.id
        ? await grok.dapi.projects.find(args.id)
        : await grok.dapi.projects.filter(`name = "${args.name}"`).first();
      if (!p) return {shared: false, reason: 'project not found for share', login: null, groupId: null};
      await grok.dapi.permissions.grant(p, target.group, false);
      if (args.full) await grok.dapi.permissions.grant(p, target.group, true);

      let confirmed = true;
      try {
        const perms = await grok.dapi.permissions.get(p);
        const groups = [...(perms.view || []), ...(perms.edit || [])];
        confirmed = groups.some((g: any) => g && g.id === target.group.id);
      } catch (_) {  }
      return {shared: confirmed, reason: confirmed ? undefined : 'recipient group not in permissions.get', login: target.login, groupId: target.group.id};
    } catch (e) {
      return {shared: false, reason: String(e).slice(0, 200), login: null, groupId: null};
    }
  }, {id: project.id ?? null, name: project.name, wanted: secondLogin, full: options?.full ?? false});

  const result: ShareSecondUserResult = {
    shared: grant.shared,
    recipientLogin: grant.login,
    recipientVisible: null,
    reason: grant.reason,
  };

  if (grant.shared && grant.login && grant.login === secondLogin) {
    await loginAsSecondUser(page);
    try {
      let visible = false;
      await expect.poll(async () => {

        visible = await page.evaluate(async (args: {id: string | null; name: string}) => {
          const grok = (window as any).grok;
          if (args.id) {
            try {
              if ((await grok.dapi.projects.find(args.id)) != null) return true;
            } catch (_) {  }
          }
          return (await grok.dapi.projects.filter(`name = "${args.name}"`).first()) != null;
        }, {id: project.id ?? null, name: project.name});
        return visible;
      }, {timeout: 30_000, intervals: [1000, 2000, 5000]}).toBe(true).catch(() => {});
      result.recipientVisible = visible;
    } finally {
      await loginToDatagrok(page);
    }
  }
  return result;
}

export async function saveProject(
  page: Page,
  name: string,
  options?: {dataSyncMode?: 'on' | 'off'; cancelAutoShare?: boolean},
): Promise<void> {
  await page.evaluate(async ({n, sync}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const project = DG.Project.create();
    project.name = n;
    if (sync !== undefined)
      project.description = `[saveProject helper hint: dataSyncMode=${sync}]`;
    await grok.dapi.projects.save(project);
  }, {n: name, sync: options?.dataSyncMode});

  void options?.cancelAutoShare;
}

export async function waitForSaved(
  page: Page,
  name: string,
  options?: {timeout?: number},
): Promise<void> {
  const timeout = options?.timeout ?? 60000;
  await expect.poll(
    async () => page.evaluate(async (n) => {
      const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
      return p != null;
    }, name),
    {timeout, intervals: [500, 1000, 2000, 5000]},
  ).toBe(true);
}

export async function saveAndReopen(
  page: Page,
  name: string,
  options?: {dataSyncMode?: 'on' | 'off'},
): Promise<{tablesAfterReopen: number; viewersAfterReopen: number}> {
  await saveProject(page, name, {dataSyncMode: options?.dataSyncMode});
  await waitForSaved(page, name);

  await page.evaluate(async (n) => {
    (window as any).grok.shell.closeAll();
    const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
    await p.open();
  }, name);

  await page.waitForTimeout(2000);

  const tablesAfterReopen: number = await page.evaluate(() => {
    try { return Number((window as any).grok?.shell?.tables?.length) || 0; }
    catch (_) { return 0; }
  });
  const viewersAfterReopen: number = await page.evaluate(() => {
    try { return Number((window as any).grok?.shell?.tv?.viewers?.length) || 0; }
    catch (_) { return 0; }
  });
  return {tablesAfterReopen, viewersAfterReopen};
}

export async function openProjectFromDashboards(
  page: Page,
  name: string,
): Promise<void> {

  const slug = name.toLowerCase();
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});

  await tile.dblclick({force: true});

  await expect.poll(async () => page.evaluate((expected) => {
    const grok = (window as any).grok;
    return grok.shell.project != null && grok.shell.project.name === expected;
  }, name), {timeout: 30000, intervals: [500, 1000, 2000, 3000]}).toBe(true);
}

export async function buildVariantsComposite(
  page: Page,
  sourceProjectName: string,
  variants: Array<'link' | 'clone' | 'pvc'> = ['link', 'clone', 'pvc'],
): Promise<{original: string; link?: string; clone?: string; pvc?: string}> {
  const result: {original: string; link?: string; clone?: string; pvc?: string} =
    {original: sourceProjectName};

  for (const variant of variants) {
    const variantName = `${sourceProjectName}-${variant}`;
    await page.evaluate(async ({src, name, variant}) => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      const source = await grok.dapi.projects.filter(`name = "${src}"`).first();
      if (!source)
        throw new Error(`buildVariantsComposite: source "${src}" not found`);

      const p = DG.Project.create();
      p.name = name;
      p.description = `c1a buildVariantsComposite ${variant} variant of "${src}" (JS API approximation; saveCopy UI flow deferred to C1b)`;
      await grok.dapi.projects.save(p);
    }, {src: sourceProjectName, name: variantName, variant});
    result[variant] = variantName;
  }

  return result;
}

export async function collectChainProducedProjects(
  page: Page,
  namePrefix: string,
): Promise<Array<{name: string; id: string}>> {
  return await page.evaluate(async (prefix) => {
    const grok = (window as any).grok;

    const list = await grok.dapi.projects.filter(`name like "${prefix}%"`).list();
    return list.map((p: any) => ({name: p.name, id: p.id}));
  }, namePrefix);
}

export async function rename(
  page: Page,
  projectId: string,
  newName: string,
): Promise<{ok: boolean; newName: string | null; verifiedVia: 'find-by-id'}> {
  return await page.evaluate(async ({id, name}) => {
    const grok = (window as any).grok;
    const p = await grok.dapi.projects.find(id);
    if (!p)
      return {ok: false, newName: null, verifiedVia: 'find-by-id' as const};
    p.name = name;
    await grok.dapi.projects.save(p);

    const verify = await grok.dapi.projects.find(id);
    return {
      ok: verify != null && verify.name === name,
      newName: verify ? verify.name : null,
      verifiedVia: 'find-by-id' as const,
    };
  }, {id: projectId, name: newName});
}

export async function deleteProjectViaContextMenu(
  page: Page,
  name: string,
  mode: 'right-click' | 'context-panel-dropdown' = 'right-click',
): Promise<void> {
  const slug = name.toLowerCase();
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});

  if (mode === 'right-click') {

    await page.evaluate((s) => {
      const el = document.querySelector(`[name="div-${s}"]`) as HTMLElement | null;
      if (!el) throw new Error(`tile [name="div-${s}"] not found`);
      const r = el.getBoundingClientRect();
      el.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + r.width / 2,
        clientY: r.top + r.height / 2,
      }));
    }, slug);
  } else {

    await tile.click({force: true});
    await page.locator('.grok-prop-panel').waitFor({timeout: 15000});
    await page.locator('.grok-prop-panel [name="icon-context-arrow-down"]')
      .waitFor({state: 'attached', timeout: 15000});
    await page.evaluate(() => {
      const el = document.querySelector('.grok-prop-panel [name="icon-context-arrow-down"]') as HTMLElement | null;
      if (!el) throw new Error('chevron not in DOM');
      el.click();
    });
  }

  const deleteItem = page.locator('[name="div-Delete-Project"]');
  await deleteItem.waitFor({timeout: 10000});
  await deleteItem.click({force: true});

  const confirmDialog = page.locator('.d4-dialog')
    .filter({has: page.locator('[name="button-DELETE"]')});
  await confirmDialog.waitFor({timeout: 15000});
  await confirmDialog.locator('[name="button-DELETE"]').click();
  await expect(confirmDialog).toBeHidden({timeout: 15000});
}

export async function shareProjectViaContextMenu(
  page: Page,
  name: string,
  options: {
    recipient: string;
    accessLevel?: 'View and use' | 'Full access';
    sendNotifications?: boolean;
    message?: string;
  },
): Promise<void> {

  const candidates = [name, name.toLowerCase()];
  let slug = name;
  for (const cand of candidates) {
    const t = page.locator(`[name="div-${cand}"]`);
    if (await t.isVisible({timeout: 5_000}).catch(() => false)) {
      slug = cand;
      break;
    }
  }
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});

  await page.evaluate((s) => {
    const el = document.querySelector(`[name="div-${s}"]`) as HTMLElement | null;
    if (!el) throw new Error(`tile [name="div-${s}"] not found`);
    const r = el.getBoundingClientRect();
    el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: r.left + r.width / 2,
      clientY: r.top + r.height / 2,
    }));
  }, slug);

  const shareItem = page.locator('[name="div-Share..."]');
  await shareItem.waitFor({timeout: 10000});
  await shareItem.click({force: true});

  const dialog = page.locator('.d4-dialog').filter({hasText: `Share ${name}`});
  await dialog.waitFor({timeout: 15000});

  const recipientInput = dialog.locator('input[placeholder="User, group, or email"]');
  await recipientInput.focus();
  await page.keyboard.type(options.recipient);

  await page.waitForTimeout(1500);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(500);

  if (options.accessLevel && options.accessLevel !== 'View and use') {
    const accessSelect = dialog.locator('[name="div-share-selector"] select.ui-input-editor').first();
    if (await accessSelect.isVisible({timeout: 3000}).catch(() => false))
      await accessSelect.selectOption(options.accessLevel);
  }

  if (options.sendNotifications !== undefined) {
    const checkbox = dialog.locator('[name="input-Send-notifications"]');
    const current = await checkbox.isChecked().catch(() => null);
    if (current !== null && current !== options.sendNotifications)
      await dialog.locator('[name="input-host-Send-notifications"]').click();
  }

  if (options.message) {
    const messageArea = dialog.locator('textarea[placeholder="Type in message here"]');
    if (await messageArea.isVisible({timeout: 3000}).catch(() => false)) {
      await messageArea.focus();
      await page.keyboard.type(options.message);
    }
  }

  await dialog.locator('[name="button-OK"]').click();
  await expect(dialog).toBeHidden({timeout: 15000});
}

export async function saveCopy(
  page: Page,
  options: {
    sourceName: string;
    mode: 'original' | 'copy' | 'personal-view-customizations';
    name?: string;
    dataSyncMode?: 'on' | 'off';
    perTableLinkOrClone?: 'link' | 'clone';

    expectIdChange?: boolean;
  },
): Promise<void> {
  void options.sourceName; 

  const saveBtn = page.locator('button[name="button-Save"]:visible').first();
  await saveBtn.waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500); 
  try {
    await saveBtn.click({timeout: 10_000});
  } catch (_) {
    await page.evaluate(() => {
      const candidates = Array.from(document.querySelectorAll('button[name="button-Save"]'));
      const visible = candidates.find((b) => (b as HTMLElement).offsetParent !== null);
      if (visible) (visible as HTMLElement).click();
    });
  }

  const dialog = page.locator('.d4-dialog').filter({hasText: 'Save project'});
  await dialog.waitFor({timeout: 15000});

  const modeLabels: Record<typeof options.mode, string> = {
    'original': 'Save original project',
    'copy': 'Save a copy',
    'personal-view-customizations': 'Save personal view customizations',
  };
  await page.evaluate((targetText) => {
    const dialogEl = document.querySelector('.d4-dialog');
    if (!dialogEl) throw new Error('Save Project dialog not in DOM');
    const labels = Array.from(dialogEl.querySelectorAll('label'));
    const target = labels.find((l) => (l.textContent || '').trim() === targetText);
    if (!target) throw new Error(`mode label "${targetText}" not found`);

    const inputId = (target as HTMLLabelElement).getAttribute('for');
    if (inputId) {
      const input = document.getElementById(inputId) as HTMLInputElement | null;
      if (input) {
        input.checked = true;
        input.dispatchEvent(new Event('change', {bubbles: true}));
        input.dispatchEvent(new Event('input', {bubbles: true}));
      }
    }
    (target as HTMLLabelElement).click();
  }, modeLabels[options.mode]);
  await page.waitForTimeout(800); 

  if (options.name !== undefined) {

    const nameInput = dialog.locator('input#name, input[type="text"].ui-input-editor').first();
    await nameInput.waitFor({state: 'visible', timeout: 15000});
    await page.evaluate((newName) => {
      const input = (document.querySelector('.d4-dialog input#name')
        ?? document.querySelector('.d4-dialog input[type="text"].ui-input-editor')) as HTMLInputElement | null;
      if (!input) throw new Error('saveCopy: Name input not found at fill-time (tried input#name and input[type="text"].ui-input-editor)');
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, newName);
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    }, options.name);
  } else if (options.mode !== 'copy') {
    throw new Error(`saveCopy: name is required for mode '${options.mode}' (only mode='copy' supports omitted-name → "Copy of <sourceName>" default)`);
  }

  if (options.dataSyncMode !== undefined) {
    const desiredOn = options.dataSyncMode === 'on';
    const dsCheckbox = dialog.locator('[name="input-host-Data-sync"] input[type="checkbox"]');
    const current = await dsCheckbox.isChecked().catch(() => null);
    if (current !== null && current !== desiredOn)
      await dialog.locator('[name="input-host-Data-sync"]').click();
  }

  if (options.mode === 'copy' && options.perTableLinkOrClone) {
    const choiceContainers = dialog.locator('.grok-project-move-entity-debug.ui-input-choice');
    const count = await choiceContainers.count();
    for (let i = 0; i < count; i++) {
      const container = choiceContainers.nth(i);

      const choiceLabel = container.locator('label').filter({
        hasText: new RegExp(`^${options.perTableLinkOrClone === 'link' ? 'Link' : 'Clone'}$`, 'i'),
      });
      if (await choiceLabel.isVisible({timeout: 2000}).catch(() => false))
        await choiceLabel.click({force: true});
    }
  }

  const preOkSourceId = await page.evaluate(() => {
    const grok = (window as any).grok;
    return grok.shell.project?.id ?? null;
  });

  await dialog.locator('[name="button-OK"]').click();

  const shareDialog = page.locator('.d4-dialog').filter({hasText: /^Share /});
  if (await shareDialog.isVisible({timeout: 10000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }

  const expectIdChange = options.expectIdChange ?? (options.mode !== 'original');

  const verified = await page.evaluate(async ({preId, expectChange, typedName}: {preId: string | null; expectChange: boolean; typedName: string | null}) => {
    const grok = (window as any).grok;
    const deadline = Date.now() + 20000;
    while (Date.now() < deadline) {

      const shellProj = grok.shell.project;
      if (shellProj?.id) {
        const okShell = !expectChange || shellProj.id !== preId;
        if (okShell) {
          const found = await grok.dapi.projects.find(shellProj.id);
          if (found?.name)
            return {id: shellProj.id, name: found.name, ok: true, via: 'shell'};
        }
      }

      if (expectChange && typedName) {
        try {
          const list = await grok.dapi.projects.filter(`name like "${typedName}"`).list();
          const candidate = list.find((p: any) => p.id !== preId);
          if (candidate)
            return {id: candidate.id, name: candidate.name, ok: true, via: 'friendly-like'};
        } catch (_) {  }
      }
      await new Promise((r) => setTimeout(r, 1000));
    }
    return {
      ok: false,
      reason: expectChange
        ? `no copy found via shell update or friendlyName like "${typedName}" within 20s (preId=${preId})`
        : 'timeout waiting for shell.project to settle with non-null find result',
    };
  }, {preId: preOkSourceId, expectChange: expectIdChange, typedName: options.name ?? null});

  if (!verified.ok)
    throw new Error(`saveCopy verification failed: ${JSON.stringify(verified)}`);
}

export interface SavedProject {

  projectId: string;

  tableInfoId: string;

  layoutId: string | null;

  resolvedName: string;
}

export async function saveProjectWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedProject> {
  return await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const tv = grok.shell.tv;
    if (!tv?.dataFrame)
      throw new Error('saveProjectWithProvenance: no active TableView');
    const df = tv.dataFrame;
    const project = DG.Project.create();
    project.name = n;
    const tableInfo = df.getTableInfo();
    project.addChild(tableInfo);
    const layout = tv.saveLayout?.() ?? null;
    if (layout) project.addChild(layout);
    if (layout) await grok.dapi.layouts.save(layout);
    await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.projects.save(project);
    return {
      projectId: project.id,
      tableInfoId: tableInfo.id,
      layoutId: layout?.id ?? null,
      resolvedName: project.name,
    };
  }, projectName);
}

export interface SavedAllTables {

  projectId: string;

  tableInfoIds: string[];

  primaryTableInfoId: string;

  layoutId: string | null;

  resolvedName: string;
}

export async function saveAllTablesWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedAllTables> {
  return await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const tables = Array.from(grok.shell.tables);
    if (tables.length === 0)
      throw new Error('saveAllTablesWithProvenance: no tables in shell');
    const project = DG.Project.create();
    project.name = n;
    const tableInfoIds: string[] = [];
    for (const df of tables as any[]) {
      const ti = df.getTableInfo();
      project.addChild(ti);
      await grok.dapi.tables.uploadDataFrame(df);
      await grok.dapi.tables.save(ti);
      tableInfoIds.push(ti.id);
    }
    const tv = grok.shell.tv;
    const layout = tv?.saveLayout?.() ?? null;
    if (layout) {
      project.addChild(layout);
      await grok.dapi.layouts.save(layout);
    }
    await grok.dapi.projects.save(project);
    return {
      projectId: project.id,
      tableInfoIds,
      primaryTableInfoId: tableInfoIds[0],
      layoutId: layout?.id ?? null,
      resolvedName: project.name,
    };
  }, projectName);
}

export interface ReopenResult {

  reopenMs: number;

  tablesAfter: number;

  reopenedName: string;

  reopenedRowCount: number;

  reopenedScript: string;
}

export async function reopenAndAssertProvenance(
  page: Page,
  projectId: string,
  expectedScriptPattern?: RegExp,
): Promise<ReopenResult> {
  const result = await page.evaluate(async (id) => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 1000));
    const t0 = performance.now();
    const proj = await grok.dapi.projects.find(id);
    if (!proj) throw new Error(`Project not found by id: ${id}`);
    await proj.open();
    const reopenMs = Math.round(performance.now() - t0);

    await new Promise((r) => setTimeout(r, 2000));
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!df) throw new Error(`Reopen of ${id}: no df after open + 2s settle`);
    const tablesAfter = (() => {
      try { return Number(grok.shell.tables?.length) || 0; }
      catch (_) { return 0; }
    })();
    return {
      reopenMs,
      tablesAfter,
      reopenedName: df.name,
      reopenedRowCount: df.rowCount,
      reopenedScript: df.tags?.get?.('.script') ?? '',
    };
  }, projectId);

  if (expectedScriptPattern && !expectedScriptPattern.test(result.reopenedScript)) {
    throw new Error(
      `reopenAndAssertProvenance: post-reopen .script = ` +
      `"${result.reopenedScript.slice(0, 200)}" does not match ` +
      `${expectedScriptPattern}. The project saved successfully but ` +
      `provenance did not survive — likely cause: prior save path skipped ` +
      `tables.uploadDataFrame / tables.save / project.addChild(tableInfo).`,
    );
  }
  return result;
}

/**
 * Saves the active TableView as a project entirely through the JS API — no ribbon click,
 * no Save/Share dialogs, no polling. `tv.getInfo()` captures the view state (viewer configs
 * included) the same way the ribbon Save does internally; `await`ing `dapi.projects.save()`
 * IS the completion signal, so there is nothing left to poll for afterwards.
 */
export async function saveProjectViaApi(
  page: Page,
  name: string,
): Promise<{projectId: string; resolvedName: string}> {
  const result = await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const errAt = (step: string, e: any) =>
      `at ${step}: ${e?.message ?? e?.toString?.() ?? e} | stack=${String(e?.stack ?? '').slice(0, 300)}`;
    const tv = grok.shell.tv;
    if (!tv?.dataFrame)
      return {error: 'no active TableView'};
    const df = tv.dataFrame;
    let project: any, tableInfo: any, viewInfo: any;
    try { project = DG.Project.create(); project.name = n; }
    catch (e) { return {error: errAt('Project.create', e)}; }
    try { tableInfo = df.getTableInfo(); }
    catch (e) { return {error: errAt('df.getTableInfo', e)}; }
    try { viewInfo = tv.getInfo(); }
    catch (e) { return {error: errAt('tv.getInfo', e)}; }
    try { project.addChild(tableInfo); }
    catch (e) { return {error: errAt('addChild(tableInfo)', e)}; }
    try { project.addChild(viewInfo); }
    catch (e) { return {error: errAt('addChild(viewInfo)', e)}; }
    try { await grok.dapi.tables.uploadDataFrame(df); }
    catch (e) { return {error: errAt('tables.uploadDataFrame', e)}; }
    try { await grok.dapi.tables.save(tableInfo); }
    catch (e) { return {error: errAt('tables.save', e)}; }
    // a project relation must point at an entity that already exists server-side —
    // the ViewInfo needs its own save, the same way a linked ViewLayout does
    try { await grok.dapi.views.save(viewInfo); }
    catch (e) { return {error: errAt('views.save', e)}; }
    try { await grok.dapi.projects.save(project); }
    catch (e) { return {error: errAt('projects.save', e)}; }
    return {projectId: String(project.id), resolvedName: String(project.name)};
  }, name);
  if ('error' in result)
    throw new Error(`saveProjectViaApi failed ${result.error}`);
  return result;
}

export async function saveProjectViaUI(
  page: Page,
  name: string,
): Promise<{projectId: string; resolvedName: string}> {
  await page.locator('[name="button-Save"]:visible').first().click();
  await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
  await page.locator('.d4-dialog input[type="text"]').first().fill(name);
  await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});
  await page.waitForTimeout(3000);

  const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
  if (await cancel.count() > 0)
    await cancel.click({force: true});
  await page.waitForTimeout(800);

  const found = await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    for (let a = 0; a < 25; a++) {
      try {
        const p = await grok.dapi.projects.filter(`name = "${n}"`).first();
        if (p) return {id: String(p.id), name: String(p.name)};
      } catch (_) {  }
      try {

        const recent = await grok.dapi.projects.list({pageSize: 50});
        const hit = (recent || []).find((p: any) =>
          p && (p.friendlyName === n || p.name === n));
        if (hit) return {id: String(hit.id), name: String(hit.name)};
      } catch (_) {  }
      await new Promise((r) => setTimeout(r, 1200));
    }
    return null;
  }, name);
  if (!found)
    throw new Error(`saveProjectViaUI: project "${name}" not visible server-side after ribbon save (polled filter()-by-name + list() scan for ~30s)`);
  return {projectId: found.id, resolvedName: found.name};
}

export async function deleteProjectWithCleanup(
  page: Page,
  ids: {projectId?: string; tableInfoId?: string; scriptId?: string},
): Promise<void> {
  await page.evaluate(async (i) => {
    const grok = (window as any).grok;
    if (i.projectId) {
      try {
        const p = await grok.dapi.projects.find(i.projectId);
        if (p) await grok.dapi.projects.delete(p);
      } catch (_) {  }
    }
    if (i.tableInfoId) {
      try {
        const ti = await grok.dapi.tables.find(i.tableInfoId);
        if (ti) await grok.dapi.tables.delete(ti);
      } catch (_) {  }
    }
    if (i.scriptId) {
      try {
        const s = await grok.dapi.scripts.find(i.scriptId);
        if (s) await grok.dapi.scripts.delete(s);
      } catch (_) {  }
    }
  }, ids).catch(() => {});
}
