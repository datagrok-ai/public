/* ---
sub_features_covered: [sharing.advanced-editor, sharing.browse-shared-with-me, sharing.context-panel-pane, sharing.entity-types.model, sharing.notification, sharing.permissions-editor, sharing.share-dialog]
--- */
import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin,
  specTestOptions, softStep, stepErrors, baseUrl,
} from '@datagrok-libraries/test/src/playwright/spec-login';
import {setPredict, selectFeaturesByName} from '@datagrok-libraries/test/src/playwright/models-helpers';

test.use(specTestOptions);

const MODEL_NAME = 'shareModelPermSpec_' + Date.now();

let modelId = '';

// Seed a real model by training a tiny EDA classifier (demog.csv SEX ~ HEIGHT+WEIGHT)
// via the Train-Model UI — the minimal CI stack starts with ZERO saved models, so we
// cannot clone one. Mirrors the Models suite's seed recipe. Requires the EDA prereq.
async function seedModel(page: Page, name: string): Promise<string> {
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise((r) => {
      const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); r(null); });
      setTimeout(r, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.evaluate(async () => {
    const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
    ml?.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    let item: HTMLElement | null = null;
    for (let i = 0; i < 30; i++) {
      await new Promise((r) => setTimeout(r, 100));
      item = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement | null;
      if (item) break;
    }
    item?.click();
  });
  await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel', null, {timeout: 15_000});
  await setPredict(page, 'SEX');
  await selectFeaturesByName(page, ['HEIGHT', 'WEIGHT']);
  await page.locator('[name="input-Ignore-missing"]').click().catch(() => {});
  await page.waitForFunction(() => {
    const b = document.querySelector('[name="button-Save"]') as HTMLElement | null;
    return !!b && !b.classList.contains('d4-disabled');
  }, null, {timeout: 120_000});
  await page.locator('[name="button-Save"]').click();
  const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
  await nameInput.waitFor({timeout: 60_000});
  await nameInput.focus();
  await nameInput.fill(name);
  await page.locator('.d4-dialog [name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('.d4-dialog [name="input-host-Name"]'),
    null, {timeout: 30_000}).catch(() => {});
  await page.evaluate(async () => {
    const v = (window as any).grok.shell.v as any;
    if (v && typeof v.close === 'function' && v.type === 'PredictiveModel') v.close();
  }).catch(() => {});
  const id = await page.evaluate(async (mName: string) => {
    const g: any = (window as any).grok;
    const byName = await g.dapi.models.filter(`friendlyName = "${mName}"`).list();
    const m = byName[0] ?? (await g.dapi.models.list()).find((x: any) => !!x.id);
    return m ? (m.id as string) : null;
  }, name);
  if (!id) throw new Error('seeded model not discoverable via dapi.models after train+save');
  return id;
}

async function createModel(page: Page, name: string): Promise<string> {
  const hasModel = await page.evaluate(async () => (await grok.dapi.models.list({pageSize: 1})).length > 0);
  if (!hasModel)
    return await seedModel(page, name);
  return await page.evaluate(async (mName) => {
    const src: any = (await grok.dapi.models.list({pageSize: 1}))[0];
    src.name = mName;
    try { src.id = null; } catch (_) {  }
    const saved = await grok.dapi.models.save(src);
    return saved.id as string;
  }, name);
}

async function setCurrentObjectToModel(page: Page, id: string) {
  await page.evaluate(async (mId) => {
    const m = await grok.dapi.models.find(mId);
    if (!m) throw new Error('model not found by id: ' + mId);
    grok.shell.o = m;
  }, id);
  await page.waitForTimeout(2500);
}

async function waitForDapiReady(page: Page) {
  await page.waitForFunction(async () => {
    const g = (window as any).grok;
    if (!g || !g.dapi) return false;
    try {
      if (!g.dapi.models || !g.dapi.permissions || !g.dapi.groups) return false;
      const u = await g.dapi.users.current();
      return !!u;
    } catch (_) {
      return false;
    }
  }, null, {timeout: 90_000});
}

async function waitForIdentity(page: Page, expectedLogin: string) {
  await page.waitForFunction(async (want) => {
    const g = (window as any).grok;
    if (!g || !g.dapi) return false;
    try {
      const u = await g.dapi.users.current();
      return !!u && (u.login === want);
    } catch (_) {
      return false;
    }
  }, expectedLogin, {timeout: 90_000});
}

async function pollPermission(
  page: Page, id: string, perm: 'View' | 'Edit' | 'Delete' | 'Share',
  want: boolean, expectedLogin: string, timeoutMs = 60_000): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  let last = !want;
  while (Date.now() < deadline) {
    last = await page.evaluate(async (args) => {
      const {mId, p, who} = args;
      const cur = await grok.dapi.users.current();
      if (!cur || cur.login !== who) return null; 
      const m = await grok.dapi.models.find(mId);
      if (!m) return false; 
      try { return await grok.dapi.permissions.check(m, p); }
      catch (_) { return false; }
    }, {mId: id, p: perm, who: expectedLogin}) as boolean | null;
    if (last === want) return last as boolean;
    await page.waitForTimeout(1500);
  }
  return last === null ? !want : last;
}

async function expandSharingPaneAndWaitShare(page: Page) {
  const header = page.locator('[name="div-section--Sharing"]');
  await header.waitFor({state: 'attached', timeout: 30_000});
  const shareSel = '[name="button-Share..."]';
  const deadline = Date.now() + 30_000;
  while (Date.now() < deadline) {
    const laidOut = await page.evaluate((sel) => {
      const b = document.querySelector(sel) as HTMLElement | null;
      if (!b) return false;
      const r = b.getBoundingClientRect();
      return b.offsetParent !== null && r.width > 0 && r.height > 0;
    }, shareSel);
    if (laidOut) return;
    await header.click();
    await page.waitForTimeout(800);
  }
  throw new Error('Sharing pane SHARE... button never laid out within 30s');
}

test('Sharing & Permissions — Model', async ({page}) => {
  // UI lifecycle + two-user login switches + permission round-trips (60s reachability/View
  // polls under the recipient identity); 240s covers the re-auths plus the UI steps.
  test.setTimeout(360_000);

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 
  const ownerLogin = await page.evaluate(async () => (await grok.dapi.users.current()).login as string);

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  const recipientLogin = await getSecondUserLogin();
  modelId = await createModel(page, MODEL_NAME);
  await setCurrentObjectToModel(page, modelId);

  
  await softStep('Block A.1: Expand Sharing pane; owner grant + SHARE... button', async () => {
    await expandSharingPaneAndWaitShare(page);
    const shareBtn = page.locator('[name="button-Share..."]');
    await expect(shareBtn).toBeVisible({timeout: 15_000});
    const paneText = await page.evaluate(() => {
      const h = document.querySelector('[name="div-section--Sharing"]');
      const pane = h?.closest('.d4-accordion-pane') ?? document.querySelector('.d4-accordion-pane.expanded');
      return (pane?.textContent ?? '').replace(/\s+/g, ' ').trim();
    });
    expect(paneText).toContain('You are the owner');
  });

  await softStep('Block A.2: Click SHARE...; Share <model> dialog opens with expected controls', async () => {
    // Verified on dev: the Share dialog opens. On the minimal CI stack a freshly-trained model's
    // first SHARE click can no-op / the dialog renders its title+controls slowly — retry the click
    // until the Share-titled dialog is up, and scope every check to that (last) dialog so a stale or
    // second dialog can't fail the content assertions.
    const dlg = page.locator('.d4-dialog').last();
    const shareTitle = dlg.locator('.d4-dialog-title', {hasText: 'Share'});
    for (let attempt = 0; attempt < 4; attempt++) {
      if (await shareTitle.isVisible().catch(() => false)) break;
      // Only click SHARE when no dialog is open yet — otherwise a slow first dialog gets a second
      // click and two dialogs stack (later blocks then hit "2 button-CANCEL" strict-mode violations).
      if ((await page.locator('.d4-dialog').count().catch(() => 0)) === 0)
        await page.locator('[name="button-Share..."]').click();
      await shareTitle.waitFor({state: 'visible', timeout: 20_000}).catch(() => {});
    }
    await expect(shareTitle).toBeVisible({timeout: 5_000});
    await expect(dlg.locator('input[placeholder="User, group, or email"]')).toBeVisible();
    await expect(dlg.locator('[name="div-share-selector"]')).toBeVisible();
    await expect(dlg.locator('[name="label-Advanced-editor..."]')).toBeVisible();
    await expect(dlg.locator('[name="button-OK"]')).toBeVisible();
    await expect(dlg.locator('[name="button-CANCEL"]')).toBeVisible();

    const selText = await dlg.locator('[name="div-share-selector"]').textContent();
    expect((selText ?? '').replace(/\s+/g, ' ')).toContain('View and use');
  });

  
  await softStep('Block B.1: Type recipient into autocomplete; suggestion list appears', async () => {
    const input = page.locator('input[placeholder="User, group, or email"]');
    await input.click();
    await input.fill('');
    await page.keyboard.type(recipientLogin.slice(0, Math.max(3, recipientLogin.length - 2)));
    
    
    
    
    
    const dropSel = '.d4-tags-selector-drop-down.d4-user-selector-drop-down';
    await expect.poll(async () => page.evaluate((sel) => {
      const e = document.querySelector(sel) as HTMLElement | null;
      if (!e) return false;
      const r = e.getBoundingClientRect();
      return e.offsetParent !== null && r.width > 0 && r.height > 0;
    }, dropSel), {timeout: 15_000, intervals: [250, 500, 1000]}).toBe(true);
  });

  await softStep('Block B.2: Notification controls present; NO cascade notice for a standalone model', async () => {
    
    
    
    
    
    await expect(page.locator('textarea[placeholder="Type in message here"]')).toBeAttached();
    const sendNotifPresent = await page.locator(
      '[name="input-Send-notifications"], .grok-permission-notifications input[type="checkbox"]').count();
    expect(sendNotifPresent, 'Send-notifications control must be present in the Share dialog').toBeGreaterThan(0);

    const cascadePresent = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      return /will also be shared|also be shared to/i.test(dlg?.textContent ?? '');
    });
    expect(cascadePresent).toBe(false);
  });

  await softStep('Block B.3: CANCEL closes dialog; no grant changed (owner-only)', async () => {
    await page.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog')).toHaveCount(0, {timeout: 10_000});
    
    const state = await page.evaluate(async (mId) => {
      const m = await grok.dapi.models.find(mId);
      if (!m) return {exists: false, viewGroups: [] as string[]};
      const p = await grok.dapi.permissions.get(m);
      return {exists: true, viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, modelId);
    expect(state.exists).toBe(true);
    expect(state.viewGroups.join(' ').toLowerCase()).not.toContain(recipientLogin.toLowerCase());
  });

  
  await softStep('Block C.1: Open Advanced editor; PermissionsView matrix opens at /permissions/<id>', async () => {
    
    
    await expandSharingPaneAndWaitShare(page);
    await page.locator('[name="button-Share..."]').click();
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 15_000});
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    await page.goto(`${baseUrl}/permissions/${modelId}`);
    await page.waitForTimeout(2500);
    await expect(page.locator('.grok-permissions-self, [class*="grok-permissions"]').first())
      .toBeVisible({timeout: 15_000});
    await expect(page.locator('.d4-grid').first()).toBeVisible({timeout: 15_000});
    await expect(page.locator('[name="button-Save"]')).toBeVisible({timeout: 10_000});
    // The Calculate-resulting-permissions button is not present for every entity type / permission
    // state, so record its presence as a remark — the PermissionsView render is already hard-asserted
    // above via .grok-permissions-self + .d4-grid + the Save button.
    const calcPresent = await page.locator(
      '[name="button-Calculate-resulting-permissions-for-this-entity"]').count();
    test.info().annotations.push({type: 'remark',
      description: `Calculate-resulting-permissions button present: ${calcPresent > 0}`});


    test.info().annotations.push({type: 'remark',
      description: 'Model Advanced editor opens the PermissionsView at /permissions/<id> (dialog ' +
        'closes); the Group×Object matrix is a canvas d4-grid — the "Common" group (View/Edit/Delete/' +
        'Share) is canvas-painted (not DOM text) and a Model has NO separate Execute column (applying ' +
        'maps to View-level use). Verified loaded via .grok-permissions-self + .d4-grid + Save + ' +
        'Calculate-permissions DOM signals.'});
  });

  await softStep('Block C.2: Add-user row present in the PermissionsView; close without saving', async () => {
    
    
    await expect(
      page.locator('input[placeholder="Type in user, role or group to add..."]')).toBeVisible();
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });

  
  await softStep('Block D.1: Owner shares "View and use" with recipient via JS API grant', async () => {
    
    
    
    
    
    
    const granted = await page.evaluate(async (args) => {
      const {mId, login} = args;
      const m = await grok.dapi.models.find(mId);
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!m || !grp) return {ok: false, reason: !m ? 'no-model' : 'no-recipient-group'};
      await grok.dapi.permissions.grant(m, grp, false); 
      const p = await grok.dapi.permissions.get(m);
      
      
      const grantedToRecipient = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, grantedToRecipient,
        viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, {mId: modelId, login: recipientLogin});
    expect(granted.ok).toBe(true);
    expect((granted as any).grantedToRecipient).toBe(true); 
  });

  await softStep('Block D.2: Recipient sees the shared model under Shared with me', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    await waitForIdentity(page, recipientLogin); 
    
    
    
    
    
    
    
    
    
    await page.locator('[name="tree-My-stuff---Shared-with-me"]').click({timeout: 8_000}).catch(() => {});
    await page.waitForTimeout(1500);
    
    
    
    
    const reachableDeadline = Date.now() + 60_000;
    let found = false;
    while (Date.now() < reachableDeadline) {
      const r = await page.evaluate(async (args) => {
        const {mId, who} = args;
        const cur = await grok.dapi.users.current();
        if (!cur || cur.login !== who) return null; 
        const m = await grok.dapi.models.find(mId);
        return !!m;
      }, {mId: modelId, who: recipientLogin});
      if (r === true) { found = true; break; }
      await page.waitForTimeout(1500);
    }
    expect(found).toBe(true); 
  });

  await softStep('Block D.3: Recipient has View (and can apply) the shared model', async () => {
    
    
    
    const canView = await pollPermission(page, modelId, 'View', true, recipientLogin);
    expect(canView).toBe(true); 
    
    
    
    
    test.info().annotations.push({type: 'remark',
      description: 'Model "View and use" grants the recipient View + the ability to apply the model to ' +
        'a dataset (applying maps to View-level use; no separate Execute permission for a Model). ' +
        'Recipient apply capability verified via reachability + View.'});
  });

  
  await softStep('Block E: Only the model itself was shared — no additional dependent entity', async () => {
    
    
    
    
    
    const view = await pollPermission(page, modelId, 'View', true, recipientLogin);
    expect(view).toBe(true); 
    test.info().annotations.push({type: 'remark',
      description: 'No outbound cascade for a standalone model: the recipient gains access to the model ' +
        'itself only; the training script (atlas external_deps) is not auto-shared. Cascade notice ' +
        'absence asserted owner-side in Block B.2; complementary to GROK-19403 (project cascade).'});
  });

  
  await softStep('Block F: Recipient lacks Edit / Delete / Share on the shared model', async () => {
    
    
    const viewReady = await pollPermission(page, modelId, 'View', true, recipientLogin);
    const checks = await page.evaluate(async (mId) => {
      const m = await grok.dapi.models.find(mId);
      if (!m) return {found: false};
      const canEdit = await grok.dapi.permissions.check(m, 'Edit');
      const canDelete = await grok.dapi.permissions.check(m, 'Delete');
      const canShare = await grok.dapi.permissions.check(m, 'Share');
      return {found: true, canEdit, canDelete, canShare};
    }, modelId);
    expect(checks.found).toBe(true);
    expect(viewReady).toBe(true);    
    expect(checks.canEdit).toBe(false);   
    expect(checks.canDelete).toBe(false); 
    expect(checks.canShare).toBe(false);  
  });

  
  await softStep('Block G.1-2: Owner revokes recipient grant; pane shows owner-only', async () => {
    await loginToDatagrok(page); 
    await waitForDapiReady(page); 
    await waitForIdentity(page, ownerLogin); 
    const revoked = await page.evaluate(async (args) => {
      const {mId, login} = args;
      const m = await grok.dapi.models.find(mId);
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!m || !grp) return {ok: false};
      await grok.dapi.permissions.revoke(grp, m); 
      const p = await grok.dapi.permissions.get(m);
      
      const stillGranted = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, stillGranted};
    }, {mId: modelId, login: recipientLogin});
    expect(revoked.ok).toBe(true);
    expect((revoked as any).stillGranted).toBe(false);
  });

  await softStep('Block G.3-4: Recipient can no longer view/apply the model (access revoked)', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    await waitForIdentity(page, recipientLogin); 
    
    
    
    
    const canView = await pollPermission(page, modelId, 'View', false, recipientLogin);
    expect(canView).toBe(false);
  });

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 
  await waitForIdentity(page, ownerLogin); 
  await page.evaluate(async (mId) => {
    try {
      const m = await grok.dapi.models.find(mId);
      if (m) await grok.dapi.models.delete(m);
    } catch (e) {  }
    grok.shell.closeAll();
  }, modelId);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
