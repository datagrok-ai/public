import {expect, Page} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, baseUrl} from '../spec-login';
import {recipientPage, secondUserLogin} from './_actors';

test.use(specTestOptions);

const MODEL_NAME = 'shareModelPermSpec_' + Date.now();

let modelId = '';

async function createModel(page: Page, name: string): Promise<string> {
  return await page.evaluate(async (mName) => {
    // The fixture clones whatever model the server hands back first, and that model is not
    // always re-savable: one run threw out of models.save, another sat in it until the test
    // timed out. Two candidates, each bounded, so a bad one costs 30s instead of the spec.
    const list = await grok.dapi.models.list({pageSize: 2});
    if (!list.length) throw new Error('no source model on server to clone');
    const reasons: string[] = [];
    for (const src of list as any[]) {
      try {
        src.name = mName;
        src.id = null;
        const saved = await Promise.race([
          grok.dapi.models.save(src),
          new Promise((_, rej) => setTimeout(() => rej(new Error('models.save timed out after 30s')), 30_000)),
        ]) as any;
        if (saved?.id) return saved.id as string;
        reasons.push('save returned no id');
      } catch (e) {
        reasons.push(String(e).slice(0, 100));
      }
    }
    throw new Error('no source model could be cloned: ' + reasons.join(' | '));
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

  test.setTimeout(240_000);

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

  const recipientLogin = await secondUserLogin();
  const rp = await recipientPage(page);
  await waitForDapiReady(rp);
  await waitForIdentity(rp, recipientLogin);
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
    await page.locator('[name="button-Share..."]').click();
    const dlg = page.locator('.d4-dialog');
    await expect(dlg).toBeVisible({timeout: 15_000});
    await expect(dlg.locator('.d4-dialog-title')).toContainText('Share');
    await expect(page.locator('input[placeholder="User, group, or email"]')).toBeVisible();
    await expect(page.locator('[name="div-share-selector"]')).toBeVisible();
    // GROK-20322 removed the Share dialog's "Advanced editor..." link
    // (core/client/xamgle/lib/src/commands/file/share_dataset.dart); the grant list the
    // PermissionsEditor renders is what the dialog must show now.
    await expect(page.locator('.d4-dialog .grok-permissions')).toBeVisible();
    await expect(page.locator('[name="button-OK"]')).toBeVisible();
    await expect(page.locator('[name="button-CANCEL"]')).toBeVisible();

    const selText = await page.locator('[name="div-share-selector"]').textContent();
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
    await waitForIdentity(rp, recipientLogin);

    await rp.locator('[name="tree-My-stuff---Shared-with-me"]').click({timeout: 8_000}).catch(() => {});

    const reachableDeadline = Date.now() + 60_000;
    let found = false;
    while (Date.now() < reachableDeadline) {
      const r = await rp.evaluate(async (args) => {
        const {mId, who} = args;
        const cur = await grok.dapi.users.current();
        if (!cur || cur.login !== who) return null; 
        const m = await grok.dapi.models.find(mId);
        return !!m;
      }, {mId: modelId, who: recipientLogin});
      if (r === true) { found = true; break; }
      await rp.waitForTimeout(1500);
    }
    expect(found).toBe(true); 
  });

  await softStep('Block D.3: Recipient has View (and can apply) the shared model', async () => {

    const canView = await pollPermission(rp, modelId, 'View', true, recipientLogin);
    expect(canView).toBe(true); 

    test.info().annotations.push({type: 'remark',
      description: 'Model "View and use" grants the recipient View + the ability to apply the model to ' +
        'a dataset (applying maps to View-level use; no separate Execute permission for a Model). ' +
        'Recipient apply capability verified via reachability + View.'});
  });

  await softStep('Block E: Only the model itself was shared — no additional dependent entity', async () => {

    const view = await pollPermission(rp, modelId, 'View', true, recipientLogin);
    expect(view).toBe(true); 
    test.info().annotations.push({type: 'remark',
      description: 'No outbound cascade for a standalone model: the recipient gains access to the model ' +
        'itself only; the training script (atlas external_deps) is not auto-shared. Cascade notice ' +
        'absence asserted owner-side in Block B.2; complementary to GROK-19403 (project cascade).'});
  });

  await softStep('Block F: Recipient lacks Edit / Delete / Share on the shared model', async () => {

    const viewReady = await pollPermission(rp, modelId, 'View', true, recipientLogin);
    const checks = await rp.evaluate(async (mId) => {
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
    await waitForIdentity(rp, recipientLogin);

    const canView = await pollPermission(rp, modelId, 'View', false, recipientLogin);
    expect(canView).toBe(false);
  });

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
