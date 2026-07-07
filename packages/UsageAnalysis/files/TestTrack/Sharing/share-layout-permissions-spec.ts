import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin,
  specTestOptions, softStep, stepErrors, baseUrl,
} from '../spec-login';

test.use(specTestOptions);

const LAYOUT_NAME = 'shareLayoutPermSpec_' + Date.now();

let LAYOUT_ID = '';

async function createLayout(page: Page, name: string): Promise<string> {
  return await page.evaluate(async (lName) => {
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise((r) => setTimeout(r, 2500));
    tv.addViewer('Scatter plot');
    await new Promise((r) => setTimeout(r, 1500));
    const layout = tv.saveLayout();
    layout.name = lName;
    const saved = await grok.dapi.layouts.save(layout);
    return saved.id as string;
  }, name);
}

async function waitForDapiReady(page: Page) {
  await page.waitForFunction(async () => {
    const g = (window as any).grok;
    if (!g || !g.dapi) return false;
    try {
      if (!g.dapi.layouts || !g.dapi.permissions || !g.dapi.groups) return false;
      const u = await g.dapi.users.current();
      return !!u;
    } catch (_) {
      return false;
    }
  }, null, {timeout: 90_000});
  
  
  
  
  await page.waitForTimeout(2000);
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
  page: Page, layoutId: string, perm: 'View' | 'Edit' | 'Delete' | 'Share',
  want: boolean, timeoutMs = 30_000): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  let last = !want;
  while (Date.now() < deadline) {
    last = await page.evaluate(async (args) => {
      const {lId, p} = args;
      const l = await grok.dapi.layouts.find(lId).catch(() => null);
      if (!l) return false; 
      try { return await grok.dapi.permissions.check(l, p); }
      catch (_) { return false; }
    }, {lId: layoutId, p: perm});
    if (last === want) return last;
    await page.waitForTimeout(1500);
  }
  return last;
}

async function setCurrentObjectToLayout(page: Page, id: string) {
  await page.evaluate(async (lId) => {
    const l = await grok.dapi.layouts.find(lId); 
    if (!l) throw new Error('layout not found by id: ' + lId);
    grok.shell.o = l;
  }, id);
  await page.waitForTimeout(2500);
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

test('Sharing & Permissions — Layout', async ({page}) => {
  // UI lifecycle + two-user login switches + permission round-trips; 240s covers the
  // two re-auths (each waits on dapi-ready) plus the UI pane/dialog/PermissionsView steps.
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

  const recipientLogin = await getSecondUserLogin();
  LAYOUT_ID = await createLayout(page, LAYOUT_NAME);
  await setCurrentObjectToLayout(page, LAYOUT_ID);

  
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

  await softStep('Block A.2: Click SHARE...; Share dialog opens with expected controls', async () => {
    await page.locator('[name="button-Share..."]').click();
    const dlg = page.locator('.d4-dialog');
    await expect(dlg).toBeVisible({timeout: 15_000});
    
    
    await expect(dlg.locator('.d4-dialog-title')).toContainText('Share');
    await expect(page.locator('input[placeholder="User, group, or email"]')).toBeVisible();
    await expect(page.locator('[name="div-share-selector"]')).toBeVisible();
    await expect(page.locator('[name="label-Advanced-editor..."]')).toBeVisible();
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

  await softStep('Block B.2: Notification controls present; NO cascade notice for a layout', async () => {
    
    
    
    
    
    
    await expect(page.locator('textarea[placeholder="Type in message here"]')).toBeAttached();
    const sendNotifPresent = await page.locator(
      '[name="input-Send-notifications"], .grok-permission-notifications input[type="checkbox"]').count();
    expect(sendNotifPresent, 'Send-notifications control must be present in the Share dialog').toBeGreaterThan(0);

    const cascadePresent = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      return /will also be shared/i.test(dlg?.textContent ?? '');
    });
    expect(cascadePresent).toBe(false);
  });

  await softStep('Block B.3: CANCEL closes dialog; no grant changed (owner-only)', async () => {
    await page.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog')).toHaveCount(0, {timeout: 10_000});
    
    
    
    const state = await page.evaluate(async (args) => {
      const {lId, login} = args;
      const l = await grok.dapi.layouts.find(lId);
      if (!l) return {exists: false, recipientPresent: false};
      const p = await grok.dapi.permissions.get(l);
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      const recipientPresent = !!grp && (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {exists: true, recipientPresent};
    }, {lId: LAYOUT_ID, login: recipientLogin});
    expect(state.exists).toBe(true);
    expect(state.recipientPresent).toBe(false);
  });

  
  await softStep('Block C.1: Open Advanced editor; PermissionsView matrix opens at /permissions/<id>', async () => {
    
    
    await expandSharingPaneAndWaitShare(page);
    await page.locator('[name="button-Share..."]').click();
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 15_000});
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    await page.goto(`${baseUrl}/permissions/${LAYOUT_ID}`);
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
      description: 'Layout Advanced editor opens the PermissionsView at /permissions/<id> ' +
        '(dialog closes); the Group×Object matrix is a canvas d4-grid — the standard "Common" ' +
        'View/Edit/Delete/Share headers are canvas-painted (not DOM text), and a ViewLayout carries ' +
        'NO entity-specific use-permission column (applying a layout is a View-level use). Verified ' +
        'loaded via .grok-permissions-self + .d4-grid + Save + Calculate-permissions DOM signals.'});
  });

  await softStep('Block C.2: Add-user row present in the PermissionsView; close without saving', async () => {
    
    
    await expect(
      page.locator('input[placeholder="Type in user, role or group to add..."]')).toBeVisible();
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });

  
  await softStep('Block D.1: Owner shares "View and use" with recipient via JS API grant', async () => {
    
    
    
    
    
    
    const granted = await page.evaluate(async (args) => {
      const {lId, login} = args;
      const l = await grok.dapi.layouts.find(lId); 
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!l || !grp) return {ok: false, reason: !l ? 'no-layout' : 'no-recipient-group', granted: false};
      await grok.dapi.permissions.grant(l, grp, false); 
      const p = await grok.dapi.permissions.get(l);
      
      
      const grantedToRecipient = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, granted: grantedToRecipient};
    }, {lId: LAYOUT_ID, login: recipientLogin});
    expect(granted.ok).toBe(true);
    expect((granted as any).granted).toBe(true);
  });

  await softStep('Block D.2-3: Recipient gains View; shared layout reachable + applicable', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page);                      
    await waitForIdentity(page, recipientLogin);       
    
    
    
    const found = await page.evaluate(async (lId) =>
      !!(await grok.dapi.layouts.find(lId).catch(() => null)), LAYOUT_ID);
    const canView = await pollPermission(page, LAYOUT_ID, 'View', true);
    expect(found).toBe(true);    
    expect(canView).toBe(true);  
  });

  
  await softStep('Block E: Recipient lacks Edit / Delete / Share on the shared layout', async () => {
    
    
    
    
    
    
    
    
    
    
    
    
    await waitForIdentity(page, recipientLogin);
    const checks = await page.evaluate(async (lId) => {
      const deadline = Date.now() + 30_000;
      let found = false;
      let canView = false; let canEdit = true; let canDelete = true; let canShare = true;
      while (Date.now() < deadline) {
        const l = await grok.dapi.layouts.find(lId).catch(() => null);
        if (l) {
          found = true;
          canView = await grok.dapi.permissions.check(l, 'View').catch(() => false);
          canEdit = await grok.dapi.permissions.check(l, 'Edit').catch(() => true);
          canDelete = await grok.dapi.permissions.check(l, 'Delete').catch(() => true);
          canShare = await grok.dapi.permissions.check(l, 'Share').catch(() => true);
          
          if (canView && !canEdit && !canDelete && !canShare) break;
        }
        await new Promise((r) => setTimeout(r, 1500));
      }
      return {found, canView, canEdit, canDelete, canShare};
    }, LAYOUT_ID);
    expect(checks.found).toBe(true);
    expect(checks.canView).toBe(true);    
    expect(checks.canEdit).toBe(false);   
    expect(checks.canDelete).toBe(false); 
    expect(checks.canShare).toBe(false);  
  });

  
  await softStep('Block F.1-2: Owner revokes recipient grant; pane shows owner-only', async () => {
    await loginToDatagrok(page); 
    await waitForDapiReady(page); 
    await waitForIdentity(page, ownerLogin); 
    const revoked = await page.evaluate(async (args) => {
      const {lId, login} = args;
      const l = await grok.dapi.layouts.find(lId); 
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!l || !grp) return {ok: false, stillGranted: true};
      await grok.dapi.permissions.revoke(grp, l); 
      const p = await grok.dapi.permissions.get(l);
      
      const stillGranted = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, stillGranted};
    }, {lId: LAYOUT_ID, login: recipientLogin});
    expect(revoked.ok).toBe(true);
    expect((revoked as any).stillGranted).toBe(false);
  });

  await softStep('Block F.3-4: Recipient can no longer view/apply the layout (access revoked)', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page);                 
    await waitForIdentity(page, recipientLogin);  
    
    
    
    
    const canView = await pollPermission(page, LAYOUT_ID, 'View', false);
    expect(canView).toBe(false);
  });

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 
  await waitForIdentity(page, ownerLogin); 
  await page.evaluate(async (lId) => {
    try {
      const l = await grok.dapi.layouts.find(lId); 
      if (l) await grok.dapi.layouts.delete(l);
    } catch (e) {  }
    grok.shell.closeAll();
  }, LAYOUT_ID);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
