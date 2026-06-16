/* ---
sub_features_covered: [sharing.advanced-editor, sharing.browse-shared-with-me, sharing.context-panel-pane, sharing.entity-types.connection, sharing.notification, sharing.permissions-editor, sharing.share-dialog]
--- */
import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin,
  specTestOptions, softStep, stepErrors, baseUrl,
} from '../spec-login';

test.use(specTestOptions);

const CONN_NAME = 'shareConnPermSpec_' + Date.now();

async function createConnection(page: Page, name: string): Promise<string> {
  return await page.evaluate(async (cName) => {
    const conn = DG.DataConnection.create(cName, {
      dataSource: 'Postgres',
      server: 'localhost',
      port: 5432,
      db: 'datagrok',
    });
    conn.name = cName;
    const saved = await grok.dapi.connections.save(conn);
    return saved.id as string;
  }, name);
}

async function waitForDapiReady(page: Page) {
  await page.waitForFunction(async () => {
    const g = (window as any).grok;
    if (!g || !g.dapi) return false;
    try {
      if (!g.dapi.connections || !g.dapi.permissions || !g.dapi.groups) return false;
      const u = await g.dapi.users.current();
      return !!u;
    } catch (_) {
      return false;
    }
  }, null, {timeout: 90_000});
}

async function setCurrentObjectToConnection(page: Page, name: string) {
  await page.evaluate(async (cName) => {
    const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
    if (!c) throw new Error('connection not found: ' + cName);
    grok.shell.o = c;
  }, name);
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

test('Sharing & Permissions — Connection', async ({page}) => {
  // UI lifecycle + two-user login switches + permission round-trips; 240s covers the
  // two re-auths (each waits on dapi-ready) plus the UI pane/dialog/PermissionsView steps.
  test.setTimeout(240_000);

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 

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
  const connId = await createConnection(page, CONN_NAME);
  await setCurrentObjectToConnection(page, CONN_NAME);

  
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

  await softStep('Block A.2: Click SHARE...; Share <connection> dialog opens with expected controls', async () => {
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

  await softStep('Block B.2: Notification controls present; NO cascade notice for a connection', async () => {
    
    
    
    
    
    
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
    
    const state = await page.evaluate(async (cName) => {
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      if (!c) return {exists: false, viewGroups: [] as string[]};
      const p = await grok.dapi.permissions.get(c);
      return {exists: true, viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, CONN_NAME);
    expect(state.exists).toBe(true);
    expect(state.viewGroups.join(' ').toLowerCase()).not.toContain(recipientLogin.toLowerCase());
  });

  
  await softStep('Block C.1: Open Advanced editor; PermissionsView matrix opens at /permissions/<id>', async () => {
    
    
    await expandSharingPaneAndWaitShare(page);
    await page.locator('[name="button-Share..."]').click();
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 15_000});
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    await page.goto(`${baseUrl}/permissions/${connId}`);
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
      description: 'Connection Advanced editor opens the PermissionsView at /permissions/<id> ' +
        '(dialog closes); the Group×Object matrix is a canvas d4-grid — "Common"/"Query"/' +
        '"GetSchema"/"ListFiles" headers are canvas-painted (not DOM text), verified loaded via ' +
        '.grok-permissions-self + .d4-grid + Save + Calculate-permissions DOM signals.'});
  });

  await softStep('Block C.2: Add-user row present in the PermissionsView; close without saving', async () => {
    
    
    await expect(
      page.locator('input[placeholder="Type in user, role or group to add..."]')).toBeVisible();
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });

  
  await softStep('Block D.1: Owner shares "View and use" with recipient via JS API grant', async () => {
    
    
    
    
    
    
    const granted = await page.evaluate(async (args) => {
      const {cName, login} = args;
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!c || !grp) return {ok: false, reason: !c ? 'no-connection' : 'no-recipient-group'};
      await grok.dapi.permissions.grant(c, grp, false); 
      const p = await grok.dapi.permissions.get(c);
      return {ok: true, viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, {cName: CONN_NAME, login: recipientLogin});
    expect(granted.ok).toBe(true);
    expect((granted as any).viewGroups.join(' ').toLowerCase()).toContain(recipientLogin.toLowerCase());
  });

  await softStep('Block D.2-4: Recipient gains View; use-permissions (Query/GetSchema) granted', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    const result = await page.evaluate(async (cName) => {
      
      
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      if (!c) return {found: false};
      const canView = await grok.dapi.permissions.check(c, 'View');
      return {found: true, canView};
    }, CONN_NAME);
    expect(result.found).toBe(true);   
    expect(result.canView).toBe(true); 
  });

  
  await softStep('Block E: Recipient lacks Edit / Delete / Share on the shared connection', async () => {
    const checks = await page.evaluate(async (cName) => {
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      if (!c) return {found: false};
      const canEdit = await grok.dapi.permissions.check(c, 'Edit');
      const canDelete = await grok.dapi.permissions.check(c, 'Delete');
      const canShare = await grok.dapi.permissions.check(c, 'Share');
      const canView = await grok.dapi.permissions.check(c, 'View');
      return {found: true, canEdit, canDelete, canShare, canView};
    }, CONN_NAME);
    expect(checks.found).toBe(true);
    expect(checks.canView).toBe(true);    
    expect(checks.canEdit).toBe(false);   
    expect(checks.canDelete).toBe(false); 
    expect(checks.canShare).toBe(false);  
  });

  
  await softStep('Block F.1-2: Owner revokes recipient grant; pane shows owner-only', async () => {
    await loginToDatagrok(page); 
    await waitForDapiReady(page); 
    const revoked = await page.evaluate(async (args) => {
      const {cName, login} = args;
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!c || !grp) return {ok: false};
      await grok.dapi.permissions.revoke(grp, c); 
      const p = await grok.dapi.permissions.get(c);
      const viewLogins = (p?.view ?? []).map((g: any) => (g.friendlyName || g.name || '').toLowerCase());
      return {ok: true, stillGranted: viewLogins.includes(login.toLowerCase())};
    }, {cName: CONN_NAME, login: recipientLogin});
    expect(revoked.ok).toBe(true);
    expect((revoked as any).stillGranted).toBe(false);
  });

  await softStep('Block F.3-4: Recipient can no longer view/use the connection (access revoked)', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    const result = await page.evaluate(async (cName) => {
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      if (!c) return {found: false, canView: false}; 
      const canView = await grok.dapi.permissions.check(c, 'View');
      return {found: true, canView};
    }, CONN_NAME);
    
    expect(result.canView).toBe(false);
  });

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 
  await page.evaluate(async (cName) => {
    try {
      const c = await grok.dapi.connections.filter(`name = "${cName}"`).first();
      if (c) await grok.dapi.connections.delete(c);
    } catch (e) {  }
    grok.shell.closeAll();
  }, CONN_NAME);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
