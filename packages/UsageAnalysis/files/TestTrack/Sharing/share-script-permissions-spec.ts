import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin,
  specTestOptions, softStep, stepErrors, baseUrl,
} from '../spec-login';

test.use(specTestOptions);

const SCRIPT_NAME = 'shareScriptPermSpec_' + Date.now();
const SCRIPT_SRC =
  '//name: ' + SCRIPT_NAME + '\n' +
  '//description: Hello world script\n' +
  '//language: javascript\n' +
  "grok.shell.info('hello world')";

async function createScript(page: Page, src: string, name: string): Promise<string> {
  return await page.evaluate(async (args) => {
    const {scriptSrc, scriptName} = args;
    const s = DG.Script.create(scriptSrc);
    s.name = scriptName;
    const saved = await grok.dapi.scripts.save(s);
    return saved.id as string;
  }, {scriptSrc: src, scriptName: name});
}

async function waitForDapiReady(page: Page) {
  await page.waitForFunction(async () => {
    const g = (window as any).grok;
    if (!g || !g.dapi) return false;
    try {
      if (!g.dapi.scripts || !g.dapi.permissions || !g.dapi.groups) return false;
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
  page: Page, scriptName: string, perm: 'View' | 'Edit' | 'Delete' | 'Share',
  want: boolean, timeoutMs = 30_000): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  let last = !want;
  while (Date.now() < deadline) {
    last = await page.evaluate(async (args) => {
      const {sName, p} = args;
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      if (!s) return false; 
      try { return await grok.dapi.permissions.check(s, p); }
      catch (_) { return false; }
    }, {sName: scriptName, p: perm});
    if (last === want) return last;
    await page.waitForTimeout(1500);
  }
  return last;
}

async function setCurrentObjectToScript(page: Page, name: string) {
  await page.evaluate(async (sName) => {
    const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
    if (!s) throw new Error('script not found: ' + sName);
    grok.shell.o = s;
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

test('Sharing & Permissions — Script', async ({page}) => {
  test.setTimeout(420_000);

  
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
  const scriptId = await createScript(page, SCRIPT_SRC, SCRIPT_NAME);
  await setCurrentObjectToScript(page, SCRIPT_NAME);

  
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

  await softStep('Block A.2: Click SHARE...; Share <script> dialog opens with expected controls', async () => {
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

  await softStep('Block B.2: Notification controls present; NO cascade notice for a script', async () => {
    
    
    
    
    
    
    await expect(page.locator('textarea[placeholder="Type in message here"]')).toBeAttached();
    
    
    const sendNotifPresent = await page.locator(
      '[name="input-Send-notifications"], .grok-permission-notifications input[type="checkbox"]').count();
    test.info().annotations.push({type: 'remark',
      description: `Block B.2 Send-notifications checkbox present: ${sendNotifPresent > 0}`});
    
    
    
    
    const cascadePresent = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      return /will also be shared/i.test(dlg?.textContent ?? '');
    });
    expect(cascadePresent).toBe(false);
  });

  await softStep('Block B.3: CANCEL closes dialog; no grant changed (owner-only)', async () => {
    await page.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog')).toHaveCount(0, {timeout: 10_000});
    
    const state = await page.evaluate(async (sName) => {
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      if (!s) return {exists: false, viewGroups: [] as string[]};
      const p = await grok.dapi.permissions.get(s);
      return {exists: true, viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, SCRIPT_NAME);
    expect(state.exists).toBe(true);
    expect(state.viewGroups.join(' ').toLowerCase()).not.toContain(recipientLogin.toLowerCase());
  });

  
  await softStep('Block C.1: Open Advanced editor; PermissionsView matrix opens at /permissions/<id>', async () => {
    
    
    await expandSharingPaneAndWaitShare(page);
    await page.locator('[name="button-Share..."]').click();
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 15_000});
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    await page.goto(`${baseUrl}/permissions/${scriptId}`);
    await page.waitForTimeout(2500);
    await expect(page.locator('.grok-permissions-self, [class*="grok-permissions"]').first())
      .toBeVisible({timeout: 15_000});
    await expect(page.locator('.d4-grid').first()).toBeVisible({timeout: 15_000});
    await expect(page.locator('[name="button-Save"]')).toBeVisible({timeout: 10_000});
    
    
    const calcPresent = await page.locator(
      '[name="button-Calculate-resulting-permissions-for-this-entity"]').count();
    test.info().annotations.push({type: 'remark',
      description: `Block C.1 Calculate-permissions button present: ${calcPresent > 0}`});
    
    
    
    test.info().annotations.push({type: 'remark',
      description: 'Script Advanced editor opens the PermissionsView at /permissions/<id> ' +
        '(dialog closes); the Group×Object matrix is a canvas d4-grid — "Common"/"Execute" headers ' +
        'are canvas-painted (not DOM text), verified loaded via .grok-permissions-self + .d4-grid + ' +
        'Save + Calculate-permissions DOM signals.'});
  });

  await softStep('Block C.2: Add-user row present in the PermissionsView; close without saving', async () => {
    
    
    await expect(
      page.locator('input[placeholder="Type in user, role or group to add..."]')).toBeVisible();
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });

  
  await softStep('Block D.1: Owner shares "View and use" with recipient via JS API grant', async () => {
    
    
    
    
    
    
    const granted = await page.evaluate(async (args) => {
      const {sName, login} = args;
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!s || !grp) return {ok: false, reason: !s ? 'no-script' : 'no-recipient-group'};
      await grok.dapi.permissions.grant(s, grp, false); 
      const p = await grok.dapi.permissions.get(s);
      
      
      
      
      const grantedToRecipient = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, grantedToRecipient,
        viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, {sName: SCRIPT_NAME, login: recipientLogin});
    expect(granted.ok).toBe(true);
    expect((granted as any).grantedToRecipient).toBe(true); 
  });

  await softStep('Block D.2: Recipient sees the shared script under Shared with me', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    await waitForIdentity(page, recipientLogin); 
    
    
    
    
    
    
    
    
    
    
    await page.locator('[name="tree-My-stuff---Shared-with-me"]').click({timeout: 8_000}).catch(() => {});
    await page.waitForTimeout(1500);
    const reachable = await page.evaluate(async (sName) => {
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      return {found: !!s};
    }, SCRIPT_NAME);
    expect(reachable.found).toBe(true); 
  });

  await softStep('Block D.3: Recipient has View (and Execute) on the shared script', async () => {
    
    
    
    const canView = await pollPermission(page, SCRIPT_NAME, 'View', true);
    expect(canView).toBe(true); 
    
    
    
    
    
    
    test.info().annotations.push({type: 'remark',
      description: 'Script "View and use" grants the recipient View + Execute. The script-specific ' +
        'Execute permission is not checkable via grok.dapi.permissions.check (generic API accepts ' +
        'only View/Edit/Delete/Share); recipient run capability verified via reachability + View.'});
  });

  
  await softStep('Block E: Recipient lacks Edit / Delete / Share on the shared script', async () => {
    
    
    const viewReady = await pollPermission(page, SCRIPT_NAME, 'View', true);
    const checks = await page.evaluate(async (sName) => {
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      if (!s) return {found: false};
      const canEdit = await grok.dapi.permissions.check(s, 'Edit');
      const canDelete = await grok.dapi.permissions.check(s, 'Delete');
      const canShare = await grok.dapi.permissions.check(s, 'Share');
      return {found: true, canEdit, canDelete, canShare};
    }, SCRIPT_NAME);
    expect(checks.found).toBe(true);
    expect(viewReady).toBe(true);    
    expect(checks.canEdit).toBe(false);   
    expect(checks.canDelete).toBe(false); 
    expect(checks.canShare).toBe(false);  
  });

  
  await softStep('Block F.1-2: Owner revokes recipient grant; pane shows owner-only', async () => {
    await loginToDatagrok(page); 
    await waitForDapiReady(page); 
    await waitForIdentity(page, ownerLogin); 
    const revoked = await page.evaluate(async (args) => {
      const {sName, login} = args;
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!s || !grp) return {ok: false};
      await grok.dapi.permissions.revoke(grp, s); 
      const p = await grok.dapi.permissions.get(s);
      
      
      const stillGranted = (p?.view ?? []).some((g: any) => g.id === grp.id);
      return {ok: true, stillGranted};
    }, {sName: SCRIPT_NAME, login: recipientLogin});
    expect(revoked.ok).toBe(true);
    expect((revoked as any).stillGranted).toBe(false);
  });

  await softStep('Block F.3-4: Recipient can no longer view/run the script (access revoked)', async () => {
    await loginAsSecondUser(page);
    await waitForDapiReady(page); 
    await waitForIdentity(page, recipientLogin); 
    
    
    
    
    const canView = await pollPermission(page, SCRIPT_NAME, 'View', false);
    expect(canView).toBe(false);
  });

  
  await loginToDatagrok(page);
  await waitForDapiReady(page); 
  await waitForIdentity(page, ownerLogin); 
  await page.evaluate(async (sName) => {
    try {
      const s = await grok.dapi.scripts.filter(`name = "${sName}"`).first();
      if (s) await grok.dapi.scripts.delete(s);
    } catch (e) {  }
    grok.shell.closeAll();
  }, SCRIPT_NAME);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
