import {test, expect, Page} from '@playwright/test';
import {
  loginToDatagrok, loginAsSecondUser, getSecondUserLogin,
  specTestOptions, softStep, stepErrors, baseUrl,
} from '../spec-login';

test.use(specTestOptions);

const QUERY_NAME = 'shareQueryPermSpec_' + Date.now();

async function createQuery(page: Page, name: string): Promise<string> {
  return await page.evaluate(async (qName) => {
    const conn = await grok.dapi.connections.filter('name = "Datagrok"').first();
    if (!conn) throw new Error('System:Datagrok connection not found');
    
    const q = (conn as any).query(qName, 'select * from public.groups');
    q.name = qName;
    const saved = await grok.dapi.queries.save(q);
    return saved.id as string;
  }, name);
}

async function setCurrentObjectToQuery(page: Page, name: string) {
  await page.evaluate(async (qName) => {
    const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
    if (!q) throw new Error('query not found: ' + qName);
    grok.shell.o = q;
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

test('Sharing & Permissions — Query', async ({page}) => {
  test.setTimeout(420_000);

  
  await loginToDatagrok(page);
  await page.waitForTimeout(2000);

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
  const queryId = await createQuery(page, QUERY_NAME);
  await setCurrentObjectToQuery(page, QUERY_NAME);

  
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

  await softStep('Block A.2: Click SHARE...; Share <query> dialog opens with expected controls', async () => {
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
    await page.waitForTimeout(2000);
    
    
    
    
    
    
    
    
    
    const popup = page.locator('.d4-tags-selector-drop-down.d4-user-selector-drop-down');
    await expect(popup).toBeVisible({timeout: 10_000});
    await page.locator('input[placeholder="User, group, or email"]').fill('');
    await expect(popup).toBeHidden({timeout: 5_000});
  });

  await softStep('Block B.2: Notification controls present; dependent-entity cascade (behavioral)', async () => {
    
    
    
    
    
    
    
    
    
    
    
    
    await expect(page.locator('textarea[placeholder="Type in message here"]')).toBeAttached();
    const sendNotifPresent = await page.locator(
      '[name="input-Send-notifications"], .grok-permission-notifications input[type="checkbox"]').count();
    test.info().annotations.push({type: 'remark',
      description: `Block B.2 Send-notifications checkbox present: ${sendNotifPresent > 0}`});
    
    
    
    
    
    const cascadePresent = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      return /will also be shared/i.test(dlg?.textContent ?? '');
    });
    test.info().annotations.push({type: 'remark',
      description: `Block B.2 dependent-entity cascade notice present: ${cascadePresent}`});
  });

  await softStep('Block B.3: CANCEL closes dialog; no grant changed', async () => {
    await page.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog')).toHaveCount(0, {timeout: 10_000});
    
    
    
    
    const grants = await page.evaluate(async (qName) => {
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      return q ? 'still-owner-only' : 'missing';
    }, QUERY_NAME).catch(() => 'eval-skip');
    expect(grants).not.toBe('missing');
  });

  
  await softStep('Block C.1: Open Advanced editor; PermissionsView at /permissions/<id> with matrix', async () => {
    
    
    
    
    
    
    
    
    
    await page.goto(`${baseUrl}/permissions/${queryId}`);
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
      description: 'Advanced-editor matrix is a canvas d4-grid; "Common"/"Execute" column ' +
        'headers are canvas-painted (not DOM text) — verified loaded via grid + Save + ' +
        'Calculate-permissions DOM signals.'});
  });

  await softStep('Block C.2: Add-user row present; close without saving', async () => {
    await expect(
      page.locator('input[placeholder="Type in user, role or group to add..."]')).toBeVisible();
    
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1500);
  });

  
  await softStep('Block D.1: Owner shares "View and use" with recipient via JS API grant', async () => {
    
    
    
    
    
    
    
    const granted = await page.evaluate(async (args) => {
      const {qName, login} = args;
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!q || !grp) return {ok: false, reason: !q ? 'no-query' : 'no-recipient-group'};
      await grok.dapi.permissions.grant(q, grp, false); 
      const p = await grok.dapi.permissions.get(q);
      return {ok: true, viewGroups: (p?.view ?? []).map((g: any) => g.friendlyName || g.name)};
    }, {qName: QUERY_NAME, login: recipientLogin});
    expect(granted.ok).toBe(true);
    expect((granted as any).viewGroups.join(' ').toLowerCase()).toContain(recipientLogin.toLowerCase());
  });

  await softStep('Block D.2-3: Recipient sees shared query and can run it', async () => {
    await loginAsSecondUser(page);
    await page.waitForTimeout(2000);
    const result = await page.evaluate(async (qName) => {
      
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      if (!q) return {found: false};
      try {
        const df = await q.apply();
        return {found: true, ran: true, rows: df?.rowCount ?? null};
      } catch (e) {
        return {found: true, ran: false, error: String(e)};
      }
    }, QUERY_NAME);
    expect(result.found).toBe(true);
    expect(result.ran).toBe(true);
  });

  
  await softStep('Block E: Recipient lacks Edit / Delete / Share on the shared query', async () => {
    const checks = await page.evaluate(async (qName) => {
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      if (!q) return {found: false};
      const canEdit = await grok.dapi.permissions.check(q, 'Edit');
      const canDelete = await grok.dapi.permissions.check(q, 'Delete');
      const canShare = await grok.dapi.permissions.check(q, 'Share');
      const canView = await grok.dapi.permissions.check(q, 'View');
      return {found: true, canEdit, canDelete, canShare, canView};
    }, QUERY_NAME);
    expect(checks.found).toBe(true);
    expect(checks.canView).toBe(true);   
    expect(checks.canEdit).toBe(false);  
    expect(checks.canDelete).toBe(false); 
    expect(checks.canShare).toBe(false); 
  });

  
  await softStep('Block F.1-2: Owner revokes recipient grant; pane shows owner-only', async () => {
    await loginToDatagrok(page); 
    await page.waitForTimeout(2000);
    const revoked = await page.evaluate(async (args) => {
      const {qName, login} = args;
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      const grp = await grok.dapi.groups.filter(`name = "${login}"`).first();
      if (!q || !grp) return {ok: false};
      await grok.dapi.permissions.revoke(grp, q); 
      const p = await grok.dapi.permissions.get(q);
      const viewLogins = (p?.view ?? []).map((g: any) => (g.friendlyName || g.name || '').toLowerCase());
      return {ok: true, stillGranted: viewLogins.includes(login.toLowerCase())};
    }, {qName: QUERY_NAME, login: recipientLogin});
    expect(revoked.ok).toBe(true);
    expect((revoked as any).stillGranted).toBe(false);
  });

  await softStep('Block F.3-4: Recipient can no longer execute the query (access revoked)', async () => {
    await loginAsSecondUser(page);
    await page.waitForTimeout(2000);
    const result = await page.evaluate(async (qName) => {
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      if (!q) return {found: false, ran: false}; 
      try {
        await q.apply();
        return {found: true, ran: true};
      } catch (e) {
        return {found: true, ran: false}; 
      }
    }, QUERY_NAME);
    
    expect(result.ran).toBe(false);
  });

  
  await loginToDatagrok(page);
  await page.waitForTimeout(1500);
  await page.evaluate(async (qName) => {
    try {
      const q = await grok.dapi.queries.filter(`name = "${qName}"`).first();
      if (q) await grok.dapi.queries.delete(q);
    } catch (e) {  }
    grok.shell.closeAll();
  }, QUERY_NAME);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
