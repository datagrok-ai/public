import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, expr: string): Promise<any> {
  return page.evaluate((e) => (0, eval)(e), expr);
}

async function setupSession(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });
}

async function openShareDialogViaPane(page: Page, projId: string) {
  await evalJs(page, `(async () => {
    const g = window.grok;
    const project = await g.dapi.projects.find('${projId}');
    g.shell.o = project;
    g.shell.windows.showContextPanel = true;
  })()`);
  await page.locator('[name="div-section--Sharing"]').waitFor({state: 'attached', timeout: 20_000});
  const shareBtn = page.locator('[name="button-Share..."]');
  const laidOut = await shareBtn.evaluate((el: any) => el && el.offsetParent !== null).catch(() => false);
  if (!laidOut) {
    await page.locator('[name="div-section--Sharing"]').click();
    await page.waitForTimeout(800);
  }
  await expect(shareBtn, 'SHARE... button must be present in the Sharing pane').toBeAttached({timeout: 10_000});
  await shareBtn.click();
  await page.locator('.d4-dialog').waitFor({state: 'visible', timeout: 15_000});
}

test('Sharing — UI Smoke (Single-Actor): share dialog, context panel, advanced editor', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await setupSession(page);

  const stamp = Date.now();
  const projectName = 'AutoTest-SharingSmoke-' + stamp;
  let projId: string | null = null;

  try {
    
    
    
    
    await softStep('Setup: create an owned, persisted project to share', async () => {
      const created = await evalJs(page, `(async () => {
        const g = window.grok, DG = window.DG;
        try {
          const df = grok.data.demo.demog(30);
          df.name = 'demog_' + ${stamp};
          const tableId = await g.dapi.tables.uploadDataFrame(df); // returns a STRING id
          const isBound = (e) => !!e && e.dart !== undefined && e.id === tableId;
          let tableInfo = null;
          for (let i = 0; i < 30; i++) {
            tableInfo = await g.dapi.tables.find(tableId).catch(() => null);
            if (isBound(tableInfo)) break;
            await new Promise(r => setTimeout(r, 500));
          }
          if (!isBound(tableInfo))
            return {ok: false, reason: 'tables.find did not return a Dart-bound TableInfo within 15s (cold-init race)'};
          const project = DG.Project.create();
          project.name = ${JSON.stringify(projectName)};
          project.isPackage = false;
          project.addChild(tableInfo);
          const saved = await g.dapi.projects.save(project);
          if (!saved || !saved.id) return {ok: false, reason: 'projects.save returned no id'};
          return {ok: true, projId: saved.id, ownerCanShare: await g.dapi.permissions.check(saved, 'Share')};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 250)}; }
      })()`);
      expect(created.ok, `project creation must succeed (${created.reason ?? ''})`).toBe(true);
      expect(created.ownerCanShare, 'owner must hold Share on the created project').toBe(true);
      projId = created.projId;
      console.log(`[sharing-smoke] Setup: project '${projectName}' (${projId}) created`);
    });

    
    
    
    await softStep('Scenario 3: Context Panel Sharing pane renders grant list + SHARE... button', async () => {
      await evalJs(page, `(async () => {
        const g = window.grok;
        const project = await g.dapi.projects.find('${projId}');
        g.shell.o = project;
        g.shell.windows.showContextPanel = true;
      })()`);
      const sharingHeader = page.locator('[name="div-section--Sharing"]');
      await expect(sharingHeader, 'Sharing context-panel pane must be present for the owned project').toBeAttached({timeout: 20_000});

      
      const shareBtn = page.locator('[name="button-Share..."]');
      const laidOut = await shareBtn.evaluate((el: any) => el && el.offsetParent !== null).catch(() => false);
      if (!laidOut) {
        await sharingHeader.click();
        await page.waitForTimeout(800);
      }
      await expect(shareBtn, 'SHARE... button must render in the Sharing pane').toBeAttached({timeout: 10_000});
      
      await expect(page.locator('[name="button-Calculate-resulting-permissions-for-this-entity"]'),
        'Calculate-permissions button confirms the Sharing pane fully rendered').toBeAttached({timeout: 10_000});
      console.log('[sharing-smoke] Scenario 3: Sharing pane + SHARE... + Calculate buttons present (DOM-driven, class-1)');
    });

    
    await softStep('Scenario 3 (cont.): SHARE... button opens the Share dialog modal', async () => {
      await page.locator('[name="button-Share..."]').click();
      const dlg = page.locator('.d4-dialog');
      await expect(dlg, 'Share dialog modal must open from the Sharing pane SHARE... button').toBeVisible({timeout: 15_000});
      await expect(page.locator('.d4-dialog .d4-dialog-title'),
        'dialog title must read "Share <EntityName>"').toContainText('Share', {timeout: 5_000});
      
      await page.locator('[name="button-CANCEL"]').click();
      await expect(dlg, 'dialog must close on CANCEL').toBeHidden({timeout: 10_000});
      console.log('[sharing-smoke] Scenario 3 (cont.): SHARE... opened Share dialog; CANCEL closed it');
    });

    
    
    
    await softStep('Scenario 1: Share dialog hosts the full PermissionsEditor widget', async () => {
      await openShareDialogViaPane(page, projId!);
      const dlg = page.locator('.d4-dialog');
      
      await expect(dlg.locator('input[placeholder="User, group, or email"]'),
        'recipient autocomplete input must be present').toBeVisible({timeout: 10_000});
      
      await expect(dlg.locator('[name="div-share-selector"]'),
        'access-level dropdown (View and use / Full access) must be present').toBeAttached({timeout: 10_000});
      
      await expect(dlg.locator('[name="label-Advanced-editor..."]'),
        'Advanced editor... link must be present').toBeAttached({timeout: 10_000});
      
      
      
      
      await expect(dlg.locator('textarea[placeholder="Type in message here"]'),
        'notification message textarea must be present').toBeAttached({timeout: 10_000});
      await expect(dlg.locator('[name="input-Send-notifications"]'),
        'Send notifications checkbox must be present').toBeAttached({timeout: 10_000});
      
      await expect(dlg.locator('[name="button-OK"]'), 'OK button must be present').toBeAttached({timeout: 10_000});
      await expect(dlg.locator('[name="button-CANCEL"]'), 'CANCEL button must be present').toBeAttached({timeout: 10_000});
      console.log('[sharing-smoke] Scenario 1: PermissionsEditor widget fully rendered (input, dropdown, advanced link, notify, OK/CANCEL)');
      
    });

    
    
    
    await softStep('Scenario 2: type into autocomplete, see suggestions, CANCEL applies no change', async () => {
      const dlg = page.locator('.d4-dialog');
      await expect(dlg, 'Share dialog must still be open from Scenario 1').toBeVisible({timeout: 5_000});

      
      const before = await evalJs(page, `(async () => {
        const g = window.grok;
        const p = await g.dapi.projects.find('${projId}');
        const perms = await g.dapi.permissions.get(p);
        return {groupCount: [...(perms.view||[]), ...(perms.edit||[])].length};
      })()`);

      
      const input = dlg.locator('input[placeholder="User, group, or email"]');
      await input.click();
      await input.fill('');
      await page.keyboard.type('adm', {delay: 60});
      
      
      const popup = page.locator('.d4-tags-selector-drop-down.d4-user-selector-drop-down');
      await expect(popup, 'autocomplete suggestion popup must appear after typing').toBeVisible({timeout: 10_000});
      const suggestionCount = await page.locator('span.d4-user-selector-user-name').count();
      expect(suggestionCount, 'at least one user/group suggestion must be listed').toBeGreaterThan(0);
      console.log(`[sharing-smoke] Scenario 2: autocomplete popup visible with ${suggestionCount} suggestion(s)`);

      
      await dlg.locator('[name="button-CANCEL"]').click();
      await expect(dlg, 'dialog must close on CANCEL').toBeHidden({timeout: 10_000});

      
      const after = await evalJs(page, `(async () => {
        const g = window.grok;
        const p = await g.dapi.projects.find('${projId}');
        const perms = await g.dapi.permissions.get(p);
        return {groupCount: [...(perms.view||[]), ...(perms.edit||[])].length};
      })()`);
      expect(after.groupCount, 'CANCEL must apply no permission change (grant count unchanged)').toBe(before.groupCount);
      console.log(`[sharing-smoke] Scenario 2: CANCEL applied no change (grants before=${before.groupCount}, after=${after.groupCount})`);
    });

    
    
    
    
    
    
    await softStep('Scenario 4: Advanced editor... opens the PermissionsView matrix', async () => {
      await openShareDialogViaPane(page, projId!);
      const dlg = page.locator('.d4-dialog');
      const advLabel = dlg.locator('[name="label-Advanced-editor..."]');
      await expect(advLabel, 'Advanced editor... link must be present').toBeAttached({timeout: 10_000});
      await advLabel.click();
      
      
      
      const loaded = page.locator(
        '[name="button-Calculate-resulting-permissions-for-this-entity"], ' +
        'input[placeholder="Type in user, role or group to add..."], ' +
        '[name="button-Save"], .d4-grid',
      ).first();
      await expect(loaded, 'PermissionsView (Advanced editor) must load with its matrix surface').toBeVisible({timeout: 25_000});
      console.log('[sharing-smoke] Scenario 4: Advanced editor PermissionsView loaded (DOM signal, class-1)');

      
      await evalJs(page, `(async () => { try { window.grok.shell.closeAll(); } catch(_){} })()`).catch(() => {});
      await page.waitForTimeout(500);
      const after = await evalJs(page, `(async () => {
        const g = window.grok;
        const p = await g.dapi.projects.find('${projId}').catch(() => null);
        if (!p) return {groupCount: -1};
        const perms = await g.dapi.permissions.get(p);
        return {groupCount: [...(perms.view||[]), ...(perms.edit||[])].length};
      })()`);
      expect(after.groupCount, 'navigating away from Advanced editor without SAVE must alter no permissions').toBe(0);
      console.log('[sharing-smoke] Scenario 4: left Advanced editor without SAVE; permissions unchanged');
    });
  } finally {
    
    await evalJs(page, `(async () => {
      const g = window.grok;
      try {
        if ('${projId}' !== 'null' && '${projId}' !== '') {
          const project = await g.dapi.projects.find('${projId}').catch(() => null);
          if (project) {
            const perms = await g.dapi.permissions.get(project).catch(() => ({view: [], edit: []}));
            for (const grp of [...(perms.view||[]), ...(perms.edit||[])])
              await g.dapi.permissions.revoke(grp, project).catch(() => {});
            await g.dapi.projects.delete(project).catch(() => {});
          }
        }
      } catch (_) { /* best-effort cleanup */ }
    })()`).catch(() => {});
    await page.evaluate(() => (window as any).grok?.shell?.closeAll?.());
  }

  if (stepErrors.length > 0)
    throw new Error('Spec step failures:\n' + stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
});
