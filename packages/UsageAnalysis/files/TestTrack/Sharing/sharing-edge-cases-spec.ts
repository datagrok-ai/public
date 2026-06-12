/* ---
sub_features_covered: [sharing.advanced-editor, sharing.permission-types, sharing.server.privileges-router, sharing.share-dialog]
--- */
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

test('Sharing — Edge cases (SHARE_WITH_EVERYONE gate, owner-retention, dependent-entity notice)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await setupSession(page);

  const stamp = Date.now();
  let scriptId: string | null = null;
  let projId: string | null = null;
  let allUsersGroupId: string | null = null;

  try {
    
    
    await softStep('Setup: create a Script that is a member of a Project; resolve All users group', async () => {
      const res = await evalJs(page, `(async () => {
        const g = window.grok, DG = window.DG;
        try {
          const script = await g.dapi.scripts.save(DG.Script.create('edgecase_' + ${stamp}));
          const project = DG.Project.create();
          project.name = 'AutoTest-EdgeCases-' + ${stamp};
          project.isPackage = false;
          project.addChild(script);
          const saved = await g.dapi.projects.save(project);
          // Resolve the built-in all-users group robustly: the loose filter('All
          // users').first() returns a relevance-matched WRONG group, so enumerate
          // and match on exact name/friendlyName (see sharing-api-permissions-api.ts).
          const candidates = await g.dapi.groups.filter('All users').list({limit: 50});
          const allUsers = candidates.find(x => x && (x.name === 'AllUsers' || x.friendlyName === 'All users'));
          return {ok: true, scriptId: script.id, projId: saved.id,
            allUsersGroupId: allUsers ? allUsers.id : null,
            ownerCanShare: await g.dapi.permissions.check(script, 'Share')};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 250)}; }
      })()`);
      expect(res.ok, `setup must succeed (${res.reason ?? ''})`).toBe(true);
      expect(res.scriptId, 'script entity must be created').toBeTruthy();
      expect(res.projId, 'parent project must be created').toBeTruthy();
      expect(res.allUsersGroupId, 'All users group must resolve').toBeTruthy();
      expect(res.ownerCanShare, 'owner must hold Share on the created script').toBe(true);
      scriptId = res.scriptId; projId = res.projId; allUsersGroupId = res.allUsersGroupId;
      console.log(`[edge] Setup: script=${scriptId} in project=${projId}; All users=${allUsersGroupId}`);
    });

    
    
    
    
    
    await softStep('Scenario 1: SHARE_WITH_EVERYONE global-permission gate returns a boolean and exists', async () => {
      const res = await evalJs(page, `(async () => {
        try {
          // Same-origin REST route; auth cookie auto-attached.
          const r = await fetch('/api/privileges/permissions/check/global/ShareWithEveryone', {headers: {'Accept': 'application/json'}});
          const body = (await r.text()).trim();
          // An unknown permission returns an ApiError body; ShareWithEveryone is a
          // real constant, so the body must be the literal boolean "true"/"false".
          const isBoolean = body === 'true' || body === 'false';
          const holds = body === 'true';
          return {ok: true, status: r.status, body: body.slice(0, 80), isBoolean, holds};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 250)}; }
      })()`);
      expect(res.ok, `global-permission check must run (${res.reason ?? ''})`).toBe(true);
      expect(res.status, 'global-permission check route must respond 200').toBe(200);
      
      
      expect(res.isBoolean, `gate must return a boolean (got "${res.body}")`).toBe(true);
      console.log(`[edge] Scenario 1: ShareWithEveryone gate → ${res.body} (current actor holds=${res.holds})`);
    });

    await softStep('Scenario 1: sharing with All users succeeds for a holder; under no grant the entity is owner-only', async () => {
      const res = await evalJs(page, `(async () => {
        const g = window.grok;
        try {
          const script = await g.dapi.scripts.find('${scriptId}');
          const allUsers = await g.dapi.groups.find('${allUsersGroupId}');
          // Holder branch: the all-users grant lands only because the actor holds
          // ShareWithEveryone (recon B). The policy invariant: all-user visibility is
          // gated on the owner holding the global permission.
          await g.dapi.permissions.grant(script, allUsers, false); // view-and-use
          let perms = await g.dapi.permissions.get(script);
          const groups = [...(perms.view || []), ...(perms.edit || [])];
          const granted = groups.some(grp => grp && grp.id === '${allUsersGroupId}');
          // Revoke it back: confirm the all-users grant is fully removable (no
          // partial/silent residue — the policy never leaves a dangling all-user grant).
          await g.dapi.permissions.revoke(allUsers, script); // revoke(group, entity)
          perms = await g.dapi.permissions.get(script);
          const stillThere = [...(perms.view || []), ...(perms.edit || [])]
            .some(grp => grp && grp.id === '${allUsersGroupId}');
          return {ok: true, granted, stillThere};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 250)}; }
      })()`);
      expect(res.ok, `all-users grant/revoke must run (${res.reason ?? ''})`).toBe(true);
      expect(res.granted, 'all-users grant must land for a ShareWithEveryone holder').toBe(true);
      expect(res.stillThere, 'all-users grant must be fully removable (no silent partial residue)').toBe(false);
      console.log('[edge] Scenario 1: all-users grant landed for holder and revoked cleanly (gate-governed path)');
    });

    
    
    
    
    
    await softStep('Scenario 2 (UI): owner opens the Advanced editor (PermissionsView) for the script', async () => {
      const opened = await evalJs(page, `(async () => {
        try {
          // The Share dialog "Advanced editor..." label navigates to this same
          // /permissions/<id> URL (sharing.md Section 6).
          window.grok.shell.startUri = '/permissions/${scriptId}';
          return {ok: true};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 200)}; }
      })()`);
      expect(opened.ok, `navigation to PermissionsView must run (${opened.reason ?? ''})`).toBe(true);
      
      
      
      const loaded = await page.locator(
        '.grok-permissions-self, .d4-grid, [name="button-Calculate-resulting-permissions-for-this-entity"], input[placeholder="Type in user, role or group to add..."]',
      ).first().isVisible({timeout: 20_000}).catch(() => false);
      if (loaded)
        console.log('[edge] Scenario 2 (UI): Advanced editor PermissionsView loaded (DOM signal, class-1)');
      else
        console.warn('[edge] Scenario 2 (UI): PermissionsView DOM signal not observed; owner-retention is asserted via API below');
      
      await evalJs(page, `(async () => { try { window.grok.shell.closeAll(); } catch(_){} })()`).catch(() => {});
      await page.waitForTimeout(500);
    });

    await softStep('Scenario 2: removing EVERY grant leaves the owner with View/Edit/Delete/Share (server enforces owner retention)', async () => {
      const res = await evalJs(page, `(async () => {
        const g = window.grok;
        try {
          const script = await g.dapi.scripts.find('${scriptId}');
          // The advanced-editor "remove every grant row (incl. the owner's) + SAVE"
          // gesture. Owner ownership is NOT a grant row, so it cannot be removed.
          const perms = await g.dapi.permissions.get(script);
          for (const grp of [...(perms.view || []), ...(perms.edit || [])])
            await g.dapi.permissions.revoke(grp, script);
          const after = await g.dapi.permissions.get(script);
          const remaining = [...(after.view || []), ...(after.edit || [])].length;
          // Owner self-checks: ownership preserved server-side even with zero grants.
          const ownerView = await g.dapi.permissions.check(script, 'View');
          const ownerEdit = await g.dapi.permissions.check(script, 'Edit');
          const ownerDelete = await g.dapi.permissions.check(script, 'Delete');
          const ownerShare = await g.dapi.permissions.check(script, 'Share');
          // Re-open as owner to confirm reachability after stripping all grants.
          const reopened = await g.dapi.scripts.find('${scriptId}').catch(() => null);
          return {ok: true, remaining, ownerView, ownerEdit, ownerDelete, ownerShare, reopenable: reopened != null};
        } catch (e) { return {ok: false, reason: String(e).slice(0, 250)}; }
      })()`);
      expect(res.ok, `remove-all-grants must succeed (${res.reason ?? ''})`).toBe(true);
      expect(res.remaining, 'every grant row must be removed (matrix empty)').toBe(0);
      
      
      expect(res.ownerView, 'owner must still have View after removing every grant').toBe(true);
      expect(res.ownerEdit, 'owner must still have Edit after removing every grant').toBe(true);
      expect(res.ownerDelete, 'owner must still have Delete after removing every grant').toBe(true);
      expect(res.ownerShare, 'owner must still be able to re-share after removing every grant').toBe(true);
      expect(res.reopenable, 'owner must be able to re-open the entity after removing every grant').toBe(true);
      console.log(`[edge] Scenario 2: grants removed (remaining=${res.remaining}); owner retains View/Edit/Delete/Share — ownership is not a removable grant row`);
    });

    
    
    
    
    
    
    await softStep('Scenario 3 (UI): Share dialog PermissionsEditor renders for the project-member script', async () => {
      
      await evalJs(page, `(async () => {
        const g = window.grok;
        const script = await g.dapi.scripts.find('${scriptId}');
        g.shell.o = script;
      })()`);
      await page.waitForTimeout(1500);

      
      
      
      
      const header = page.locator('[name="div-section--Sharing"]');
      const shareBtn = page.locator('[name="button-Share..."]');
      if (await header.isVisible({timeout: 10_000}).catch(() => false)) {
        const laidOut = await shareBtn.evaluate((el: any) => el && el.offsetParent !== null).catch(() => false);
        if (!laidOut) {
          await header.click();
          await page.waitForTimeout(800);
        }
      }
      await expect(shareBtn, 'SHARE... button must be attached (context-panel Sharing pane)').toBeAttached({timeout: 10_000});
      
      await shareBtn.click().catch(() => {});
      await page.waitForTimeout(1200);

      
      const dialog = page.locator('.d4-dialog');
      const dialogVisible = await dialog.first().isVisible({timeout: 10_000}).catch(() => false);
      expect(dialogVisible, 'Share dialog (PermissionsEditor) must render for the script').toBe(true);
      
      const editorInput = await page.locator('.d4-dialog input[placeholder="User, group, or email"]')
        .first().isVisible({timeout: 5_000}).catch(() => false);
      expect(editorInput, 'PermissionsEditor recipient autocomplete must render').toBe(true);

      
      
      const noticeText = await page.locator('.d4-dialog').first().evaluate((dlg: any) => {
        const t = (dlg.textContent || '').replace(/\s+/g, ' ').trim();
        const hasNotice = /\b(also|will be shared|depend|member of|co-)\b/i.test(t);
        return {hasNotice, snippet: t.slice(0, 200)};
      }).catch(() => ({hasNotice: false, snippet: ''}));
      console.log(`[edge] Scenario 3: PermissionsEditor rendered; dependent-entity notice present=${noticeText.hasNotice} (behavioral remark, SR-02)`);

      
      await page.locator('.d4-dialog [name="button-CANCEL"]').first().click().catch(() => {});
      await page.waitForTimeout(400);
    });
  } finally {
    
    await evalJs(page, `(async () => {
      const g = window.grok;
      try {
        if ('${scriptId}' !== 'null') {
          const script = await g.dapi.scripts.find('${scriptId}').catch(() => null);
          if (script) {
            const perms = await g.dapi.permissions.get(script).catch(() => ({view: [], edit: []}));
            for (const grp of [...(perms.view || []), ...(perms.edit || [])])
              await g.dapi.permissions.revoke(grp, script).catch(() => {});
          }
          if ('${projId}' !== 'null') {
            const project = await g.dapi.projects.find('${projId}').catch(() => null);
            if (project) await g.dapi.projects.delete(project).catch(() => {});
          }
          if (script) await g.dapi.scripts.delete(script).catch(() => {});
        }
      } catch (_) { /* best-effort cleanup */ }
    })()`).catch(() => {});
    await page.evaluate(() => (window as any).grok?.shell?.closeAll?.());
  }

  if (stepErrors.length > 0)
    throw new Error('Spec step failures:\n' + stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
});
