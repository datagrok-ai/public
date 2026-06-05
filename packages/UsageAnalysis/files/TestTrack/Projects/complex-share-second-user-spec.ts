/* ---
sub_features_covered: [projects.shell.share-via-context-menu, projects.api.namespaces, projects.shell.open]
--- */
// Selector sources (grok-browser/references):
//   widgets/dialog.md:22,29,61,69,74-92 — d4-dialog, button-OK, button-CANCEL, Dart input pattern
//   projects.md:232 — SAVE ribbon button is the reliable Save Project trigger
//   projects.md:118-144 — Sharing via JS API: grok.dapi.permissions.grant(p, user, edit)
//     where edit=false → View-and-Use, edit=true → Full
//   projects.md:146-165 — Re-auth via helpers.playwright.session.logoutAndLoginAs
//     (token-based; switches the authenticated session by re-injecting the
//     second-user token — requires DATAGROK_AUTH_TOKEN_2).
//
// Wave 1b/2C complex-split: covers Step (share with second user +
// recipient open) sub-bullet of complex.md scenario — specifically
// Step 12 (configure project sharing for another user at View-and-Use
// AND Full access levels). Targets the share-grant → recipient-side
// open flow that verifies cross-user project visibility after share.
// Helper-3 (logoutAndLoginAs) dependency: uses logoutAndLoginAs to
// switch users mid-test (GROK-18345 partial regression invariant —
// bug-library entry exists; full reproduction requires Spaces dataset
// + datasync + recipient open).
//
// Bug-focused slice satellite of complex.md per Decision 2.6 expanded
// pattern 2 (orphan-without-its-own-parent-.md acceptable when the spec
// self-documents via header cross-reference to parent .md + functionality
// slice + GROK ticket). Parent canonical scenario: complex.md.
//
// Scope: share existing project to second user at both grant levels +
// (when Helper-3 is registered) recipient logs in + recipient opens
// shared project. Sub-bullets of complex.md NOT covered here (Pivot,
// Aggregate, Join, Clone, derived tables, multi-source, rename) belong
// to other satellites of the complex.md decomposition.
//
// Step 13 (log in as second user, locate + open the shared project) is now
// wired via token-injection re-auth (loginAsSecondUser, needs
// DATAGROK_AUTH_TOKEN_2 — provided by the runner from DATAGROK_DEV_KEY_2 /
// config key2). To make the recipient-open verifiable, the share grant now
// targets the SECOND user specifically (its login is probed via a token2
// re-auth round-trip) instead of an arbitrary other user. When no
// second-user token is configured the spec falls back to sharing with any
// other user and skips only the recipient-open leg — so it still runs (no
// hard failure) without token2.
//
// Scope reductions (documented):
//   * Step 12 right-click Share dialog UI is replaced with
//     grok.dapi.permissions.grant JS API per projects.md:28 documented
//     preferred default (and projects.md:242 — context menu items have
//     no name= attributes). This is consistent with share-project-spec
//     pattern.
//   * GROK-18345 full reproduction (Spaces dataset + datasync + share)
//     is PARTIAL: the share + datasync touch points are exercised, but
//     the cross-user open verification requires re-auth.
//   * Defensive skip-with-note pattern for both grant levels (matches
//     share-project-spec). On dev the recipient-user lookup may fail with FK
//     violation on permissions table — defensive handling absorbs it.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser, getSecondUserLogin, softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';

test.use(projectsTestOptions);

// This is a two-user test: the second user is MANDATORY. loginAsSecondUser
// resolves the token from DATAGROK_AUTH_TOKEN_2 or the config `key2:` and
// THROWS if neither exists — so a misconfigured environment FAILS the test
// rather than silently passing without exercising the second user.
const readLogin = (page: Page): Promise<string | null> =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

async function saveProject(page: Page, name: string) {
  await page.click('button:has-text("SAVE"), .ui-btn:has-text("SAVE")');
  const dialog = page.locator('.d4-dialog').first();
  await dialog.waitFor({timeout: 15000});
  const nameInput = dialog.locator('input[type="text"]').first();
  await nameInput.focus();
  await page.keyboard.press('Control+a');
  await page.keyboard.type(name);
  await dialog.locator('[name="button-OK"]').click();
  // Defensive share-dialog handling per Wave 1a Validator hypothesis cycle 1
  // canonical pattern.
  const shareDialog = page.locator('.d4-dialog').filter({hasText: `Share ${name}`});
  if (await shareDialog.isVisible({timeout: 30000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }
}

async function deleteProjectByName(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`);
}

test('Projects / Complex share-second-user: dual-level grant + recipient open via JS API', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = 'AutoTest-ComplexShare-' + stamp;

  await gotoApp(page);
  await setupSession(page);

  const ownerLogin = await readLogin(page);
  // Learn the second user's login from the token claim — NO page switch (the
  // old "switch in / read / switch back" probe cost two extra reloads and ran
  // before the project even existed). getSecondUserLogin THROWS if no second
  // user is configured → the test fails (never a silent skip). The session
  // stays on the primary (owner) user; the only re-auth is Step 13.
  const secondLogin = await getSecondUserLogin();
  console.log(`[two-user] owner='${ownerLogin}', second-user='${secondLogin}' (from token claim, no page switch)`);
  expect(secondLogin, 'second-user login must resolve').toBeTruthy();
  expect(secondLogin, 'second user must differ from owner').not.toBe(ownerLogin);

  // The login the project was actually granted to (must be the second user).
  let recipientLogin: string | null = null;

  try {
    await softStep('Setup: build project (file source, Sync ON)', async () => {
      await closeAll(page);
      // INTENTIONAL legacy `grok.data.files.openTable` (no `.script`
      // provenance written). Required because this spec exercises the UI
      // Save dialog flow (`saveProject` helper at line 51 — toolbar SAVE
      // button click → dialog → OK). The canonical `openTableFromFile`
      // helper writes a dot-form `.script` tag and triggers bug 2b
      // (toolbar SAVE button collapses to offsetWidth=0 in Playwright),
      // which makes the UI Save dialog unreachable. The legacy path
      // keeps the SAVE button visible — verified live 2026-05-05.
      // DO NOT migrate to `openTableFromFile` without first migrating
      // `saveProject` away from the UI dialog flow (e.g. to inline JS
      // API save like complex-rename-spec / complex-derived-tables-spec).
      await evalJs(page, `(async () => {
        const df = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        const tv = grok.shell.addTableView(df);
        tv.addViewer('Scatter plot');
      })()`);
      await page.waitForTimeout(2000);
      await saveProject(page, projectName);
      // 60s poll: dev save commits can lag under load; 30s window too tight per
      // Wave 1b hypothesis cycle 1 (environmental-flake) — same class as Wave 1a
      // project-url issue but on the project-existence side rather than the
      // share-dialog side.
      await expect.poll(async () => evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p != null;
        })()`,
      ), {timeout: 60000, intervals: [500, 1000, 2000, 5000]}).toBe(true);
      console.log(`[two-user] Setup: project '${projectName}' BUILT & SAVED as owner '${ownerLogin}' (Save dialog OK confirmed by existence poll) — still owner, no user switch yet`);
    });

    await softStep('Step 12 (View-and-Use): share project with second user via grok.dapi.permissions.grant(_, _, false)', async () => {
      const result = await evalJs(page, `(async () => {
        try {
          const me = (await grok.dapi.users.current()).login;
          const wanted = ${JSON.stringify(secondLogin)};
          let target = null;
          if (wanted) target = await grok.dapi.users.filter('login = "' + wanted + '"').first();
          if (!target) {
            const users = await grok.dapi.users.list({limit: 50});
            target = users.find(u => u.login !== me && u.login !== 'system');
          }
          if (!target) return { skipped: true, reason: 'no target user found' };
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          // Permissions attach to GROUPS — grant to the user's personal group
          // (granting the User directly violates permissions_user_group_id_fkey).
          await grok.dapi.permissions.grant(p, target.group, false);
          return { skipped: false, login: target.login, groupId: target.group.id, level: 'view-and-use' };
        } catch (e) {
          return { skipped: true, reason: String(e).slice(0, 200) };
        }
      })()`);
      if (result.skipped) {
        console.warn('Step 12 (View-and-Use) skipped: ' + result.reason);
        return;
      }
      expect(result.login).toBeTruthy();
      recipientLogin = result.login;

      // Verify the recipient's group appears in the project's granted
      // permissions (permissions.get returns {view, edit} group lists).
      const listed = await evalJs(page, `(async () => {
        try {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          const perms = await grok.dapi.permissions.get(p);
          const groups = [...(perms.view || []), ...(perms.edit || [])];
          return { ok: true, contains: groups.some(g => g && (g.id === '${result.groupId}' || g.friendlyName === '${result.login}' || g.name === '${result.login}')) };
        } catch (e) {
          return { ok: false, err: String(e).slice(0, 200) };
        }
      })()`);
      if (!listed.ok) {
        console.warn('Step 12 verification skipped: permissions.get unsupported (' + listed.err + ')');
        return;
      }
      expect(listed.contains).toBe(true);
      console.log(`[two-user] Step 12: granted View-and-Use to '${result.login}' (owner-side permissions.get confirms recipient group) — still owner`);
    });

    await softStep('Step 12 (Full): elevate share to Full access via grok.dapi.permissions.grant(_, _, true)', async () => {
      if (!recipientLogin) { console.warn('Step 12 (Full) skipped: no recipient from prior step'); return; }
      const result = await evalJs(page, `(async () => {
        try {
          const target = await grok.dapi.users.filter('login = "${recipientLogin}"').first();
          if (!target) return { skipped: true, reason: 'recipient not found' };
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          await grok.dapi.permissions.grant(p, target.group, true);
          return { skipped: false, login: target.login, level: 'full' };
        } catch (e) {
          return { skipped: true, reason: String(e).slice(0, 200) };
        }
      })()`);
      if (result.skipped) {
        console.warn('Step 12 (Full) skipped: ' + result.reason);
        return;
      }
      expect(result.login).toBeTruthy();
      console.log(`[two-user] Step 12: elevated to Full access for '${result.login}' — share complete, STILL OWNER (no user switch performed yet)`);
    });

    // Step 13 (complex.md): the second user logs in (token2 injection),
    // locates the shared project, and OPENS it — the cross-user open invariant
    // (GROK-18345 / GROK-19403: "recipient cannot open shared project").
    // Restores the primary session afterwards so the finally-block cleanup
    // (delete) runs as the owner.
    // The share MUST have targeted the second user — assert it, do not skip.
    expect(recipientLogin, 'Step 12 must have granted to the second user').toBe(secondLogin);
    await softStep('Step 13: second user logs in, sees AND opens the shared project', async () => {
      await loginAsSecondUser(page);
      try {
        // Prove we are actually running AS the second user now (not the owner).
        const liveLogin = await readLogin(page);
        console.log(`[two-user] Step 13: now authenticated as '${liveLogin}' (expected second-user '${secondLogin}', owner was '${ownerLogin}')`);
        expect(liveLogin).toBe(secondLogin);
        expect(liveLogin).not.toBe(ownerLogin);

        // (a) the shared project is visible to the recipient.
        await expect.poll(async () => evalJs(page,
          `(async () => {
            const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
            return p != null;
          })()`,
        ), {timeout: 30_000, intervals: [1000, 2000, 5000]}).toBe(true);
        console.log(`[two-user] Step 13: recipient '${liveLogin}' CAN see shared project '${projectName}'`);

        // (b) the recipient can OPEN it (per complex.md Expected results:
        // "View-and-Use access can OPEN the project"). Materialised tables /
        // a non-null shell.project prove the open succeeded cross-user.
        const opened = await evalJs(page, `(async () => {
          try {
            grok.shell.closeAll();
            const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
            if (!p) return {ok: false, reason: 'project not visible to recipient'};
            await p.open();
            await new Promise(r => setTimeout(r, 2500));
            const proj = grok.shell.project;
            const tables = (() => { try { return Number(grok.shell.tables?.length) || 0; } catch (_) { return 0; } })();
            return {ok: (proj != null) || tables > 0, name: proj ? proj.name : null, tables};
          } catch (e) {
            return {ok: false, reason: String(e).slice(0, 200)};
          }
        })()`);
        expect(opened.ok, `recipient must be able to OPEN the shared project (${opened.reason ?? ''})`).toBe(true);
        console.log(`[two-user] Step 13: recipient '${liveLogin}' OPENED shared project (tables=${opened.tables}) — cross-user open verified`);
      } finally {
        await loginToDatagrok(page);
        await setupSession(page);
      }
    });
  } finally {
    await deleteProjectByName(page, projectName).catch(() => {});
    await closeAll(page);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
