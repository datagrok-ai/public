/* ---
sub_features_covered: [projects.api.namespaces, projects.shell.open, projects.shell.share-via-context-menu]
--- */
// GROK-18345: share a project with a second user at View-and-Use + Full, then verify recipient can open it.
// Share grant uses grok.dapi.permissions.grant; recipient re-auth via loginAsSecondUser (needs DATAGROK_AUTH_TOKEN_2).
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser, getSecondUserLogin, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';

test.use(projectsTestOptions);

// Two-user test: second user is MANDATORY. loginAsSecondUser THROWS if no token2 is configured.
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
  // Defensive auto-Share dialog handling.
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
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = 'AutoTest-ComplexShare-' + stamp;

  await gotoApp(page);
  await setupSession(page);

  const ownerLogin = await readLogin(page);
  // Learn the second user's login from the token claim (no page switch). THROWS if no second user configured.
  const secondLogin = await getSecondUserLogin();
  console.log(`[two-user] owner='${ownerLogin}', second-user='${secondLogin}' (from token claim, no page switch)`);
  expect(secondLogin, 'second-user login must resolve').toBeTruthy();
  expect(secondLogin, 'second user must differ from owner').not.toBe(ownerLogin);

  // The login the project was actually granted to (must be the second user).
  let recipientLogin: string | null = null;

  try {
    await softStep('Setup: build project (file source, Sync ON)', async () => {
      await closeAll(page);
      // INTENTIONAL legacy `grok.data.files.openTable` (no `.script` tag) — keeps the toolbar SAVE button
      // visible for the UI Save dialog flow. `openTableFromFile` writes a .script tag that triggers bug 2b
      // (SAVE button collapses to offsetWidth=0). Don't migrate without first moving saveProject off the UI dialog.
      await evalJs(page, `(async () => {
        const df = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        const tv = grok.shell.addTableView(df);
        tv.addViewer('Scatter plot');
      })()`);
      await page.waitForTimeout(2000);
      await saveProject(page, projectName);
      // 60s poll: dev save commits can lag under load.
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
          // Permissions attach to GROUPS — grant to the user's personal group (User grant violates FK).
          await grok.dapi.permissions.grant(p, target.group, false);
          return { skipped: false, login: target.login, groupId: target.group.id, level: 'view-and-use' };
        } catch (e) {
          return { skipped: true, reason: String(e).slice(0, 200) };
        }
      })()`);
      // The grant must succeed — a thrown error (FK violation, missing project,
      // permissions API change) is a real regression, not a skip.
      expect(result.skipped, result.skipped ? `Step 12 grant failed: ${result.reason}` : '').toBe(false);
      expect(result.login).toBeTruthy();
      recipientLogin = result.login;

      // Verify the recipient's group appears in the project's granted permissions.
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
      // permissions.get returning an error is a real failure to verify the grant.
      expect(listed.ok, listed.ok ? '' : `permissions.get failed: ${listed.err}`).toBe(true);
      expect(listed.contains, 'recipient group must appear in the project granted permissions').toBe(true);
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
      // Elevation to Full must succeed — a thrown error is a real regression.
      expect(result.skipped, result.skipped ? `Step 12 (Full) elevation failed: ${result.reason}` : '').toBe(false);
      expect(result.login).toBeTruthy();
      console.log(`[two-user] Step 12: elevated to Full access for '${result.login}' — share complete, STILL OWNER (no user switch performed yet)`);
    });

    // Step 13: second user logs in (token2), locates + opens the shared project (GROK-18345 / GROK-19403).
    // The share MUST have targeted the second user — assert it, do not skip.
    expect(recipientLogin, 'Step 12 must have granted to the second user').toBe(secondLogin);
    await softStep('Step 13: second user logs in, sees AND opens the shared project', async () => {
      await loginAsSecondUser(page);
      try {
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

        // (b) the recipient can OPEN it — materialized tables / non-null shell.project prove cross-user open.
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

  finishSpec();
});
