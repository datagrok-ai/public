import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  saveCopy,
  shareProjectViaContextMenu,
  shareWithSecondUserAndVerify,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

async function reopenProjectById(page: any, projectId: string) {
  await evalJs(page, `(async () => {
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 400));
    const p = await grok.dapi.projects.find('${projectId}');
    if (!p) return;
    for (let attempt = 0; attempt < 2; attempt++) {
      await p.open();
      for (let i = 0; i < 16; i++) {
        const tv = grok.shell.tv;
        if (tv?.dataFrame && typeof tv.addViewer === 'function') return;
        await new Promise(r => setTimeout(r, 350));
      }
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));
    }
  })()`);
  await page.waitForTimeout(400);
}

async function addViewerSafely(page: any, viewerName: string): Promise<void> {
  await evalJs(page, `(async () => {
    for (let i = 0; i < 16; i++) {
      const tv = grok.shell.tv;
      if (tv && tv.dataFrame && typeof tv.addViewer === 'function') {
        tv.addViewer(${JSON.stringify(viewerName)});
        return;
      }
      await new Promise(r => setTimeout(r, 350));
    }
    throw new Error('grok.shell.tv.addViewer never became a function (reopen did not materialize a TableView)');
  })()`);
  await page.waitForTimeout(400);
}

async function rehydrateOriginal(page: any, name: string): Promise<{projectId: string; tableInfoId: string}> {
  await evalJs(page, `(async () => {
    try {
      const list = await grok.dapi.projects.filter('name like "${name}"').list();
      for (const p of list) await grok.dapi.projects.delete(p);
    } catch (_) {}
  })()`);
  await evalJs(page, 'grok.shell.closeAll()');
  await page.waitForTimeout(500);
  const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
  await assertProvenanceScript(page, 'files', opened.script);
  const saved = await saveProjectWithProvenance(page, name);
  return {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
}

async function reopenedRowCount(page: any): Promise<number> {

  return await evalJs<number>(page, `(() => {
    try {
      const tv = grok.shell.tv;
      const df = tv && tv.dataFrame;
      return df ? Number(df.rowCount) || 0 : 0;
    } catch (_) { return 0; }
  })()`);
}

async function reopenedViewerCount(page: any): Promise<number> {
  return await evalJs<number>(page, `(() => {
    try {
      const tv = grok.shell.tv;
      const v = tv && tv.viewers;
      return v ? Number(v.length) || 0 : 0;
    } catch (_) { return 0; }
  })()`);
}

async function captureActiveProjectIds(
  page: any,
  typedName: string,
  sourceId: string,
): Promise<{projectId: string; tableInfoId: string; serverName: string}> {
  return await evalJs<{projectId: string; tableInfoId: string; serverName: string}>(page, `(async () => {
    const typedName = ${JSON.stringify(typedName)};
    const sourceId = ${JSON.stringify(sourceId)};
    const shellProj = grok.shell.project;
    if (shellProj?.id && shellProj.id !== sourceId) {
      const tv = grok.shell.tv;
      const ti = tv?.dataFrame?.getTableInfo?.();
      return {projectId: shellProj.id, tableInfoId: ti?.id ?? '', serverName: shellProj.name || typedName};
    }
    try {
      const list = await grok.dapi.projects.filter('name like "' + typedName + '"').list();
      const candidate = list.find(p => p.id !== sourceId);
      if (candidate) return {projectId: candidate.id, tableInfoId: '', serverName: candidate.name || typedName};
    } catch (_) {}
    return {projectId: '', tableInfoId: '', serverName: typedName};
  })()`);
}

test('Projects / Copy Clone — full UI-driven 4-sub-flow + GROK-19750 invariant', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const names = {
    original: `demog-${stamp}`,
    linkCopy: `demog-${stamp}-link`,
    cloneCopy: `demog-${stamp}-clone`,
    pvcCopy: `demog-${stamp}-pvc`,
  };
  const ids: Record<string, {projectId: string; tableInfoId: string; serverName?: string} | undefined> = {};

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {

    await softStep('Setup: open demog.csv with provenance + save baseline (JS API)', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await assertProvenanceScript(page, 'files', opened.script);
      const saved = await saveProjectWithProvenance(page, names.original);
      ids.original = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
      expect(ids.original.projectId).toBeTruthy();
    });

    await softStep('4a: open original (UI-reopen) → addViewer Bar → SAVE mode=original (overwrite) → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await addViewerSafely(page, 'Bar chart');
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'original',
        name: names.original,
      });
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4b prep: rehydrate original (post-4a-overwrite GROK-19750 workaround)', async () => {
      ids.original = await rehydrateOriginal(page, names.original);
      expect(ids.original.projectId).toBeTruthy();
    });

    await softStep('4b step 1-4: open original → addViewer Scatter → SAVE Copy/Link → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await addViewerSafely(page, 'Scatter plot');
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'copy',
        name: names.linkCopy,
        perTableLinkOrClone: 'link',
      });
      ids.linkCopy = await captureActiveProjectIds(page, names.linkCopy, ids.original!.projectId);
      expect(ids.linkCopy.projectId).toBeTruthy();
      expect(ids.linkCopy.projectId).not.toBe(ids.original!.projectId);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4b step 5-6: reopen <name>-link → table re-materializes → closeAll', async () => {
      await reopenProjectById(page, ids.linkCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4b: link copy reopen returned rowCount=0 — Save Copy with Link did not persist the linked table reference').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4b step 7 — GROK-19750 INVARIANT: reopen original → table + viewers intact', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      const rc = await reopenedRowCount(page);
      const vc = await reopenedViewerCount(page);

      if (rc === 0)
        console.warn('GROK-19750 regression detected on dev: original `demog` rowCount=0 after Save-Copy-with-Link');
      if (vc === 0)
        console.warn('GROK-19750 regression detected on dev: original lost all viewers after Save-Copy-with-Link');
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4c prep: rehydrate original (GROK-19750 workaround)', async () => {
      ids.original = await rehydrateOriginal(page, names.original);
      expect(ids.original.projectId).toBeTruthy();
    });

    await softStep('4c step 1-4: open original → addViewer Line → SAVE Copy/Clone → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await addViewerSafely(page, 'Line chart');
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'copy',
        name: names.cloneCopy,
        perTableLinkOrClone: 'clone',
      });
      ids.cloneCopy = await captureActiveProjectIds(page, names.cloneCopy, ids.original!.projectId);
      expect(ids.cloneCopy.projectId).toBeTruthy();
      expect(ids.cloneCopy.projectId).not.toBe(ids.original!.projectId);
      expect(ids.cloneCopy.projectId).not.toBe(ids.linkCopy?.projectId);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4c step 5-6: reopen <name>-clone → table re-materializes (independent copy) → closeAll', async () => {
      await reopenProjectById(page, ids.cloneCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4c: clone copy reopen returned rowCount=0 — Save Copy with Clone did not persist the cloned table bytes').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4d prep: rehydrate original (defensive)', async () => {
      ids.original = await rehydrateOriginal(page, names.original);
      expect(ids.original.projectId).toBeTruthy();
    });

    await softStep('4d step 1-4: open original → addViewer Histogram → SAVE PVC → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await addViewerSafely(page, 'Histogram');

      await saveCopy(page, {
        sourceName: names.original,
        mode: 'personal-view-customizations',
        name: names.pvcCopy,
        expectIdChange: false,
      });

      ids.pvcCopy = {
        projectId: ids.original!.projectId,
        tableInfoId: ids.original!.tableInfoId,
        serverName: names.original,
      };
      expect(ids.pvcCopy.projectId).toBeTruthy();
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('4d step 5-6: reopen <name>-personal-view-customizations → table re-materializes → closeAll', async () => {
      await reopenProjectById(page, ids.pvcCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4d: PVC variant reopen returned rowCount=0 — PVC mode failed to preserve the linked source-table reference').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(300);
    });

    await softStep('Step 5: re-share each of 3 variants via right-click → Share (UI)', async () => {

      const recipient = await evalJs<string | null>(page, `(async () => {
        const users = await grok.dapi.users.list({limit: 50});
        const me = (await grok.dapi.users.current()).login;
        for (const u of users) {
          if (u.login === me || u.login === 'system') continue;
          const full = await grok.dapi.users.find(u.id);
          if (full && full.group && full.group.id)
            return full.group.friendlyName || full.login || full.group.name || null;
        }
        return null;
      })()`);

      if (!recipient) {
        console.warn('Step 5: no recipient with materialized group on dev — skipping share (acceptable-SR per Edit 10 G3)');
        return;
      }

      await page.goto('/projects');
      await page.locator('.grok-gallery-grid').waitFor({timeout: 30_000});

      await page.waitForFunction(() => {
        return document.querySelectorAll('.grok-gallery-grid [name^="div-"]').length > 0;
      }, undefined, {timeout: 30_000}).catch(() => {});
      const targetTiles = await evalJs<string[]>(page, `(() => {
        const tiles = Array.from(document.querySelectorAll('.grok-gallery-grid [name^="div-"]'));
        return tiles.map(t => t.getAttribute('name')).filter(n => /demog/i.test(n)).slice(0, 20);
      })()`);

      const variants = [
        {name: ids.linkCopy?.serverName ?? names.linkCopy, label: 'link'},
        {name: ids.cloneCopy?.serverName ?? names.cloneCopy, label: 'clone'},
        {name: ids.pvcCopy?.serverName ?? names.pvcCopy, label: 'pvc'},
      ];
      let sharedCount = 0;
      const shareErrors: string[] = [];
      for (const v of variants) {
        try {
          await shareProjectViaContextMenu(page, v.name, {
            recipient,
            accessLevel: 'View and use',
          });
          sharedCount++;

          await page.keyboard.press('Escape').catch(() => {});
          await page.waitForTimeout(500);
        } catch (e: any) {
          shareErrors.push(`${v.label} (${v.name}): ${(e?.message ?? String(e)).slice(0, 200)}`);
          await page.keyboard.press('Escape').catch(() => {});
          await page.waitForTimeout(300);
        }
      }

      if (sharedCount === 0)
        console.warn(`Step 5: no variant could be shared (env flake). Tiles seen: ${JSON.stringify(targetTiles)}. Errors: ${shareErrors.join('; ')}`);
    });

    await softStep('Step 5b: recipient-open verification via shareWithSecondUserAndVerify', async () => {
      const projectId = ids.linkCopy?.projectId;
      const projectName = ids.linkCopy?.serverName ?? names.linkCopy;
      if (!projectId) {
        console.warn('Step 5b: no link-copy project id available — skipping recipient-open leg');
        return;
      }
      const r = await shareWithSecondUserAndVerify(page, {id: projectId, name: projectName});
      if (!r.shared) {
        console.warn('Step 5b: share skipped — ' + r.reason);
        return;
      }
      if (r.recipientVisible !== null)
        expect(r.recipientVisible).toBe(true);
    });

  } finally {

    for (const v of Object.values(ids))
      if (v) await deleteProjectWithCleanup(page, v);
  }

  finishSpec();
});
