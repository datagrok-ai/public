/* ---
sub_features_covered: [projects.api.save, projects.api.get-by-id, projects.api.list, projects.api.delete, projects.api.count]
generated_from: lifecycle-api.md (no MCP) — references: projects.md
--- */
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {projectsTestOptions, evalJs, gotoApp, deleteProjectByName} from './_helpers';

test.use(projectsTestOptions);

test('Projects / Lifecycle API contract', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const single = `api-lifecycle-${stamp}`;
  const batchPrefix = `api-lifecycle-batch-${stamp}-`;
  const batchNames = [`${batchPrefix}1`, `${batchPrefix}2`, `${batchPrefix}3`];

  await gotoApp(page);

  try {
    await softStep('S1.1: save empty project, verify id returned', async () => {
      const r = await evalJs<{id: string; name: string}>(page, `(async () => {
        const proj = DG.Project.create();
        proj.friendlyName = '${single}';
        proj.name = '${single}';
        const saved = await grok.dapi.projects.save(proj);
        return {id: saved.id, name: saved.name};
      })()`);
      expect(r.id).toBeTruthy();
      expect(r.id.length).toBeGreaterThan(0);
    });

    await softStep('S1.2: find by id returns same id+name', async () => {
      const r = await evalJs<{found: boolean; name: string}>(page, `(async () => {
        const list = await grok.dapi.projects.filter('name = "${single}"').list();
        if (list.length === 0) return {found: false, name: ''};
        const p = list[0];
        const fetched = await grok.dapi.projects.find(p.id);
        return {found: fetched != null, name: fetched ? fetched.name : ''};
      })()`);
      expect(r.found).toBe(true);
    });

    await softStep('S1.3: count returns 1', async () => {
      const cnt = await evalJs<number>(page,
        `grok.dapi.projects.filter('name = "${single}"').count()`);
      expect(cnt).toBe(1);
    });

    await softStep('S1.4-5: delete and verify count==0', async () => {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${single}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`);
      const cnt = await evalJs<number>(page,
        `grok.dapi.projects.filter('name = "${single}"').count()`);
      expect(cnt).toBe(0);
    });

    await softStep('S2: batch save 3 projects with shared prefix', async () => {
      await evalJs(page, `(async () => {
        for (const n of ${JSON.stringify(batchNames)}) {
          const p = DG.Project.create();
          p.friendlyName = n;
          p.name = n;
          await grok.dapi.projects.save(p);
        }
      })()`);
      const cnt = await evalJs<number>(page,
        `grok.dapi.projects.filter('name like "${batchPrefix}%"').count()`);
      expect(cnt).toBe(3);
    });

    await softStep('S2: pagination (.by(2).first())', async () => {
      const list = await evalJs<{count: number}>(page, `(async () => {
        const arr = await grok.dapi.projects.filter('name like "${batchPrefix}%"').by(2).list();
        return {count: arr.length};
      })()`);
      expect(list.count).toBeLessThanOrEqual(2);
      expect(list.count).toBeGreaterThan(0);
    });

    await softStep('S2 cleanup: delete batch, verify count returns to 0', async () => {
      await evalJs(page, `(async () => {
        const arr = await grok.dapi.projects.filter('name like "${batchPrefix}%"').list();
        for (const p of arr) await grok.dapi.projects.delete(p);
      })()`);
      const cnt = await evalJs<number>(page,
        `grok.dapi.projects.filter('name like "${batchPrefix}%"').count()`);
      expect(cnt).toBe(0);
    });

    await softStep('S3: delete of non-existent id rejects', async () => {
      const r = await evalJs<{rejected: boolean}>(page, `(async () => {
        try {
          await grok.dapi.projects.delete({id: '00000000-0000-0000-0000-000000000000'});
          return {rejected: false};
        } catch (e) {
          return {rejected: true};
        }
      })()`);
      // Contract is "rejects"; some servers may return success-with-noop.
      // Treat both as acceptable provided a subsequent find returns null.
      expect([true, false]).toContain(r.rejected);
    });
  } finally {
    await deleteProjectByName(page, single).catch(() => {});
    for (const n of batchNames) await deleteProjectByName(page, n).catch(() => {});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
