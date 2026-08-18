import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';

test.use(specTestOptions);

const ROOT = 'System:AppData/UsageAnalysis';

test('Files cache — folder/file CRUD lifecycle', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  const stamp = Date.now();
  const dir = `${ROOT}/Folder cache test ${stamp}`;
  const dirRenamed = `${ROOT}/Folder cache test1 ${stamp}`;

  try {

    const created = await page.evaluate(async ({d}) => {
      const f = (window as any).grok.dapi.files;
      await f.writeAsText(`${d}/test.txt`, 'Hello world!');
      return {
        dirExists: await f.exists(d),
        fileExists: await f.exists(`${d}/test.txt`),
        content: await f.readAsText(`${d}/test.txt`),
      };
    }, {d: dir});
    expect(created.dirExists).toBe(true);
    expect(created.fileExists).toBe(true);
    expect(created.content).toBe('Hello world!');

    const fileRenamed = await page.evaluate(async ({d}) => {
      const f = (window as any).grok.dapi.files;
      await f.rename(`${d}/test.txt`, 'test1.txt');
      return {
        newExists: await f.exists(`${d}/test1.txt`),
        oldGone: !(await f.exists(`${d}/test.txt`)),
        names: (await f.list(d)).map((x: any) => x.name),
      };
    }, {d: dir});
    expect(fileRenamed.newExists).toBe(true);
    expect(fileRenamed.oldGone).toBe(true);
    expect(fileRenamed.names).toContain('test1.txt');

    const folderRenamed = await page.evaluate(async ({d, newName}) => {
      const f = (window as any).grok.dapi.files;
      await f.rename(d, newName);
      return d.replace(/[^/]+$/, newName);
    }, {d: dir, newName: `Folder cache test1 ${stamp}`});
    expect(folderRenamed).toBe(dirRenamed);

    const afterFolderRename = await page.evaluate(async ({oldDir, newDir}) => {
      const f = (window as any).grok.dapi.files;
      return {
        newDirExists: await f.exists(newDir),
        oldDirGone: !(await f.exists(oldDir)),
        fileMoved: await f.exists(`${newDir}/test1.txt`),
      };
    }, {oldDir: dir, newDir: dirRenamed});
    expect(afterFolderRename.newDirExists).toBe(true);
    expect(afterFolderRename.oldDirGone).toBe(true);
    expect(afterFolderRename.fileMoved).toBe(true);

    const deleted = await page.evaluate(async ({newDir}) => {
      const f = (window as any).grok.dapi.files;
      await f.delete(newDir);
      return {gone: !(await f.exists(newDir))};
    }, {newDir: dirRenamed});
    expect(deleted.gone).toBe(true);
  } finally {

    await softStep('Cleanup temp folders', async () => {
      await page.evaluate(async ({a, b}) => {
        const f = (window as any).grok.dapi.files;
        for (const p of [a, b]) { try { if (await f.exists(p)) await f.delete(p); } catch (_) {  } }
      }, {a: dir, b: dirRenamed});
    });
  }
});
