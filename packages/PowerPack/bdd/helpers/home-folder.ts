/* Shared by the bindings that put a file into the account's own "My files" share, and outside
   bindings/ because a binding module imported by another is loaded twice and its steps then resolve
   to no export. */
import type {Page} from '@playwright/test';
import {atFeatureEnd, expect, fixtureFamilies, isStaleFixture, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** Writes the file into the "My files" share of the signed-in account and deletes it when the feature
 * ends, read back; `write` gets the file's path ("<project>:Home/<name>") and writes it in the page. A file
 * of the same {run} or {time} family over an hour old is a killed run's, and goes first. */
export async function putIntoHomeFolder(page: Page, name: string,
  write: (path: string) => Promise<void>): Promise<void> {
  const home: string = await page.evaluate(async () => {
    const project = grok.shell.user.project.name;
    const share = (await grok.dapi.connections.list())
      .find((c: any) => c.dataSource === 'Files' && c.nqName === `${project}:Home`);
    if (!share)
      throw new Error(`the account has no home folder ${project}:Home on this stand`);
    return String(share.nqName);
  });
  const extension = /\.[^.]+$/.exec(name)?.[0] ?? '';
  const base = (file: string): string =>
    file.endsWith(extension) ? file.slice(0, file.length - extension.length) : file;
  const families = fixtureFamilies([base(name)]);
  if (families.length > 0) {
    const files: {name: string; changed: number}[] = await page.evaluate(async (dir) =>
      (await grok.dapi.files.list(`${dir}/`, false)).filter((f: any) => !f.isDirectory)
        .map((f: any) => ({name: String(f.name), changed: f.updatedOn?.valueOf() ?? 0})), home);
    const stale = files.filter((f) => f.name.endsWith(extension) &&
      isStaleFixture({name: base(f.name), friendlyName: base(f.name), createdOn: f.changed}, families));
    await page.evaluate(async ([dir, names]) => {
      for (const n of names)
        await grok.dapi.files.delete(`${dir}/${n}`);
    }, [home, stale.map((f) => f.name)] as [string, string[]]);
  }
  const path = `${home}/${name}`;
  atFeatureEnd(page, async () => {
    await page.evaluate(async (p) => {
      if (await grok.dapi.files.exists(p))
        await grok.dapi.files.delete(p);
    }, path);
    await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
      {message: `${path} still in My files`, timeout: pollMs(15000)}).toBe(false);
  });
  await write(path);
}
