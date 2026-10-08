/* The Oligo Pattern app (SequenceTranslator) keeps a saved pattern in the shared user settings
   `OligoToolkit`, one record per pattern hash: {patternConfig: {patternName, …}, authorID, date}. A
   feature that saves one removes it by name now and when the feature ends, and a killed run's
   pattern of the same family (the name with another `-<time>` suffix, over an hour old) with it. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, fixtureFamilies, isStaleFixture, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const STORAGE = 'OligoToolkit';

type Saved = {hash: string; name: string; createdOn: number; mine: boolean};

function savedPatterns(page: Page): Promise<Saved[]> {
  return page.evaluate(async (storage) => {
    const me = (await grok.dapi.users.current()).id;
    return Object.entries(await grok.dapi.userDataStorage.get(storage, false) ?? {}).map(([hash, text]) => {
      let record: any = null;
      try {
        record = JSON.parse(text as string);
      }
      catch {
        // a record another tool wrote: not a pattern
      }
      return {hash, name: String(record?.patternConfig?.patternName ?? ''), createdOn: Date.parse(record?.date?.create ?? '') || 0,
        mine: record?.authorID === me};
    });
  }, STORAGE);
}

async function removePatterns(page: Page, name: string): Promise<void> {
  const families = fixtureFamilies([name]);
  const doomed = (all: Saved[]): Saved[] => all.filter((p) => p.mine &&
    (p.name === name || isStaleFixture({name: p.name, friendlyName: p.name, createdOn: p.createdOn}, families)));
  const hashes = doomed(await savedPatterns(page)).map((p) => p.hash);
  if (hashes.length > 0)
    await page.evaluate(async ([storage, list]) => {
      for (const hash of list)
        grok.userSettings.delete(storage, hash, false);
      await grok.userSettings.flush();
    }, [STORAGE, hashes] as [string, string[]]);
  await expect.poll(async () => doomed(await savedPatterns(page)).map((p) => p.name).join(', '),
    {message: `the account's oligo patterns named "${name}" (or a killed run's of its family) in the shared user settings`,
      timeout: pollMs(30000)}).toBe('');
}

export const noOligoPattern = Given('no oligo pattern named {string} is in the user\'s settings', async (page: Page, name: string) => {
  atFeatureEnd(page, () => removePatterns(page, name));
  await removePatterns(page, name);
}, {tier: 'api', description: 'the account\'s saved Oligo Pattern patterns of that name go now and when the feature ends, read back from the server'});
