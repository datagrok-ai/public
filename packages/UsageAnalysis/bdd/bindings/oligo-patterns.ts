/* The Oligo Pattern app (SequenceTranslator) keeps a saved pattern in the shared user settings
   `OligoToolkit`, one record per pattern hash: {patternConfig: {patternName, …}, authorID, date}. A
   feature that saves one removes it by name now and when the feature ends, and every earlier pattern of
   the same family (the name with another `-<time>` suffix) with it. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, fixtureFamilies, pollMs, RUN_SUFFIX} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const STORAGE = 'OligoToolkit';

type Saved = {hash: string; name: string; createdOn: number; mine: boolean; held: boolean};

/* The server's patterns and the page's own copy (`grok.userSettings`, read at sign-in and written back
   every 10 s), which the app reads: a pattern saved a moment ago is only in the page, one another run
   removed since the page signed in is only there too. */
function savedPatterns(page: Page): Promise<Saved[]> {
  return page.evaluate(async (storage) => {
    const me = (await grok.dapi.users.current()).id;
    const held = grok.userSettings.get(storage, false) ?? {};
    const saved = {...held, ...await grok.dapi.userDataStorage.get(storage, false) ?? {}};
    return Object.entries(saved).map(([hash, text]) => {
      let record: any = null;
      try {
        record = JSON.parse(text as string);
      }
      catch {
        // a record another tool wrote: not a pattern
      }
      return {hash, name: String(record?.patternConfig?.patternName ?? ''), createdOn: Date.parse(record?.date?.create ?? '') || 0,
        mine: record?.authorID === me, held: hash in held};
    });
  }, STORAGE);
}

/* Any age, not only a killed run's: the app keys a pattern by its content and refuses a second one with
   the same strands. A removal that read only the server missed the pattern the page had not written yet,
   which the page's next write then left there for the next run to be refused by. */
async function removePatterns(page: Page, name: string): Promise<void> {
  const families = fixtureFamilies([name]);
  const doomed = (all: Saved[]): Saved[] => all.filter((p) => p.mine &&
    (p.name === name || RUN_SUFFIX.test(p.name) && families.includes(p.name.replace(RUN_SUFFIX, ''))));
  const found = doomed(await savedPatterns(page));
  if (found.length > 0)
    await page.evaluate(async ([storage, list]) => {
      for (const {hash, held} of list) {
        if (held)
          grok.userSettings.delete(storage, hash, false);
        await grok.dapi.userDataStorage.remove(storage, hash, false);
      }
      await grok.userSettings.flush();
    }, [STORAGE, found] as [string, Saved[]]);
  await expect.poll(async () => doomed(await savedPatterns(page)).map((p) => p.name).join(', '),
    {message: `the account's oligo patterns named "${name}" (or of its family) on the server or in the page`,
      timeout: pollMs(30000)}).toBe('');
}

export const noOligoPattern = Given('no oligo pattern named {string} is in the user\'s settings', async (page: Page, name: string) => {
  atFeatureEnd(page, () => removePatterns(page, name));
  await removePatterns(page, name);
}, {tier: 'api', description: 'the account\'s saved Oligo Pattern patterns of that name, and of its family, go now and when the feature ends, from the server and from the page\'s own copy, read back from both'});
