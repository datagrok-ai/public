/* Custom platform events (`grok.events.fireCustomEvent` / `onCustomEvent`) — the word a package
   gives about work it finished off-screen (Bio's `bio-monomer-lib-loaded`): listened for by id
   from a step, read by the one that claims it. */
import {type Page} from '@playwright/test';
import {expect} from './patience.js';
import {installViewerRuntime} from './viewers.js';

interface CustomRead {
  count: number;
  last: unknown;
}

export async function listenCustomEvent(page: Page, id: string): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate((i) => { (window as any).__bdd.listenCustom(i); }, id);
}

/** The event fired at least once since "listens for" or the previous read — with `accept`, a firing
 * whose arguments it accepts, so a late one left over from before is waited past; the read zeroes the
 * count and returns the last event's arguments. */
export async function expectCustomEvent(page: Page, id: string, timeoutMs = 30000,
  accept?: (args: unknown) => boolean): Promise<unknown> {
  await installViewerRuntime(page);
  const read = (take: boolean): Promise<CustomRead> => page.evaluate(([i, t]) => (window as any).__bdd.customFired(i, t), [id, take] as [string, boolean]);
  let last: CustomRead = await read(false);
  if (last.count < 0)
    throw new Error(`the "${id}" custom event is not listened for in this scenario (Given user listens for "${id}" custom event)`);
  try {
    await expect.poll(async () => (last = await read(false)).count > 0 && (accept?.(last.last) ?? true), {timeout: timeoutMs}).toBe(true);
  }
  catch {
    throw new Error(last.count > 0
      ? `the "${id}" custom event fired, but not as expected: it last carried ${JSON.stringify(last.last)} (${Math.round(timeoutMs / 1000)} s)`
      : `the "${id}" custom event has not fired since it was listened for (${Math.round(timeoutMs / 1000)} s)`);
  }
  return (await read(true)).last;
}
