/* Custom platform events (`grok.events.fireCustomEvent`): the word a package gives about work it
   finished off-screen — Bio's `bio-monomer-lib-loaded` after the monomer libraries reload. A
   scenario listens by id, acts, and claims the event fired (or did not). */
import {Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Given, Then} from '../../src/registry.js';
import {expectCustomEvent, listenCustomEvent} from '../../src/runtime/events.js';

export const listenCustom = Given('user listens for {string} custom event', (page: Page, id: string) => listenCustomEvent(page, id),
  {tier: 'api', description: 'grok.events.onCustomEvent(id), counted until the page resets; a "should have fired" read zeroes the count'});

/* The task bar shows a progress entry for the length of a job and removes it when the job ends —
   often before any check could look. "watches" records every entry text the bar is given from then
   on, so a claim afterwards reads the record, not the bar. */
export const watchTaskBar = Given('user watches the task bar', (page: Page) => page.evaluate(() => {
  const w = window as any;
  w.__bddTaskBar?.observer.disconnect();
  const record: string[] = [];
  const observer = new MutationObserver((mutations) => {
    for (const m of mutations) {
      for (const n of Array.from(m.addedNodes)) {
        const el = n.nodeType === 1 ? n as Element : n.parentElement;
        const entry = el?.closest('.d4-task-bar-entry') ?? (el as Element | null)?.querySelector?.('.d4-task-bar-entry');
        if (entry)
          record.push((entry.textContent ?? '').trim());
      }
    }
  });
  observer.observe(document.body, {childList: true, subtree: true});
  w.__bddTaskBar = {observer, record};
}), {tier: 'ui', description: 'records the text of every task bar entry added from now until the page resets'});

export const taskBarFinished = Then('the task bar should have finished {string}', async (page: Page, text: string) => {
  await expect.poll(() => page.evaluate((t) => {
    const w = window as any;
    const shown = ((w.__bddTaskBar?.record ?? []) as string[]).some((s) => s.includes(t));
    const busy = Array.from(document.querySelectorAll('.d4-task-bar-entry')).some((e) => (e.textContent ?? '').includes(t));
    return !shown ? 'never shown' : busy ? 'still shown' : 'done';
  }, text), {message: `the task bar entry "${text}" since "user watches the task bar"`, timeout: pollMs(120000)}).toBe('done');
}, {description: 'an entry containing the text was shown since the watch began and is gone now: the job it stood for (a clustering, a computation) has ended, on a 120 s budget'});

export const taskBarShown = Then('the task bar should have shown {string}', async (page: Page, text: string) => {
  await expect.poll(() => page.evaluate(() => (window as any).__bddTaskBar?.record ?? null), {message: `task bar entries since "user watches the task bar"`})
    .toEqual(expect.arrayContaining([expect.stringContaining(text)]));
}, {description: 'an entry containing the text was shown at some point since the watch began, however briefly'});

export const customFired = Then('the {string} custom event should have fired', async (page: Page, id: string) => { await expectCustomEvent(page, id); },
  {description: 'at least once since "listens for" or the previous read (up to 30 s); reading zeroes the count'});
