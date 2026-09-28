/* Shared by the queries and the scripts bindings, and outside bindings/ because a binding module
   imported by another is loaded twice and its steps then resolve to no export. */
import type {Page} from '@playwright/test';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const HOUR_MS = 60 * 60 * 1000;

/** Deletes this account's layouts whose name matches [pattern] and that were made since [from], and
 * those matching [stale] made over an hour ago (a killed run's), each with the project that wraps
 * it; then reads the listing back until none of them is left. */
async function deleteLayouts(page: Page, pattern: string, from: number, stale: string | null): Promise<void> {
  // no named function inside the page function: tsx wraps one in a `__name` call the page lacks
  const sweep = (remove: boolean): Promise<string[]> => page.evaluate(async ([source, since, old, cutoff, del]) => {
    const me = String((await grok.dapi.users.current()).id);
    const name = new RegExp(source as string);
    const family = old ? new RegExp(old as string) : null;
    const listed: [string, any][] = [
      ...(await grok.dapi.layouts.list({pageSize: 1000})).map((x: any) => ['layout', x]),
      ...(del ? (await grok.dapi.projects.list({pageSize: 1000})).map((x: any) => ['project', x]) : []),
    ];
    const doomed: {kind: string; entity: any; names: string[]}[] = [];
    for (const [kind, x] of listed) {
      // the grok name drops what the friendly name keeps ("BDD-Q-layout-1" is "BDDQLayout1")
      const names = [String(x.friendlyName ?? ''), String(x.name ?? '')].filter(Boolean);
      const at = x.createdOn ? new Date(x.createdOn.toString()).getTime() : 0;
      if (String(x.author?.id ?? '') === me && ((at >= (since as number) && names.some((n) => name.test(n))) ||
        (family != null && at < (cutoff as number) && names.some((n) => family.test(n)))))
        doomed.push({kind, entity: x, names});
    }
    const layouts = doomed.filter((d) => d.kind === 'layout');
    for (const layout of del ? layouts : []) {
      // the wrapper is found by name, so it is taken only when it is ours and doomed too: "Df" is a
      // name another account may hold on a shared stand
      const project = doomed.find((d) => d.kind === 'project' && d.names.some((n) => layout.names.includes(n)));
      if (project)
        await grok.dapi.projects.delete(project.entity);
      await grok.dapi.layouts.delete(layout.entity);
    }
    return layouts.map((d) => d.names[0]);
  }, [pattern, from, stale, Date.now() - HOUR_MS, remove] as [string, number, string | null, number, boolean]);
  if ((await sweep(true)).length > 0)
    await expect.poll(() => sweep(false), {message: `layouts matching ${pattern} still on the server`, timeout: pollMs(15000)}).toEqual([]);
}

/** A view that saves a layout leaves the layout and the project that wraps it on the server; the
 *  entity the layout belongs to does not take them with it. [pattern] matches the layout's name. */
export function deleteLayoutsAtEnd(page: Page, pattern: string): void {
  const since = Date.now() - 60 * 1000;
  atFeatureEnd(page, () => deleteLayouts(page, pattern, since, null));
}

/** The same for the layout named after a feature's entity, swept now as well: the exact name at any
 * age, and the older layouts of its family (the name with another run's time) a killed run left. */
export async function deleteNamedLayouts(page: Page, name: string): Promise<void> {
  const escape = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const exact = `^${escape(name)}$`;
  const prefix = /^(.*\D)\d{13,}$/.exec(name)?.[1];
  await deleteLayouts(page, exact, 0, prefix ? `^${escape(prefix)}\\d{13,}$` : null);
  deleteLayoutsAtEnd(page, exact);
}
