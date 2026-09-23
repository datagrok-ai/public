/* Shared by the queries and the scripts bindings, and outside bindings/ because a binding module
   imported by another is loaded twice and its steps then resolve to no export. */
import type {Page} from '@playwright/test';
import {atFeatureEnd} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** A view that saves a layout leaves the layout and the project that wraps it on the server; the
 *  entity the layout belongs to does not take them with it. [pattern] matches the layout's name. */
export function deleteLayoutsAtEnd(page: Page, pattern: string): void {
  const since = Date.now() - 60 * 1000;
  atFeatureEnd(page, async () => {
    await page.evaluate(async ([from, source]) => {
      const me = (await grok.dapi.users.current()).id;
      const name = new RegExp(source as string);
      // the grok name drops what the friendly name keeps ("BDD-Q-layout-1" is "BDDQLayout1")
      const names = (x: any) => [String(x.friendlyName ?? ''), String(x.name ?? '')].filter(Boolean);
      const ours = (x: any) => String(x.author?.id ?? '') === String(me) &&
        (x.createdOn ? new Date(x.createdOn.toString()).getTime() : 0) >= (from as number);
      for (const layout of await grok.dapi.layouts.list({pageSize: 1000})) {
        if (!ours(layout) || !names(layout).some((n) => name.test(n)))
          continue;
        // the wrapper is found by name, so it is taken only when it is ours and of this run too:
        // "Df" is a name another account may hold on a shared stand
        const project = (await grok.dapi.projects.list({pageSize: 1000}))
          .find((p: any) => ours(p) && names(p).some((n) => names(layout).includes(n)));
        if (project)
          await grok.dapi.projects.delete(project);
        await grok.dapi.layouts.delete(layout);
      }
    }, [since, pattern] as [number, string]);
  });
}
