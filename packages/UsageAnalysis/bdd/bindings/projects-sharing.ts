/* What only the Projects features need of a project's own identity: what the Copy submenu of a
   project card (or Details > Links...) puts on the clipboard, and a project saved over by another
   account. */
import {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {expect, gestures, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* --- the Copy submenu of a project ---------------------------------------------------------------
   What the menu should copy is rebuilt from the server's entity, field by field: the id; the grok
   name (<namespace>:<name>, the namespace being the owner's — "Admin" on a stand, capitalized);
   the markup #{x.<grok name>."<friendly name>"}; the link <server>/p/<namespace>.<name>. */
const COPY_FORMS = ['id', 'grok name', 'markup', 'url'] as const;

export const clipboardHoldsProjectForm = Then('the clipboard should hold the {word} of the {string} project', async (page: Page, form: string, name: string) => {
  const wanted = form.toLowerCase().replace(/-/g, ' ');
  if (!(COPY_FORMS as readonly string[]).includes(wanted))
    throw new Error(`no "${form}" of a project; one of: ${COPY_FORMS.map((f) => f.replace(' ', '-')).join(', ')}`);
  const expected: string = await page.evaluate(async ([n, f]) => {
    const listed = await grok.dapi.projects.filter(`name = "${n}"`).first();
    if (!listed)
      throw new Error(`no project "${n}" on the server`);
    const nq = String(listed.nqName);
    switch (f) {
      case 'id': return String(listed.id);
      case 'grok name': return nq;
      case 'markup': return `#{x.${nq}."${listed.friendlyName}"}`;
      default: return `${location.origin}/p/${nq.replace(':', '.')}`;
    }
  }, [name, wanted] as [string, string]);
  await expect.poll(() => gestures.readClipboard(page), {message: `the clipboard, against the ${wanted} of "${name}" (${expected})`}).toBe(expected);
}, {description: 'id, grok-name, markup or URL — the expected text built from the project the server holds under that name'});

/* --- a project saved over ---------------------------------------------------------------------------
   "Save original project" writes the same project again: it keeps its id and its update time moves
   on. A copy would be a new project under the same name, which leaves the original's time as it was. */
const savedTimes = new Map<string, {id: string; at: number}>();

async function projectSaved(page: Page, name: string): Promise<{id: string; at: number} | string> {
  return page.evaluate(async (n) => {
    const listed = await grok.dapi.projects.filter(`friendlyName = "${n}" or name = "${n}"`).list();
    if (listed.length !== 1)
      return `${listed.length} projects named "${n}" on the server`;
    const updated = (await grok.dapi.projects.find(listed[0].id)).updatedOn;
    return {id: String(listed[0].id), at: updated?.valueOf() ?? 0};
  }, name);
}

export const rememberProjectSaved = When('user remembers when the {string} project was saved', async (page: Page, name: string) => {
  const got = await projectSaved(page, name);
  if (typeof got === 'string')
    throw new Error(got);
  savedTimes.set(name, got);
}, {tier: 'api', description: 'the id and update time of the one project of that name on the server'});

export const projectSavedAgain = Then('the {string} project should have been saved again since remembered', async (page: Page, name: string) => {
  const before = savedTimes.get(name);
  if (!before)
    throw new Error(`no save time remembered for "${name}": "user remembers when the … project was saved" first`);
  await expect.poll(async () => {
    const now = await projectSaved(page, name);
    return typeof now === 'string' ? now : now.id !== before.id ? `another project (${now.id}), not ${before.id}` :
      now.at > before.at ? 'saved again' : `not saved since ${new Date(before.at).toISOString()}`;
  }, {message: `the "${name}" project against its remembered save`, timeout: pollMs(30000)}).toBe('saved again');
}, {tier: 'api', description: 'the same project (same id) with a later update time than remembered'});
