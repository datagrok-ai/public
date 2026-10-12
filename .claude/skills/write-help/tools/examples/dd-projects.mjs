// Lists the current user's projects with a given name, or deletes one project by the id a recording printed.
// Usage: node dd-projects.mjs "<project name>"     node dd-projects.mjs --delete <id>
// Deletes only by id, and only a project of the current user: another user's project with the same name, or an
// older one of yours, is never touched.
import {open} from '../rec-lib.mjs';

const del = process.argv.indexOf('--delete');
const id = del >= 0 ? process.argv[del + 1] : null;
const name = del >= 0 ? null : process.argv[2];
if (!id && !name) throw new Error('usage: node dd-projects.mjs "<project name>" | --delete <id>');

const {browser, page} = await open();
try {
  const result = await page.evaluate(async ({id, name}) => {
    const me = grok.shell.user.id;
    if (id) {
      const p = await grok.dapi.projects.include('author').find(id);
      if (!p) return `no project ${id}`;
      if (p.author?.id !== me) return `not deleted: ${id} "${p.friendlyName}" belongs to ${p.author?.login}`;
      await grok.dapi.projects.delete(p);
      return (await grok.dapi.projects.find(id).catch(() => null)) ? `still there: ${id}` : `deleted: ${id} "${p.friendlyName}"`;
    }
    const all = await grok.dapi.projects.filter(`friendlyName = "${name}"`).include('author').list();
    const mine = all.filter((p) => p.friendlyName === name && p.author?.id === me);
    return `found: ${mine.length ? mine.map((p) => `${p.id} (created ${p.createdOn})`).join(', ') : 'none'}`;
  }, {id, name});
  console.log(result);
}
finally { await browser.close(); }
