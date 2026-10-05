// Lists, or with --delete deletes, the projects with a given name that the current user saved during a recording.
// Usage: node dd-projects.mjs "<project name>" [--delete]
import {open} from '../rec-lib.mjs';

const name = process.argv[2];
const del = process.argv.includes('--delete');
if (!name) throw new Error('usage: node dd-projects.mjs "<project name>" [--delete]');

const {browser, page} = await open();
try {
  const found = await page.evaluate(async ({name, del}) => {
    const me = grok.shell.user.id;
    const all = await grok.dapi.projects.filter(`friendlyName = "${name}"`).include('author').list();
    const mine = all.filter((p) => p.friendlyName === name && p.author?.id === me);
    if (del)
      for (const p of mine) await grok.dapi.projects.delete(p);
    return mine.map((p) => `${p.id} ${p.nqName}`);
  }, {name, del});
  console.log(`${del ? 'deleted' : 'found'}: ${found.length ? found.join(', ') : 'none'}`);
}
finally { await browser.close(); }
