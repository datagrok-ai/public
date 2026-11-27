/* The Dashboards how-to clips (features/howto/06-dashboards-and-sharing): a layout the clip saves under a name a person
   would give it ("Height and weight"), which the library sweeps at feature start only for {time} or {run} names. */
import {type Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {deleteLayoutsAtEnd, escapeRegExp, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

export const ownLayoutsGone = Given('the user\'s own layouts named {string} are removed now and at feature end', async (page: Page, list: string) => {
  const names = list.split(',').map((n) => n.trim()).filter(Boolean);
  const pattern = `^(${names.map(escapeRegExp).join('|')})$`;
  await deleteLayoutsAtEnd(page, pattern, pattern);
}, {tier: 'api', description: 'the running account\'s layouts of those names (comma-separated): the ones it saves during the feature go when it ends, with their pictures and wrapper projects, and the ones an earlier run left (over an hour old) go now'});

/* The Save dialog's Presentation mode switch writes the project's ".presentation" meta param, which opening the
   project reads (shell_project.dart openProject). */
export const savedInPresentationMode = Then('the project named {string} should be saved to open in presentation mode', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(async (n) => {
    const listed = await grok.dapi.projects.filter(`friendlyName = ${JSON.stringify(n)} or name = ${JSON.stringify(n)}`).first();
    if (!listed)
      return `no project "${n}" on the server`;
    const project = await grok.dapi.projects.find(listed.id);
    return String(project.meta['.presentation'] ?? 'not set');
  }, name), {message: `the presentation-mode flag of the project "${name}" on the server`, timeout: pollMs(30000)}).toBe('true');
}, {tier: 'api', description: 'the ".presentation" meta param of the project as the server keeps it, which makes the project open in presentation mode'});
