import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** How `project` sits in the space named `space` (where the Save project dialog's Space row put it): 'owned', 'linked'
 * or 'absent'. */
function placement(page: Page, project: string, space: string): Promise<string> {
  return page.evaluate(async ([p, s]) => {
    const spaces: any[] = [...await grok.dapi.spaces.list()];
    let found: any = null;
    while (spaces.length > 0 && found === null) {
      const next = spaces.shift();
      if (next.friendlyName === s)
        found = next;
      else {
        const kids: any[] = await grok.dapi.spaces.id(next.id).children.ofTypes('Project').list();
        spaces.push(...kids.filter((c) => c.isSpace));
      }
    }
    if (found === null)
      return `no space named "${s}"`;
    const children = grok.dapi.spaces.id(found.id).children;
    const named = (list: any[]) => list.some((e) => e.friendlyName === p);
    if (named(await children.ofTypes('Project').list()))
      return 'owned';
    return named(await children.ofTypes('Project', true).list()) ? 'linked' : 'absent';
  }, [project, space] as [string, string]);
}

export const projectInSpace = Then('the project named {string} should be in the space {string} on the server',
  async (page: Page, project: string, space: string) => {
    await expect.poll(() => placement(page, project, space),
      {message: `where "${project}" sits in "${space}"`, timeout: pollMs(30000)}).toBe('owned');
  }, {tier: 'api', description: 'the space owns the project: it was saved or moved there, not linked'});

export const projectLinkedFromSpace = Then(
  'the project named {string} should be linked from the space {string} on the server',
  async (page: Page, project: string, space: string) => {
    await expect.poll(() => placement(page, project, space),
      {message: `where "${project}" sits in "${space}"`, timeout: pollMs(30000)}).toBe('linked');
  }, {tier: 'api', description: 'the space holds a link to the project, which lives elsewhere'});

export const projectNotInSpace = Then('the project named {string} should not be in the space {string} on the server',
  async (page: Page, project: string, space: string) => {
    await expect.poll(() => placement(page, project, space),
      {message: `where "${project}" sits in "${space}"`, timeout: pollMs(30000)}).toBe('absent');
  }, {tier: 'api', description: 'the space neither owns nor links the project'});
