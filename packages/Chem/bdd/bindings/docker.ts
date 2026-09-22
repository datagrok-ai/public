/* The Chem service the Descriptors dialog reads its tree from. Its state is a precondition of the
   scenarios that watch the dialog, not their subject, so it is started through the platform's own
   docker API when it is not running, and left running: other features and users rely on it. */
import {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

export const containerRunning = Given('the {string} container is running', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const containers = await grok.dapi.docker.dockerContainers.filter(`name = "${n}"`).list();
    if (containers.length === 0)
      throw new Error(`no docker container named "${n}" on the server`);
    await grok.dapi.docker.dockerContainers.run(containers[0].id, true);
  }, name);
  // the server answers the run request before the container reports itself started ("checking")
  await expect.poll(() => page.evaluate(async (n) => {
    const containers = await grok.dapi.docker.dockerContainers.filter(`name = "${n}"`).list();
    return String(containers[0].status);
  }, name), {message: `the status of the "${name}" container`}).toMatch(/started|running/i);
}, {tier: 'api', description: 'starts the docker container of a package and waits for it to report itself started, so a dialog that reads from it has something to read'});
