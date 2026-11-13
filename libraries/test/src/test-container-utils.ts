import * as grok from 'datagrok-api/grok';

import {tests} from './test';

/** Skips every test of [category] when the server reports Grok Spawner switched off or unreachable, so no container
 * of [what] can start; a test that already has a skip reason keeps it. Throws when the status cannot be read (no
 * `isActive` in it) or the platform has no `getStatus`, so the category fails instead. */
export async function skipCategoryWithoutSpawner(category: string, what: string): Promise<void> {
  if (typeof grok.dapi.docker.getStatus !== 'function')
    throw new Error('The server does not provide grok.dapi.docker.getStatus; update the platform');
  const status = await grok.dapi.docker.getStatus();
  if (status.isActive == null)
    throw new Error(`Cannot tell whether ${what} can run: the server's Grok Spawner status could not be read`);
  if (status.enabled && status.isActive)
    return;
  const reason = `${what} runs in a container, and Grok Spawner is ${status.enabled ? 'unreachable' : 'switched off'} on this server`;
  for (const t of tests[category].tests!)
    t.options!.skipReason ??= reason;
}

export async function ensureContainerRunning(containerName: string, timeout: number = 300000) {
  const container = await grok.dapi.docker.dockerContainers.filter(containerName).first();
  if (!(container.status.startsWith('started') || container.status.startsWith('checking'))) {
    console.log(`starting container ${container.name}`);
    await grok.dapi.docker.dockerContainers.run(container.id, false);
  }

  const deadline = Date.now() + timeout;
  let restartCount = 0;
  while (Date.now() < deadline) {
    const cont = await grok.dapi.docker.dockerContainers.find(container.id) as any;
    if (cont.status.startsWith('started') || cont.status.startsWith('checking'))
      return;
    if (cont.status === 'error' || cont.status.includes('stopped')) {
      if (restartCount < 3) {
        restartCount++;
        await grok.dapi.docker.dockerContainers.run(container.id, false);
      }
      else
        throw new Error(`Container ${containerName} failed: ${cont.status} after ${restartCount} restart attempts`);
    }
    await new Promise((r) => setTimeout(r, 5000));
  }
  throw new Error(`${containerName} hasn't been started after ${timeout / 60000} minutes`);
}
