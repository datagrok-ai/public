import * as grok from 'datagrok-api/grok';

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
