import * as grok from 'datagrok-api/grok';
import { awaitCheck } from './test';

export async function ensureContainerRunning(containerName: string, timeout: number = 300000) {
  const container = await grok.dapi.docker.dockerContainers.filter(containerName).first();
  if (!(container.status.startsWith('started') || container.status.startsWith('checking'))) {
    console.log(`starting container ${container.name}`);
    await grok.dapi.docker.dockerContainers.run(container.id, false);
  }
  
  let started = false;
  await awaitCheck(() => {
    grok.dapi.docker.dockerContainers.find(container.id).then((cont) => {
      started = cont.status.startsWith('started') || cont.status.startsWith('checking');
    });
    return started;
  }, `${containerName} hasn't been started after ${timeout / 60000} minutes`, timeout, 5000);
}