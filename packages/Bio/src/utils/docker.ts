import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';

import {Pepsea} from './pepsea';

/** Waits if container is not started
 * @param {number} ms - time to wait in milliseconds */
export async function awaitContainerStart(ms: number = 60000): Promise<void> {
  const dc = await Pepsea.getDockerContainer();
  await startDockerContainer(dc, ms);
}

export async function startDockerContainer(argDc: DG.DockerContainer, timeout: number = 60000): Promise<void> {
  // argDc contains current container status
  let dc: DG.DockerContainer | null = argDc;
  let end: boolean = false;
  for (let i = 0; i < timeout / 200; ++i) {
    if (dc === null) dc = await grok.dapi.docker.dockerContainers.find(argDc.id);

    if (isStarted(dc)) {
      end = true;
      break;
    }

    switch (dc.status) {
    case 'stopped': {
      // TODO: Use the new dockerContainers API
      await grok.dapi.docker.dockerContainers.run(dc.id);
      break;
    }
    case 'pending change':
    case 'changing': {
      // skip to wait
      break;
    }
    case 'error': {
      throw new Error('Docker container error state.');
    }
    }
    dc = null;
    await delay(200);
  }
  if (!end) throw new Error('Docker container start timeout.');
  // this.dc = await grok.dapi.docker.dockerContainers.find(dcId);
}

export function isStarted(dc: DG.DockerContainer): boolean {
  return dc.status === 'checking' || dc.status === 'started';
}
