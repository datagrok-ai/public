import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';
import { testEvent } from './test';

// TODO: Use type from datagrok-api
export type DockerContainerStatus = 'stopped' | 'started' | 'pending change' | 'changing' | 'error' | 'checking';

export async function awaitStatus(
  dcId: string, targetStatus: DockerContainerStatus, timeout: number = 30000, logger: DG.PackageLogger
): Promise<void> {
  const t1: number = window.performance.now();
  const event = new Subject<DG.DockerContainer>();
  const pollingHandler = async (): Promise<DG.DockerContainer> => {
    const dc = await grok.dapi.docker.dockerContainers.find(dcId);
    if (dc.status === targetStatus)
      event.next(dc);
    return dc;
  };
  let interval!: number;
  try {
    await testEvent(event, (dc: DG.DockerContainer) => {
      const t2: number = window.performance.now();
      logger.debug(`awaitStatus(), ` +
        `docker container ('${dc.name}') GET status = '${targetStatus}' in ${t2 - t1} ms.`);
    }, async () => {
      const dc = await pollingHandler(); // -> event
      logger.debug(`awaitStatus(), ` +
        `docker container ('${dc.name}') HAS status = '${dc.status}'.`);
      interval = window.setInterval(pollingHandler, 200);
    }, timeout);
  } finally {
    window.clearInterval(interval);
  }
}
