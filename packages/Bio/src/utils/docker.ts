import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {delay} from '@datagrok-libraries/utils/src/test';

export async function startDockerContainer(dcId: string, timeout: number = 30000): Promise<void> {
  // TODO: Use the new dockerContainers API
  const res = await grok.dapi.docker.dockerContainers.run(dcId /*, true */);
  let end: boolean = false;
  for (let i = 0; i < timeout / 200; ++i) {
    const dc = await grok.dapi.docker.dockerContainers.find(dcId);
    switch (dc.status) {
      case 'stopped': {
        await grok.dapi.docker.dockerContainers.run(dcId);
        break;
      }
      case 'pending change':
      case 'changing': {
        // skip to wait
        break;
      }
      case 'checking':
      case 'started': {
        end = true;
        break;
      }
      case 'error': {
        throw new Error('Docker container error state.');
      }
    }
    if (end) break;
    await delay(200);
  }
  if (!end) throw new Error('Docker container run timeout.');
  // this.dc = await grok.dapi.docker.dockerContainers.find(dcId);
}
