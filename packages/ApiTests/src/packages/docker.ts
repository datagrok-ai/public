import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import {DockerContainer} from 'datagrok-api/dg';
import {delay, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';


export async function awaitContainerStart(): Promise<DockerContainer> {
  let dockerContainer;
  do {
    dockerContainer = await grok.dapi.docker.dockerContainers.filter(_package.name).first();
    await delay(500);
  }
  while (!dockerContainer || (dockerContainer?.status !== 'started' && dockerContainer?.status !== 'checking'));
  return dockerContainer;
}

category('Packages: Docker', () => {
  test('Get response', async () => {
    const dockerContainer = await awaitContainerStart();
    console.log(`Docker container ID: ${dockerContainer.id}`);
    console.log(`Docker container status: ${dockerContainer.status}`);
    expect(dockerContainer.status === 'started' || dockerContainer.status === 'checking',
      true, 'Docker container has not started yet');

    const params: RequestInit = {
      method: 'GET',
    };
    const path = '/square?number=4';
    const response = await grok.dapi.docker.dockerContainers.request(dockerContainer.id, path, params);
    expect(response?.trim(), '{"result":16}');
  });
});
