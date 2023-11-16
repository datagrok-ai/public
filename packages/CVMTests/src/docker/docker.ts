import * as grok from 'datagrok-api/grok';
import {delay, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {DockerContainer} from 'datagrok-api/dg';

category('Packages: Docker', () => {
  const containerOnDemandName: string = 'Cvmtests-docker-test1';

  test('Get response: On demand', async () => {
    const container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    if (container.status !== 'stopped') {
      await grok.dapi.docker.dockerContainers.stop(container.id);
      await delay(5000);
    }
    await testResponse(container.id);
  }, {timeout: 60000});

  test('Container timeout', async () => {
    let container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    await grok.dapi.docker.dockerContainers.run(container.id);
    await delay(20000);
    await awaitContainerStart(containerOnDemandName, 10);
    await delay(70000);
    container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    expect(container.status, 'stopped');
  }, {timeout: 120000});
});

// Duplicates, but I don't want to depend on other package with dockerfiles

export async function awaitContainerStart(containerName: string, tries: number = 3): Promise<DockerContainer> {
  let dockerContainer;
  let count: number = 0;
  do {
    dockerContainer = await grok.dapi.docker.dockerContainers.filter(containerName).first();
    count++;
    if (dockerContainer.status == 'started' || dockerContainer.status == 'checking')
      return dockerContainer;
    await delay(500);
  }
  while (count < tries);
  throw Error('Container didn\'t start');
}


async function testResponse(containerId: string): Promise<void> {
  const params: RequestInit = {
    method: 'GET',
  };
  const path = '/square?number=4';
  const response = await grok.dapi.docker.dockerContainers.request(containerId, path, params);
  expect(response?.trim(), '{"result":16}');
}
