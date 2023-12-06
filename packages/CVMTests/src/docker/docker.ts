import * as grok from 'datagrok-api/grok';
import {delay, category, expect, test} from '@datagrok-libraries/utils/src/test';

category('Packages: Docker', () => {
  const containerOnDemandName: string = 'Cvmtests-docker-test1';

  test('Get response: On demand', async () => {
    const container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    if (container.status !== 'stopped')
      await grok.dapi.docker.dockerContainers.stop(container.id, true);
    await testResponse(container.id);
  }, {timeout: 60000});

  test('Container timeout', async () => {
    let container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    await grok.dapi.docker.dockerContainers.run(container.id, true);
    await delay(80000);
    container = await grok.dapi.docker.dockerContainers.filter(containerOnDemandName).first();
    expect(container.status, 'stopped');
  }, {timeout: 120000});
});

async function testResponse(containerId: string): Promise<void> {
  const params: RequestInit = {
    method: 'GET',
  };
  const path = '/square?number=4';
  const response = await grok.dapi.docker.dockerContainers.request(containerId, path, params);
  expect(response?.trim(), '{"result":16}');
}
