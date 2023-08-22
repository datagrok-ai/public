import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DockerContainer, DockerContainersDataSource} from 'datagrok-api/dg';

import {delay, awaitCheck, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';


export async function awaitContainerStart(ms: number = 1000): Promise<DockerContainer> {
    const dockerContainer = await grok.dapi.docker.dockerContainers.filter(_package.name).first();
    if (dockerContainer.status !== 'started' && dockerContainer.status !== 'checking')
        await delay(ms);
    return dockerContainer;
}

category('Dapi: Docker', () => {

    test('Get response', async () => {
        const dockerContainer = await awaitContainerStart();
        console.log(`Docker container ID: ${dockerContainer.id}`);
        console.log(`Docker container status: ${dockerContainer.status}`);
        expect(dockerContainer.status === 'started' || dockerContainer.status === 'checking',
            true, 'Docker container has not started yet');

        const params: RequestInit = {
            method: 'GET'
        };
        const path = '/square?number=4';
        const response = await grok.dapi.docker.dockerContainers.request(dockerContainer.id, path, params);
        expect(response?.trim(), '{"result":16}');
    });
});
