import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  delay,
  category,
  expect,
  test,
  expectObject,
  expectExceptionAsync,
  before,
} from '@datagrok-libraries/test/src/test';

category('Docker', () => {
  // Platform registers a package container as `kebab(package.name)-<dockerfileFolder>`
  // (see datlas deploy_docker.dart). Package `CvmTests` kebabs to `cvm-tests`, so the
  // dockerfiles/cvmtests-docker-test1|2 folders become these entity names.
  const containerOnDemandName: string = 'cvm-tests-cvmtests-docker-test1';
  const containerSimple: string = 'cvm-tests-cvmtests-docker-test2';
  const incorrectId: string = '00000000-0000-0000-0000-000000000000';

  // Queues the warm-up without awaiting anything: before() runs on the framework's fixed
  // 100s budget, and one awaited container start over it fails the whole category with
  // "before() failed" (build 1192). Every test starts what it needs under its own timeout.
  before(async () => {
    await grok.dapi.docker.dockerContainers.run((await findContainer(containerSimple)).id);
  });

  // Budget above datlas' own containerStatusTimeout (5 min) so a start that never
  // completes surfaces its real error instead of an opaque EXECUTION TIMEOUT.
  test('Get response: On demand', async () => {
    const container = await stopContainer(containerOnDemandName);
    await testResponse(container.id);
  }, {timeout: 360000 /*stressTest: true*/});

  test('Container timeout', async () => {
    let container = await stopContainer(containerOnDemandName);
    await grok.dapi.docker.dockerContainers.run(container.id, true);
    await delay(90000);
    container = await grok.dapi.docker.dockerContainers.filter(`name = "${containerOnDemandName}"`).first();
    expect(container.status.startsWith('stop'), true);
  }, {timeout: 240000, skipReason: 'Tool long'});

  test('Get response and logs: Incorrect', async () => {
    let response = await grok.dapi.docker.dockerContainers.fetchProxy(incorrectId, '/square?number=4');
    // container not found
    expect(response.status, 404, 'Status should be 404 indicating that container doesn\'t exist');
    const container = await startContainer(containerSimple);
    try {
      response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/foo?number=4');
      // response from docker with 404 - no such path mapping
      expect(response.status, 404, 'Status should be 404 indicating that in container there is not mapping for this url');
      response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/square?number=-0.01244');
      // response from docker
      expect(response.status, 500, 'Status should be 500 indicating that something went wrong in container');
      await grok.dapi.docker.dockerContainers.stop(container.id, true);
      response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/square?number=4');
      // response from server, container status is stopped
      expect(response.status, 400, 'Status should be 400 indicating that container is stopped and it is bad request');

      // requesting logs of container with incorrect id
      await expectExceptionAsync(async () => {
        await grok.dapi.docker.dockerContainers.getContainerLogs(incorrectId);
      });
    } finally {
      grok.dapi.docker.dockerContainers.run(container.id).then((_) => {}).catch((_) => {});
    }
  }, {timeout: 240000 /*stressTest: true*/});

  test('Proxy WebSocket', async () => {
    let ws: WebSocket | undefined;
    try {
      const container = await startContainer(containerSimple);
      await waitUntilAnswering(container.id);
      ws = await grok.dapi.docker.dockerContainers.webSocketProxy(container.id, '/ws');
      const testMessage = 'Hello World!';
      await new Promise<void>((res, rej) => {
        const onMessage = (event: MessageEvent) => {
          if (event.data === testMessage) {
            cleanup();
            res();
          } else {
            cleanup();
            rej(new Error(`First message wasn't the same as expected: ${event.data}`));
          }
        };

        const onError = (_: Event) => {
          cleanup();
          rej(new Error('WebSocket encountered an error'));
        };

        function cleanup() {
          ws!.removeEventListener('message', onMessage);
          ws!.removeEventListener('error', onError);
        }

        ws!.addEventListener('message', onMessage);
        ws!.addEventListener('error', onError);
        ws!.send(testMessage);
      });
    } finally {
      ws?.close();
    }
    // startContainer alone measured 27s on the CI stand — nothing left of the 30s default.
  }, {timeout: 240000}); // browser WebSocket global with cookie-session auth
});
// Tests are browser-only by default. This category has no node-capable tests yet:
// under the Node runtime the docker proxy returns 401 where the browser (cookie
// session) gets the real status (e.g. 404 for a missing container), and the
// WebSocket proxy relies on the browser WebSocket global.

async function findContainer(containerName: string): Promise<DG.DockerContainer> {
  const container = await grok.dapi.docker.dockerContainers.filter(`name = "${containerName}"`).first();
  if (!container)
    throw new Error(`Docker container "${containerName}" is not registered — the package's ` +
      `dockerfiles were not deployed to this server, or the name drifted from ` +
      `kebab(package.name)-<dockerfileFolder>.`);
  return container;
}

// run/stop reject with the bare response body, which reaches JS opaque (logged as just "null").
function describeFailure(action: string, container: DG.DockerContainer, name: string, e: any): Error {
  return new Error(`Failed to ${action} container "${name}" (id ${container.id}, ` +
    `status "${container.status}"): ${e?.message ?? e}`);
}

async function stopContainer(containerName: string): Promise<DG.DockerContainer> {
  const container = await findContainer(containerName);
  // DockerRouter.stop answers 400 for `error`; nothing to stop anyway, and a later run clears it.
  if (container.status === 'error')
    return container;
  //@ts-ignore
  if (!container.status.startsWith('stopped') && !(container.status.startsWith('pending') || container.status === 'stopping')) {
    try {
      await grok.dapi.docker.dockerContainers.stop(container.id, true);
    } catch (e: any) {
      throw describeFailure('stop', container, containerName, e);
    }
  }
  return container;
}

async function startContainer(containerName: string): Promise<DG.DockerContainer> {
  const container = await findContainer(containerName);
  try {
    await grok.dapi.docker.dockerContainers.run(container.id, true);
  } catch (e: any) {
    throw describeFailure('start', container, containerName, e);
  }
  return container;
}

/// `run(id, true)` returns once the platform marks the container started, which is earlier
/// than the app inside binding its port — opening the WebSocket proxy right then fails with
/// `Connection refused, errno = 111`. Poll a cheap HTTP path until it answers.
async function waitUntilAnswering(containerId: string, attempts: number = 45): Promise<void> {
  for (let i = 0; i < attempts; i++) {
    try {
      const response = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, '/square?number=4');
      if (response.status === 200)
        return;
    }
    catch (_) {}
    await delay(2000);
  }
  throw new Error(`Container ${containerId} did not answer within ${attempts * 2}s of being reported started`);
}

async function testResponse(containerId: string): Promise<void> {
  const path = '/square?number=4';
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, path);
  const body = await response.json();
  expect(response.status, 200, `Response status was ${response.status}, responded with: ${body['datagrok-error']}`);
  expectObject(body, {'result': 16});
}
