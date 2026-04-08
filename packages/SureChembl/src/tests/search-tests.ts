import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/test/src/test';

export const CONTAINER_TIMEOUT = 1200000;

async function ensureContainerRunningDebug(containerName: string, timeout: number) {
  console.log(`[docker-debug] Looking for container '${containerName}'...`);
  const containers = await grok.dapi.docker.dockerContainers.filter(containerName).list();
  console.log(`[docker-debug] Found ${containers.length} container(s): ${containers.map((c) => `${c.name}(${c.status})`).join(', ')}`);

  const container = containers[0] as any;
  if (!container)
    throw new Error(`Container '${containerName}' not found`);

  console.log(`[docker-debug] Container id=${container.id} name=${container.name} status=${container.status} image=${container.dockerImage} address=${container.address}`);

  const images = await grok.dapi.docker.dockerImages.filter(containerName).list();
  console.log(`[docker-debug] Docker images matching '${containerName}': ${images.map((i) => `${i.name}(status=${(i as any).status}, version=${(i as any).version})`).join(', ')}`);

  if (!(container.status.startsWith('started') || container.status.startsWith('checking'))) {
    console.log(`[docker-debug] Starting container ${container.name} (current status: ${container.status})...`);
    await grok.dapi.docker.dockerContainers.run(container.id, false);
    console.log(`[docker-debug] run() returned, polling for status...`);
  }

  const deadline = Date.now() + timeout;
  let pollCount = 0;
  while (Date.now() < deadline) {
    const cont = await grok.dapi.docker.dockerContainers.find(container.id) as any;
    pollCount++;
    if (pollCount <= 10 || pollCount % 10 === 0)
      console.log(`[docker-debug] Poll #${pollCount}: status=${cont.status} address=${cont.address}`);
    if (cont.status.startsWith('started') || cont.status.startsWith('checking')) {
      console.log(`[docker-debug] Container ready after ${pollCount} polls`);
      return;
    }
    if (cont.status === 'error') {
      const logs = await grok.dapi.docker.getServiceLogs(cont.serviceName, 50).catch(() => 'no logs');
      console.log(`[docker-debug] Container error. Logs: ${logs}`);
      throw new Error(`Container ${containerName} failed with error. Logs: ${logs}`);
    }
    await new Promise((r) => setTimeout(r, 5000));
  }
  throw new Error(`${containerName} hasn't been started after ${timeout / 60000} minutes (${pollCount} polls)`);
}

category('Search tests', () => {
  test('Substructure search', async () => {
    await ensureContainerRunningDebug('surechembl', CONTAINER_TIMEOUT);
    const df: DG.DataFrame | null = await grok.functions.call('Surechembl:sureChemblSubstructureSearch', {
      molecule: 'FC(F)(F)c1ccc(OC2CCNCC2)cc1',
      limit: 10,
    });
    expect(df && df.rowCount > 0, true);
    const patentIdx = df!.col('doc_id')?.toList().findIndex((it) => it === 'US-11918575-B2');
    expect(patentIdx !== undefined && patentIdx !== -1, true);
    expect(df!.get('smiles', patentIdx!), 'FC(F)(F)c1ccc(OC2CCNCC2)cc1');
    expect(df!.get('doc_surechembl_id', patentIdx!), 81759);

    if (df != null)
      grok.shell.closeTable(df);
  }, {timeout: 30000 + CONTAINER_TIMEOUT});

  test('Similarity search', async () => {
    await ensureContainerRunningDebug('surechembl', CONTAINER_TIMEOUT);
    const df: DG.DataFrame | null = await grok.functions.call('Surechembl:sureChemblSimilaritySearch', {
      molecule: 'FC(F)(F)c1ccc(OC2CCNCC2)cc1',
      limit: 10,
      similarityThreshold: 0.6,
    });
    expect(df && df.rowCount > 0, true);
    const patentIdx = df!.col('doc_id')?.toList().findIndex((it) => it === 'US-11918575-B2');
    expect(patentIdx !== undefined && patentIdx !== -1, true);
    expect(df!.get('smiles', patentIdx!), 'FC(F)(F)c1ccc(OC2CCNCC2)cc1');
    expect(df!.get('doc_surechembl_id', patentIdx!), 81759);

    if (df != null)
      grok.shell.closeTable(df);
  }, {timeout: 30000 + CONTAINER_TIMEOUT});
});
