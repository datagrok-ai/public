import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Predictive models: the engines packages register, and the /ml surface the Train Model view and the
// Apply dialog use. The JS API covers the entity (grok.dapi.models); creating one, its blobs and images,
// the zip and the suggestions have no wrapper, so they are called on the endpoints.

const ml = (path: string, init: RequestInit = {}) => fetch(`${grok.dapi.root}/ml${path}`, {credentials: 'include', ...init});

const blob = [1, 2, 3, 4, 5, 6, 7, 8];
const png = [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, 1, 2, 3];

async function bytesOf(response: Response): Promise<string> {
  return Array.from(new Uint8Array(await response.arrayBuffer())).join(',');
}

async function saveModel(): Promise<string> {
  const name = `apitests-model-${DG.Utils.randomString(8)}`;
  const response = await ml('', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({
    '#type': 'PredictiveModelInfo', name, friendlyName: name,
    input: [{name: 'AGE', type: 'int'}, {name: 'WEIGHT', type: 'double'}], output: [{name: 'outcome', type: 'double'}],
  })});
  expect(response.status, 200);
  return (await response.json()).id;
}

async function deleteModel(id: string): Promise<void> {
  const model = await grok.dapi.models.find(id).catch(() => null);
  if (model)
    await grok.dapi.models.delete(model);
}

function suggestedFor(columns: _DG.Column[]): Promise<any[]> {
  return ml('/suggestions', {method: 'POST', headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({'#type': 'TableInfo', columns: columns.map((c) => ({name: c.name, type: c.type}))})})
    .then((r) => r.json());
}

category('Dapi: models', () => {
  test('Engines: every train function has an apply function of the same mlname', async () => {
    const funcs = DG.Func.find({});
    const mlnames = (role: string) => funcs.filter((f) => f.options?.['mlrole'] === role).map((f) => f.options['mlname']);
    const train = mlnames('train');
    const apply = mlnames('apply');
    expect(train.length > 0, true, 'no function registers meta.mlrole: train');
    for (const name of train)
      expect(apply.includes(name), true, `the "${name}" engine trains but has no apply function`);
  });

  test('A saved model is found, listed and suggested for a table with its inputs', async () => {
    const id = await saveModel();
    try {
      const model = await grok.dapi.models.find(id);
      expect(model.id, id);
      expect((await grok.dapi.models.filter(`friendlyName = "${model.friendlyName}"`).list()).length, 1);
      expect((await suggestedFor(grok.data.demo.demog(10).columns.toList())).some((m) => m.id === id), true,
        'not suggested for demog, which has AGE and WEIGHT');
      expect((await suggestedFor([DG.Column.fromStrings('x', ['a'])])).some((m) => m.id === id), false,
        'suggested for a table without its inputs');
    }
    finally {
      await deleteModel(id);
    }
  });

  test('Blobs and images round-trip, and the zip holds them with the model', async () => {
    const id = await saveModel();
    try {
      expect((await ml(`/blobs/${id}?name=&ext=bin`, {method: 'POST', body: new Uint8Array(blob)})).status, 200);
      expect(await bytesOf(await ml(`/blobs/${id}?ext=bin`)), blob.join(','));
      expect((await ml(`/images/roc/${id}?ext=png`, {method: 'POST', body: new Uint8Array(png)})).status, 200);
      const images: string[] = await (await ml(`/images/${id}`)).json();
      expect(images.length, 1);
      expect(images[0].endsWith('roc.png'), true, images[0]);
      expect(await bytesOf(await ml(`/images/${id}/${images[0]}`)), png.join(','));

      const zip = new Uint8Array(await (await ml(`/zip/${id}`)).arrayBuffer());
      expect(Array.from(zip.slice(0, 4)).join(','), '80,75,3,4', 'the zip signature');
      const entries = new TextDecoder('latin1').decode(zip);
      for (const name of [`${id}.json`, images[0], '.bin'])
        expect(entries.includes(name), true, `the zip has no ${name}`);
    }
    finally {
      await deleteModel(id);
    }
  });

  test('A deleted model takes its blobs and images with it', async () => {
    const id = await saveModel();
    try {
      await ml(`/blobs/${id}?name=&ext=bin`, {method: 'POST', body: new Uint8Array(blob)});
      await ml(`/images/roc/${id}?ext=png`, {method: 'POST', body: new Uint8Array(png)});
      await grok.dapi.models.delete(await grok.dapi.models.find(id));
      expect((await grok.dapi.models.find(id)) == null, true, 'the model is still found');
      expect((await bytesOf(await ml(`/blobs/${id}?ext=bin`))) === blob.join(','), false, 'the blob is still served');
      expect((await (await ml(`/images/${id}`)).json()).length, 0);
    }
    finally {
      await deleteModel(id);
    }
  });
}, {owner: 'aparamonov@datagrok.ai'});
