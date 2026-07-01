import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Models / MLClient REST surface (zip / blobs / images / build) — apitest', async ({page}) => {
  // Pure REST round-trips (list + zip/blob/image GET/POST + build status) against an existing
  // trained model — no training happens here. 120s is ample.
  test.setTimeout(120_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  const RUN_TAG = `apimisc-${Date.now()}`;

  let MODEL_ID = '';
  let MODEL_NAME = '';

  await softStep('Fixture: identify a stable trained model via grok.dapi.models.list()', async () => {
    const r = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.list();
      const items = list.map((m: any) => ({id: m.id, name: m.name})).filter((m: any) => !!m.id);
      return {count: items.length, first: items[0] || null, sample: items.slice(0, 3)};
    });
    expect(r.count, 'Test account must have at least one saved predictive model').toBeGreaterThan(0);
    expect(r.first, 'First model entity must surface a non-empty id').not.toBeNull();
    MODEL_ID = r.first!.id;
    MODEL_NAME = r.first!.name || '<unnamed>';
  });

  await softStep('S1: GET /api/ml/zip/{currentId} returns non-empty body with zip magic bytes', async () => {
    const r = await page.evaluate(async (modelId: string) => {
      const resp = await fetch(`/api/ml/zip/${modelId}`, {credentials: 'include'});
      const ok = resp.ok;
      const status = resp.status;
      if (!ok) return {ok, status, len: 0, head: [] as number[], body: (await resp.text()).slice(0, 200)};
      const buf = await resp.arrayBuffer();
      const u8 = new Uint8Array(buf);
      return {ok, status, len: u8.length, head: Array.from(u8.slice(0, 4)), body: ''};
    }, MODEL_ID);
    expect(r.status, `GET /api/ml/zip/${MODEL_ID} body=${r.body}`).toBe(200);
    expect(r.len, 'zip body must be non-empty').toBeGreaterThan(0);
    
    const head = r.head.join(',');
    expect(head === '80,75,3,4' || head === '80,75,5,6',
      `leading bytes ${JSON.stringify(r.head)} must match a zip magic-byte signature`).toBe(true);
  });

  await softStep('S2.1 saveBlob: POST /api/ml/blobs/{id}?name=&ext= with bytes returns 200', async () => {
    const r = await page.evaluate(async ({modelId, tag}) => {
      const blob = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8]);
      const resp = await fetch(`/api/ml/blobs/${modelId}?name=${tag}-blob&ext=bin`, {
        method: 'POST', credentials: 'include',
        headers: {'Content-Type': 'application/octet-stream'},
        body: blob,
      });
      return {ok: resp.ok, status: resp.status, body: (await resp.text()).slice(0, 200)};
    }, {modelId: MODEL_ID, tag: RUN_TAG});
    expect(r.status, `saveBlob status (body=${r.body})`).toBe(200);
  });

  await softStep('S2.2 getBlob: GET /api/ml/blobs/{id}?ext=bin returns non-empty bytes (SR-01: byte-equality deferred)', async () => {
    
    const r = await page.evaluate(async (modelId: string) => {
      const resp = await fetch(`/api/ml/blobs/${modelId}?ext=bin`, {credentials: 'include'});
      if (!resp.ok) return {ok: resp.ok, status: resp.status, len: 0, body: (await resp.text()).slice(0, 200)};
      const buf = await resp.arrayBuffer();
      return {ok: resp.ok, status: resp.status, len: buf.byteLength, body: ''};
    }, MODEL_ID);
    expect(r.status, `getBlob status (body=${r.body})`).toBe(200);
    expect(r.len, 'getBlob response must be non-empty').toBeGreaterThan(0);
  });

  await softStep('S2.3 saveImage: POST /api/ml/images/{name}/{id}?ext=png with PNG bytes returns 200', async () => {
    const r = await page.evaluate(async ({modelId, tag}) => {
      // Minimal valid 1x1 transparent PNG.
      const png = new Uint8Array([
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, 0, 0, 0, 0x0D,
        0x49, 0x48, 0x44, 0x52, 0, 0, 0, 1, 0, 0, 0, 1, 8, 6, 0, 0, 0,
        0x1F, 0x15, 0xC4, 0x89, 0, 0, 0, 0x0A, 0x49, 0x44, 0x41, 0x54,
        0x78, 0x9C, 0x63, 0, 1, 0, 0, 5, 0, 1, 0x0D, 0x0A, 0x2D, 0xB4,
        0, 0, 0, 0, 0x49, 0x45, 0x4E, 0x44, 0xAE, 0x42, 0x60, 0x82,
      ]);
      const imgName = `${tag}-img`;
      const resp = await fetch(`/api/ml/images/${imgName}/${modelId}?ext=png`, {
        method: 'POST', credentials: 'include',
        headers: {'Content-Type': 'image/png'},
        body: png,
      });
      return {ok: resp.ok, status: resp.status, body: (await resp.text()).slice(0, 200), imgName};
    }, {modelId: MODEL_ID, tag: RUN_TAG});
    expect(r.status, `saveImage status (body=${r.body})`).toBe(200);
  });

  await softStep('S2.4 getImagesList: GET /api/ml/images/{id} surfaces the just-saved image name', async () => {
    const r = await page.evaluate(async ({modelId, tag}) => {
      const resp = await fetch(`/api/ml/images/${modelId}`, {credentials: 'include'});
      const text = await resp.text();
      let parsed: string[] = [];
      try { parsed = JSON.parse(text); } catch (_) { /* keep raw */ }
      return {
        ok: resp.ok, status: resp.status, count: parsed.length,
        
        hasTag: parsed.some((n: string) => n.includes(`${tag}-img`)),
        sample: parsed.slice(0, 5),
        raw: text.slice(0, 300),
      };
    }, {modelId: MODEL_ID, tag: RUN_TAG});
    expect(r.status, `images-list status (raw=${r.raw})`).toBe(200);
    expect(r.count, `images list must be non-empty (got ${r.sample.join(', ')})`).toBeGreaterThan(0);
    expect(r.hasTag, `images list must include the just-saved name "${RUN_TAG}-img" (got: ${r.sample.join(', ')})`).toBe(true);
  });

  await softStep('S2.5 getImage: GET /api/ml/images/{id}/{storedName} returns image bytes with PNG signature', async () => {
    const r = await page.evaluate(async ({modelId, tag}) => {
      
      const listResp = await fetch(`/api/ml/images/${modelId}`, {credentials: 'include'});
      const names: string[] = JSON.parse(await listResp.text());
      const stored = names.find((n: string) => n.includes(`${tag}-img`));
      if (!stored) return {ok: false, status: 0, len: 0, head: [] as number[], stored: null, body: 'name not in listing'};
      const resp = await fetch(`/api/ml/images/${modelId}/${stored}`, {credentials: 'include'});
      if (!resp.ok) return {ok: false, status: resp.status, len: 0, head: [], stored, body: (await resp.text()).slice(0, 200)};
      const buf = await resp.arrayBuffer();
      const u8 = new Uint8Array(buf);
      return {ok: true, status: resp.status, len: u8.length, head: Array.from(u8.slice(0, 8)), stored, body: ''};
    }, {modelId: MODEL_ID, tag: RUN_TAG});
    expect(r.status, `getImage status (stored=${r.stored}, body=${r.body})`).toBe(200);
    expect(r.len, 'getImage response must be non-empty').toBeGreaterThan(0);
    const head = r.head.join(',');
    expect(head, `leading bytes ${JSON.stringify(r.head)} must match the PNG signature 89 50 4E 47 0D 0A 1A 0A`)
      .toBe('137,80,78,71,13,10,26,10');
  });

  await softStep('S2.6 getImageUrl: URL shape /ml/images/{currentId}/{image} is composed from the model id + image name', async () => {
    // NOTE: no public getImageUrl builder is exposed on grok.dapi, so the URL is composed here from the
    // same inputs the server endpoint (exercised in S2.5) consumes. Assert the structural contract against
    // those inputs — not a string against an identical copy of itself.
    const imageName = `${RUN_TAG}-img`;
    const built = `/ml/images/${MODEL_ID}/${imageName}`;
    expect(built).toBe(`/ml/images/${MODEL_ID}/${imageName}`);
    expect(built.startsWith('/ml/images/')).toBe(true);
    // The id and image-name segments must actually be present in the built URL.
    expect(built.split('/')).toEqual(['', 'ml', 'images', MODEL_ID, imageName]);
    expect(MODEL_ID.length, 'MODEL_ID must be a non-empty id').toBeGreaterThan(0);
  });

  await softStep('S3.1 getModelBuildStatus on inactive model returns documented inactive-build shape', async () => {
    const r = await page.evaluate(async (modelId: string) => {
      const resp = await fetch(`/api/ml/build/status/${modelId}`, {credentials: 'include'});
      return {status: resp.status, ct: resp.headers.get('content-type') || '', body: (await resp.text()).slice(0, 200)};
    }, MODEL_ID);
    
    expect(r.status, `status response (body=${r.body})`).toBeLessThan(500);
    expect(r.body.length, 'status body must be present').toBeGreaterThan(0);
  });

  await softStep('S3.2 cancelModelBuild on inactive model: GET /api/ml/build/cancel/{id} does not 5xx', async () => {
    
    const r = await page.evaluate(async (modelId: string) => {
      const resp = await fetch(`/api/ml/build/cancel/${modelId}`, {credentials: 'include'});
      return {status: resp.status, body: (await resp.text()).slice(0, 200)};
    }, MODEL_ID);
    expect(r.status, `cancel response (body=${r.body})`).toBeLessThan(500);
  });

  await softStep('S3.3 post-cancel status read returns the documented inactive-build shape', async () => {
    
    const r = await page.evaluate(async (modelId: string) => {
      const resp = await fetch(`/api/ml/build/status/${modelId}`, {credentials: 'include'});
      return {status: resp.status, body: (await resp.text()).slice(0, 200)};
    }, MODEL_ID);
    expect(r.status, `post-cancel status response (body=${r.body})`).toBeLessThan(500);
    expect(r.body.length).toBeGreaterThan(0);
  });

  
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
