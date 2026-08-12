/** Local-file uploads: the pending-bytes registry, the `readUploadedFile` node function,
 *  server persistence, and the creation-script round-trip. Server tests need a live stand. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {
  MAX_UPLOAD_SIZE, addPendingFile, getPendingFile, removePendingFile,
  isPendingFileId, persistPendingFile, parseFileToDataFrame,
  uploadedFileIdsFromFlowBody, syncFlowFilePermissions,
} from '../utils/uploaded-files';
import {readUploadedFile} from '../package';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {buildCreationScriptGraph, applyGraphToEditor} from '../import/creation-script-importer';
import {emitCreationScript} from '../compiler/creation-script-emitter';
import {flowScriptText} from '../serialization/flow-script-format';
import {makeEditor, destroyEditor, addNode} from './test-utils';

const CSV = 'a,b\n1,2\n3,4';
const csvBytes = (): Uint8Array => new TextEncoder().encode(CSV);

category('Flow: uploaded files', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('pending registry: add, read, remove, 100 MB cap', async () => {
    const id = addPendingFile('test.csv', csvBytes());
    try {
      expect(isPendingFileId(id), true, 'fresh id is a pending id');
      expect(getPendingFile(id)?.name, 'test.csv', 'bytes retrievable by id');
    } finally {
      removePendingFile(id);
    }
    expect(getPendingFile(id) == null, true, 'removed from the registry');

    // Boundary check without allocating 100 MB — the cap only reads `length`.
    const oversized = {length: MAX_UPLOAD_SIZE + 1} as Uint8Array;
    let threw = '';
    try {
      addPendingFile('big.csv', oversized);
    } catch (e: any) {
      threw = String(e?.message ?? e);
    }
    expect(threw.includes('100 MB'), true, `cap error names the limit: "${threw}"`);
  });

  test('parseFileToDataFrame parses csv and names the table after the file', async () => {
    const df = await parseFileToDataFrame('demo data.csv', csvBytes());
    expect(df.rowCount, 2);
    expect(df.columns.length, 2);
    expect(df.name, 'demo data', 'name derived from the file name');
  });

  test('parseFileToDataFrame parses sdf through the registered file handler', async () => {
    const bytes = await grok.dapi.files.readAsBytes('System:AppData/Chem/mol1K.sdf');
    const df = await parseFileToDataFrame('mol1K.sdf', bytes);
    expect(df.rowCount > 0, true, 'rows parsed from the sdf');
    expect(df.columns.length > 0, true, 'columns parsed from the sdf');
    expect(df.name, 'mol1K', 'name derived from the file name');
  });

  test('an uploaded sdf runs through the Uploaded File node function', async () => {
    const bytes = await grok.dapi.files.readAsBytes('System:AppData/Chem/mol1K.sdf');
    const id = addPendingFile('mol1K.sdf', bytes);
    try {
      const df = await readUploadedFile(id, 'mol1K.sdf');
      expect(df.rowCount > 0, true, 'rows parsed from the pending sdf');
    } finally {
      removePendingFile(id);
    }
  });

  test('OpenFile opens a server sdf by path (the Open File node contract)', async () => {
    const df: DG.DataFrame = await grok.functions.call('OpenFile',
      {fullPath: 'System:AppData/Chem/mol1K.sdf'});
    expect(df != null, true, 'OpenFile returned a table');
    expect(df.rowCount > 0, true, 'rows parsed from the sdf');
  });

  test('readUploadedFile serves pending bytes with no server round-trip', async () => {
    const id = addPendingFile('pending.csv', csvBytes());
    try {
      const df = await readUploadedFile(id, 'pending.csv');
      expect(df.rowCount, 2);
      expect(df.columns.names().join(','), 'a,b');
    } finally {
      removePendingFile(id);
    }
  });

  test('a lost pending id produces an actionable error', async () => {
    let msg = '';
    try {
      await readUploadedFile('pending:no-such-id', 'gone.csv');
    } catch (e: any) {
      msg = String(e?.message ?? e);
    }
    expect(msg.includes('gone.csv'), true, 'names the file');
    expect(msg.toLowerCase().includes('drop it onto the canvas again'), true,
      `tells the user how to recover: "${msg}"`);
  });

  test('catalog: the Uploaded File node registers via meta.includeInFlow', async () => {
    const info = getRegisteredFuncs().find((f) => f.func.name === 'readUploadedFile');
    expect(info != null, true, 'readUploadedFile is a registered node type');
  });

  test('round-trip: persist to the server blob store, read back, clean up', async () => {
    const id = addPendingFile('flow-upload-test.csv', csvBytes());
    const fi = await persistPendingFile(id);
    const entity = await grok.dapi.entities.find(fi.id).catch(() => null);
    try {
      expect(isPendingFileId(fi.id), false, 'persisted id is a real server id');
      expect(getPendingFile(id) == null, true, 'pending entry consumed');
      // files_service.dart adds connection-less files to the author's project — permission sync depends on this.
      expect(entity != null, true, 'FileInfo entity visible to its author');

      const bytes = await grok.dapi.files.readAsBytes(fi.id);
      expect(new TextDecoder().decode(bytes), CSV, 'blob bytes survive the round-trip');

      const df = await readUploadedFile(fi.id, 'flow-upload-test.csv');
      expect(df.rowCount, 2, 'node function reads the persisted file');
    } finally {
      if (entity != null)
        await grok.dapi.entities.delete(entity);
    }
  });

  test('share sync: flow grantees get read on the uploaded files', async () => {
    const pid = addPendingFile('share-sync-test.csv', csvBytes());
    const fi = await persistPendingFile(pid);
    let script: DG.Script | null = null;
    try {
      const e = makeEditor();
      let body: string;
      try {
        const info = getRegisteredFuncs().find((f) => f.func.name === 'readUploadedFile');
        expect(info != null, true, 'Uploaded File node registered');
        const node = await addNode(e.flow, info!.nodeTypeName);
        node.inputValues['fileId'] = fi.id;
        node.inputValues['fileName'] = 'share-sync-test.csv';
        body = flowScriptText(e.flow,
          {scriptName: 'Flow share sync test', scriptDescription: '', tags: ['flow']});
      } finally {
        destroyEditor(e);
      }
      expect(uploadedFileIdsFromFlowBody(body).join(','), fi.id, 'body parser finds the file id');

      // Run the sync the flowShareSync autostart performs on d4-entity-shared.
      script = await grok.dapi.scripts.save(DG.Script.create(body));
      const allUsers = await grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
      await grok.dapi.permissions.grant(script, allUsers, false);
      await syncFlowFilePermissions(script);

      const fileEntity = await grok.dapi.entities.find(fi.id);
      const perms = await grok.dapi.permissions.get(fileEntity) as
        unknown as {view?: DG.Group[]; edit?: DG.Group[]};
      const granted = [...(perms.view ?? []), ...(perms.edit ?? [])]
        .some((g) => g.id === allUsers.id);
      expect(granted, true, 'All users got read on the file entity');
    } finally {
      if (script?.id != null)
        await grok.dapi.scripts.delete(script).catch(() => {});
      const entity = await grok.dapi.entities.find(fi.id).catch(() => null);
      if (entity != null)
        await grok.dapi.entities.delete(entity).catch(() => {});
    }
  });

  test('creation script: Flow:readUploadedFile round-trips through import/emit', async () => {
    const func = DG.Func.find({name: 'readUploadedFile'})[0];
    expect(func != null, true, 'function is published');
    const call = `T = ${func.nqName}("00000000-0000-0000-0000-000000000000", "data.csv")`;
    const e = makeEditor();
    try {
      const g = buildCreationScriptGraph(call);
      await applyGraphToEditor(g, e.flow);
      const r = emitCreationScript(e.flow);
      expect(r.warnings.length, 0, r.warnings.join(' ; '));
      expect(r.script.toLowerCase().includes('readuploadedfile('), true, `emitted: ${r.script}`);
      expect(r.script.includes('00000000-0000-0000-0000-000000000000'), true, 'file id survives');
      expect(r.script.includes('"data.csv"'), true, 'file name survives');
    } finally {
      destroyEditor(e);
    }
  });
});
