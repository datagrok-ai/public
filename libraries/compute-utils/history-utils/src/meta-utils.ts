import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// schema

type UserMeta = {
  authorName: string,
  authorId: string,
  title?: string,
  description?: string,
  tags?: string,
  approvalRequested?: boolean,
  approvalRequestTimestamp?: number,
  isDeleted?: boolean,
}

type ApprovalStatus = 'Approved' | 'Rejected' | 'Recalled';

type SMEMeta = {
  status: ApprovalStatus,
  authorName: string,
  authorId: string,
  reviewerName: string,
  reviewerId: string,
  resolutionTimestamp: number,
  title?: string,
  description?: string,
  tags?: string,
  rejectionReason?: string,
  recalledByName?: string,
  recalledById?: string,
  recalledByReason?: string,
  recallTimestamp?: number,
}

const sourceTagName = `source`;
const metaColName = `run_id`;

function getNames(modelName: string) {
  const modelUserSchemaName = `ModelHub-User-${modelName}-Schema`;
  const modelSMESchemaName = `ModelHub-SME-${modelName}-Schema`;
  const runColNameUser = `ModelHub-User-RunId-${modelName}-Col`;
  const runColNameSME = `ModelHub-SME-RunId-${modelName}-Col`;
  const runTagUser = `ModelHub-User-RunId-${modelName}-Tag`;
  const runTagSME = `ModelHub-SME-RunId-${modelName}-Tag`;
  return {modelUserSchemaName, modelSMESchemaName, runTagUser, runTagSME, runColNameUser, runColNameSME};
}

async function getOrCreateModelStickyMeta(modelName: string) {
  const {modelUserSchemaName, modelSMESchemaName, runTagUser, runTagSME, runColNameUser, runColNameSME} = getNames(modelName);

  const schemas = await grok.dapi.stickyMeta.getSchemas();
  const userSchema = schemas.find((s) => s.name === modelUserSchemaName);
  const SMESchema = schemas.find((s) => s.name === modelSMESchemaName);

  if (userSchema && SMESchema)
    return {userSchema, SMESchema};


  const modelUserSchema = await grok.dapi.stickyMeta.createSchema(modelUserSchemaName,
    [{name: runColNameUser, matchBy: `source=${runTagUser}`}],
    [
      {name: 'authorName', type: 'string'},
      {name: 'authorId', type: 'string'},

      {name: 'title', type: 'string'},
      {name: 'description', type: 'string'},
      {name: 'tags', type: 'string'},

      {name: 'approvalRequested', type: 'bool'},
      {name: 'approvalRequestTimestamp', type: 'int'},
      {name: 'isDeleted', type: 'bool'},
    ],
  );

  const modelSMESchema = await grok.dapi.stickyMeta.createSchema(modelSMESchemaName,
    [{name: runColNameSME, matchBy: `source=${runTagSME}`}],
    [
      {name: 'status', type: 'string'}, // ApprovalStatus
      {name: 'authorName', type: 'string'},
      {name: 'authorId', type: 'string'},
      {name: 'reviewerName', type: 'string'},
      {name: 'reviewerId', type: 'string'},
      {name: 'resolutionTimestamp', type: 'int'},

      {name: 'title', type: 'string'},
      {name: 'description', type: 'string'},
      {name: 'tags', type: 'string'},

      {name: 'rejectionReason', type: 'string'},
      {name: 'recalledByName', type: 'string'},
      {name: 'recalledById', type: 'string'},
      {name: 'recalledByReason', type: 'string'},
      {name: 'recallTimestamp', type: 'int'},
    ],
  );
  return {userSchema: modelUserSchema, SMESchema: modelSMESchema};
}

function getUserMetaDF(meta: UserMeta) {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'authorName', [meta.authorName]),
    DG.Column.fromList('string', 'authorId', [meta.authorId]),
    DG.Column.fromList('string', 'title', [meta.title ?? '']),
    DG.Column.fromList('string', 'description', [meta.description ?? '']),
    DG.Column.fromList('string', 'tags', [meta.tags ?? '']),
    DG.Column.fromList('bool', 'approvalRequested', [meta.approvalRequested ?? false]),
    DG.Column.fromList('int', 'approvalRequestTimestamp', [meta.approvalRequestTimestamp ?? 0]),
    DG.Column.fromList('bool', 'isDeleted', [meta.isDeleted ?? false]),
  ]);
}

function getSMEMetaDF(meta: SMEMeta) {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'status', [meta.status]),
    DG.Column.fromList('string', 'authorName', [meta.authorName]),
    DG.Column.fromList('string', 'authorId', [meta.authorId]),
    DG.Column.fromList('string', 'reviewerName', [meta.reviewerName]),
    DG.Column.fromList('string', 'reviewerId', [meta.reviewerId]),
    DG.Column.fromList('int', 'resolutionTimestamp', [meta.resolutionTimestamp]),

    DG.Column.fromList('string', 'title', [meta.title ?? '']),
    DG.Column.fromList('string', 'description', [meta.description ?? '']),
    DG.Column.fromList('string', 'tags', [meta.tags ?? '']),

    DG.Column.fromList('string', 'rejectionReason', [meta.rejectionReason ?? '']),
    DG.Column.fromList('string', 'recalledByName', [meta.recalledByName ?? '']),
    DG.Column.fromList('string', 'recalledById', [meta.recalledById ?? '']),
    DG.Column.fromList('string', 'recalledByReason', [meta.recalledByReason ?? '']),
    DG.Column.fromList('int', 'recallTimestamp', [meta.recallTimestamp ?? 0]),
  ]);
}

function getUserMetaCol(modelName: string, ids: string[]) {
  const {runTagUser} = getNames(modelName);
  const metaCol = DG.Column.fromList('string', metaColName, ids);
  metaCol.setTag(sourceTagName, runTagUser);
  return metaCol;
}

function getSMEMetaCol(modelName: string, ids: string[]) {
  const {runTagSME} = getNames(modelName);
  const metaCol = DG.Column.fromList('string', metaColName, ids);
  metaCol.setTag(sourceTagName, runTagSME);
  return metaCol;
}

async function attachUserStickyMeta(modelName: string, id: string, meta: UserMeta) {
  const {userSchema} = await getOrCreateModelStickyMeta(modelName);
  const metaCol = getUserMetaCol(modelName, [id]);
  const metaDF = getUserMetaDF(meta);
  await grok.dapi.stickyMeta.setAllValues(userSchema, metaCol, metaDF);
}

async function getUserStickyMeta(modelName: string, ids: string[]) {
  const {userSchema} = await getOrCreateModelStickyMeta(modelName);
  const metaCol = getUserMetaCol(modelName, ids);
  const meta = await grok.dapi.stickyMeta.getAllValues(userSchema, metaCol);
  return meta;
}

async function attachSMEStickyMeta(modelName: string, id: string, meta: SMEMeta) {
  const {SMESchema} = await getOrCreateModelStickyMeta(modelName);
  const metaCol = getSMEMetaCol(modelName, [id]);
  const metaDF = getSMEMetaDF(meta);
  await grok.dapi.stickyMeta.setAllValues(SMESchema, metaCol, metaDF);
}

async function getSMEStickyMeta(modelName: string, ids: string[]) {
  const {SMESchema} = await getOrCreateModelStickyMeta(modelName);
  const metaCol = getSMEMetaCol(modelName, ids);
  const meta = await grok.dapi.stickyMeta.getAllValues(SMESchema, metaCol);
  return meta;
}


// test workflow

async function makeMockFuncCall() {
  const func = await grok.functions.eval('Libtests:TestAdd2') as DG.Func;
  const fc = func.prepare({a: 1, b: 2});
  fc.newId();
  return fc;
}


async function testUserWorkflow() {
  const modelName = 'MyMockModel';
  const fc = await makeMockFuncCall();
  console.log(`created FC ${fc.id}`);
  // work around bug
  const id = `UUID: ${fc.id}`;
  await attachUserStickyMeta(modelName, id, {title: 'Some Title', authorId: 'MockId', authorName: 'MockName'});
  console.log(`saved FC user meta ${fc.id}`);
  console.log(`starting loading FC ${fc.id}`);
  // TODO: reader should handle absent stiky metas
  const meta = await getUserStickyMeta(modelName, [id]);
  console.log(`loaded FC ${fc.id} user meta:\n${meta.toCsv()}`);
}
