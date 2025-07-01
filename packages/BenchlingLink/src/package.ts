/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { queryAASequences, postAASequence } from './aaSequencesApi';
import { queryDNASequences, postDNASequence } from './dnaSequencesApi';
import { queryAssayResults, queryAssayRuns, postAssayResult, postAssayRun } from './assayApi';

export const _package = new DG.Package();

//tags: app
//name: Benchling Link
//output: view v
//meta.browsePath: Chem
export async function benchlingLinkApp(): Promise<DG.ViewBase> {

  return DG.View.fromRoot(ui.div('Pilot version of Benchling link application'));
}

//input: dynamic treeNode
//input: view browseView
export async function benchlingLinkAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  try {

    const aaSequenceNode = treeNode.group('AA Sequences', DG.toDart(DG.Func.byName('BenchlingLink:getAASequences')));
    aaSequenceNode.item('Create AA Sequence', DG.toDart(DG.Func.byName('BenchlingLink:createAASequence')));

    const dnaSequencesNode = treeNode.group('DNA Sequences', DG.toDart(DG.Func.byName('BenchlingLink:getDNASequences')));
    dnaSequencesNode.item('Create DNA Sequence', DG.toDart(DG.Func.byName('BenchlingLink:createDNASequence')));

    const assayResultsNode = treeNode.group('Assay Results', DG.toDart(DG.Func.byName('BenchlingLink:getAssayResults')));
    assayResultsNode.item('Create Assay Result', DG.toDart(DG.Func.byName('BenchlingLink:createAssayResult')));

    const assayRunsNode = treeNode.group('Assay Runs', DG.toDart(DG.Func.byName('BenchlingLink:getAssayRuns')));
    assayRunsNode.item('Create Assay Run', DG.toDart(DG.Func.byName('BenchlingLink:createAssayRun')));
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
  }
}

//name: Get AA Sequences
//input: string sort {nullable:true}
//input: string createdAt {nullable:true}
//input: string modifiedAt {nullable:true}
//input: string name {nullable:true}
//input: string nameIncludes {nullable:true}
//input: string aminoAcids {nullable:true}
//input: string folderId {nullable:true}
//input: string mentionedIn {nullable:true}
//input: string projectId {nullable:true}
//input: string registryId {nullable:true}
//input: string schemaId {nullable:true}
//input: string schemaFields {nullable:true}
//input: string archiveReason {nullable:true}
//input: string mentions {nullable:true}
//input: string ids {nullable:true}
//input: string entityRegistryIds_anyOf {nullable:true}
//input: string names_anyOf {nullable:true}
//input: string names_anyOf_caseSensitive {nullable:true}
//input: string creatorIds {nullable:true}
//input: string authorIds_anyOf {nullable:true}
//input: string returning {nullable:true}
//output: dataframe df
export async function getAASequences(
  sort?: string,
  createdAt?: string,
  modifiedAt?: string,
  name?: string,
  nameIncludes?: string,
  aminoAcids?: string,
  folderId?: string,
  mentionedIn?: string,
  projectId?: string,
  registryId?: string,
  schemaId?: string,
  schemaFields?: string,
  archiveReason?: string,
  mentions?: string,
  ids?: string,
  entityRegistryIds_anyOf?: string,
  names_anyOf?: string,
  names_anyOf_caseSensitive?: string,
  creatorIds?: string,
  authorIds_anyOf?: string,
  returning?: string
): Promise<DG.DataFrame> {
  const params = {
    sort,
    createdAt,
    modifiedAt,
    name,
    nameIncludes,
    aminoAcids,
    folderId,
    mentionedIn,
    projectId,
    registryId,
    schemaId,
    schemaFields,
    archiveReason,
    mentions,
    ids,
    entityRegistryIds_anyOf,
    names_anyOf,
    names_anyOf_caseSensitive,
    creatorIds,
    authorIds_anyOf,
    returning,
  };
  const aaSequences = await queryAASequences(params);
  return aaSequences;
}

//name: Get DNA Sequences
//input: string sort {nullable:true}
//input: string createdAt {nullable:true}
//input: string modifiedAt {nullable:true}
//input: string name {nullable:true}
//input: string nameIncludes {nullable:true}
//input: string bases {nullable:true}
//input: string folderId {nullable:true}
//input: string mentionedIn {nullable:true}
//input: string projectId {nullable:true}
//input: string registryId {nullable:true}
//input: string schemaId {nullable:true}
//input: string schemaFields {nullable:true}
//input: string archiveReason {nullable:true}
//input: string mentions {nullable:true}
//input: string ids {nullable:true}
//input: string entityRegistryIds_anyOf {nullable:true}
//output: dataframe df
export async function getDNASequences(
  sort?: string,
  createdAt?: string,
  modifiedAt?: string,
  name?: string,
  nameIncludes?: string,
  bases?: string,
  folderId?: string,
  mentionedIn?: string,
  projectId?: string,
  registryId?: string,
  schemaId?: string,
  schemaFields?: string,
  archiveReason?: string,
  mentions?: string,
  ids?: string,
  entityRegistryIds_anyOf?: string
): Promise<DG.DataFrame> {
  const params = {
    sort,
    createdAt,
    modifiedAt,
    name,
    nameIncludes,
    bases,
    folderId,
    mentionedIn,
    projectId,
    registryId,
    schemaId,
    schemaFields,
    archiveReason,
    mentions,
    ids,
    entityRegistryIds_anyOf,
  };
  const dnaSequences = await queryDNASequences(params);
  return dnaSequences;
}

//name: Get Assay Results
//input: string schemaId {nullable:true}
//input: string createdAt_lt {nullable:true}
//input: string createdAt_gt {nullable:true}
//input: string createdAt_lte {nullable:true}
//input: string createdAt_gte {nullable:true}
//input: int minCreatedTime {nullable:true}
//input: int maxCreatedTime {nullable:true}
//input: string sort {nullable:true}
//input: string entityIds {nullable:true}
//input: string storageIds {nullable:true}
//input: string assayRunIds {nullable:true}
//input: string automationOutputProcessorId {nullable:true}
//input: string ids {nullable:true}
//input: string modifiedAt_lt {nullable:true}
//input: string modifiedAt_gt {nullable:true}
//input: string modifiedAt_lte {nullable:true}
//input: string modifiedAt_gte {nullable:true}
//input: string archiveReason {nullable:true}
//output: dataframe df
export async function getAssayResults(
  schemaId?: string,
  createdAt_lt?: string,
  createdAt_gt?: string,
  createdAt_lte?: string,
  createdAt_gte?: string,
  minCreatedTime?: number,
  maxCreatedTime?: number,
  sort?: string,
  entityIds?: string,
  storageIds?: string,
  assayRunIds?: string,
  automationOutputProcessorId?: string,
  ids?: string,
  modifiedAt_lt?: string,
  modifiedAt_gt?: string,
  modifiedAt_lte?: string,
  modifiedAt_gte?: string,
  archiveReason?: string
): Promise<DG.DataFrame> {
  const params = {
    schemaId,
    createdAt_lt,
    createdAt_gt,
    createdAt_lte,
    createdAt_gte,
    minCreatedTime,
    maxCreatedTime,
    sort,
    entityIds,
    storageIds,
    assayRunIds,
    automationOutputProcessorId,
    ids,
    modifiedAt_lt,
    modifiedAt_gt,
    modifiedAt_lte,
    modifiedAt_gte,
    archiveReason,
  };
  const assayResults = await queryAssayResults(params);
  return assayResults;
}

//name: Get Assay Runs
//input: string schemaId {nullable:true}
//input: int minCreatedTime {nullable:true}
//input: int maxCreatedTime {nullable:true}
//input: string ids {nullable:true}
//output: dataframe df
export async function getAssayRuns(
  schemaId?: string,
  minCreatedTime?: number,
  maxCreatedTime?: number,
  ids?: string
): Promise<DG.DataFrame> {
  const params = {
    schemaId,
    minCreatedTime,
    maxCreatedTime,
    ids,
  };
  const assayRuns = await queryAssayRuns(params);
  return assayRuns;
}

//name: Create AA Sequence
//input: string name
//input: string aminoAcids
//input: string aliases {nullable:true}
//input: string annotations {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//input: string fields {nullable:true}
//input: string folderId {nullable:true}
//input: string schemaId {nullable:true}
//output: dataframe df
export async function createAASequence(
  name: string,
  aminoAcids: string,
  aliases?: string,
  annotations?: string,
  authorIds?: string,
  customFields?: string,
  fields?: string,
  folderId?: string,
  schemaId?: string,
): Promise<DG.DataFrame> {
  // Parse JSON string inputs for array/object fields if provided
  const body: any = { name, aminoAcids };
  if (aliases) body.aliases = JSON.parse(aliases);
  if (annotations) body.annotations = JSON.parse(annotations);
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  if (fields) body.fields = JSON.parse(fields);
  if (folderId) body.folderId = folderId;
  if (schemaId) body.schemaId = schemaId;
  const result = await postAASequence(body);
  return DG.DataFrame.fromObjects([result]) ?? DG.DataFrame.create();
}

//name: Create DNA Sequence
//input: string name
//input: string bases
//input: string aliases {nullable:true}
//input: string annotations {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//input: string fields {nullable:true}
//input: string folderId {nullable:true}
//input: string schemaId {nullable:true}
//output: dataframe df
export async function createDNASequence(
  name: string,
  bases: string,
  aliases?: string,
  annotations?: string,
  authorIds?: string,
  customFields?: string,
  fields?: string,
  folderId?: string,
  schemaId?: string,
): Promise<DG.DataFrame> {
  const body: any = { name, bases };
  if (aliases) body.aliases = JSON.parse(aliases);
  if (annotations) body.annotations = JSON.parse(annotations);
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  if (fields) body.fields = JSON.parse(fields);
  if (folderId) body.folderId = folderId;
  if (schemaId) body.schemaId = schemaId;
  const result = await postDNASequence(body);
  return DG.DataFrame.fromObjects([result]) ?? DG.DataFrame.create();
}

//name: Create Assay Result
//input: string schemaId
//input: string fields {nullable:true}
//input: string entityIds {nullable:true}
//input: string storageIds {nullable:true}
//input: string assayRunId {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//output: dataframe df
export async function createAssayResult(
  schemaId: string,
  fields?: string,
  entityIds?: string,
  storageIds?: string,
  assayRunId?: string,
  authorIds?: string,
  customFields?: string,
): Promise<DG.DataFrame> {
  const body: any = { schemaId };
  if (fields) body.fields = JSON.parse(fields);
  if (entityIds) body.entityIds = JSON.parse(entityIds);
  if (storageIds) body.storageIds = JSON.parse(storageIds);
  if (assayRunId) body.assayRunId = assayRunId;
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  const result = await postAssayResult(body);
  return DG.DataFrame.fromObjects([result]) ?? DG.DataFrame.create();
}

//name: Create Assay Run
//input: string schemaId
//input: string fields {nullable:true}
//input: string name {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//output: dataframe df
export async function createAssayRun(
  schemaId: string,
  fields?: string,
  name?: string,
  authorIds?: string,
  customFields?: string,
): Promise<DG.DataFrame> {
  const body: any = { schemaId };
  if (fields) body.fields = JSON.parse(fields);
  if (name) body.name = name;
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  const result = await postAssayRun(body);
  return DG.DataFrame.fromObjects([result]) ?? DG.DataFrame.create();
}