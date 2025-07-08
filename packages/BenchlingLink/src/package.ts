/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { queryAASequences, postAASequence } from './aaSequencesApi';
import { queryDNASequences, postDNASequence } from './dnaSequencesApi';
import { queryAssayResults, queryAssayRuns, postAssayResult, postAssayRun } from './assayApi';
import { queryMolecules, postMolecule } from './moleculesApi';
import { dataFrameFromObjects } from './utils';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { queryPlates, postPlate } from './platesApi';
import { queryMixtures, postMixture } from './mixturesApi';
import { queryDnaOligos, postDnaOligo } from './dnaOligosApi';

export const _package = new DG.Package();
const STORAGE_NAME = 'BenchlingLinkFuncEditor';
let openedView: DG.View | null = null;

//tags: app
//name: Benchling
//output: view v
//meta.browsePath: Chem
export async function benchlingLinkApp(): Promise<DG.ViewBase> {

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/benchling.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/BenchlingLink/README.md',
    description: '- Integrate with your Benchling account.\n' +
      '- Analyze assay data.\n' +
      '- Find contextual information on molecules and sequences.\n' +
      '- Create entities and post them to the registry.\n' +
      '- Browse the tenant content.\n'
  });

  const view = DG.View.fromRoot(appHeader);
  view.name = 'Benchling';
  return view;
}

//input: dynamic treeNode
//input: view browseView
export async function benchlingLinkAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  function createFuncEditorView(funcName: string, v: DG.View) {
    let func = DG.Func.byName(funcName);
    let editorDiv = ui.div();
    let gridDiv = ui.div();
    let root = ui.splitV([editorDiv, gridDiv], {style: {height: '100%', width: '100%'}}, true);
    gridDiv.setAttribute('style', 'overflow:hidden !important');
    let df: DG.DataFrame | null;

    const addToWorkspaceButton = ui.icons.add(() => {
      if (df) {
        grok.shell.addTablePreview(df);
      }
    }, 'Add results to workspace');
    v.setRibbonPanels([[addToWorkspaceButton]]);

    // Restore state if exists
    const saved = grok.userSettings.getValue(STORAGE_NAME, funcName);
    const funcCall = func.prepare(saved ? JSON.parse(saved) : null);
    renderEditorAndRun();

    async function renderEditorAndRun() {
      ui.empty(editorDiv);
      let form = await DG.InputForm.forFuncCall(funcCall);
      form.root.style.flexWrap = 'wrap';
      form.root.style.height = '100%';
      form.root.style.maxWidth = 'unset';
      DG.debounce(form.onInputChanged, 1000).subscribe(() => {
        runFunc();
      })
      editorDiv.appendChild(form.root);
      runFunc();
    }

    async function runFunc() {
      //save request settings
      const params: {[key: string]: any} = {};
      Object.keys(funcCall.inputParams).forEach((paramName) => {
        if (funcCall.inputParams[paramName].value)
          params[paramName] = funcCall.inputParams[paramName].value;
      });    
      grok.userSettings.add(STORAGE_NAME, funcName, JSON.stringify(params));

      //run function
      ui.empty(gridDiv);
      try {
        df = (await funcCall.call()).getOutputParamValue();
        let grid = df!.plot.grid();
        grid.root.style.width = '100%';
        grid.root.style.height = '95%';
        gridDiv.appendChild(grid.root);
      } catch (e) {
        gridDiv.appendChild(ui.divText('Error: ' + e));
      }
    }

    return root;
  }

  const addBenchlingView = (viewName: string, funcName: string) => {
    if (openedView)
      openedView.close();
    openedView = DG.View.create(viewName);
    openedView.name = viewName;
    openedView.root.appendChild(createFuncEditorView(funcName, openedView));
    const path = ['Home', 'Benchling', viewName];
    grok.shell.addPreview(openedView);
    setBreadcrumbs(path, treeNode);
  }

  const setBreadcrumbs = (path: string[], tree: DG.TreeViewGroup) => {
    const breadcrumbs = ui.breadcrumbs(path);
    breadcrumbs.onPathClick.subscribe(async (value) => {
      const actualItem = value[value.length - 1];
      if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
        return;
      if (actualItem === 'Benchling')
        tree.currentItem = tree;
    });
    if (grok.shell.v) {
      if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
        const homeIcon = ui.iconFA('home', () => {
          grok.shell.v.close();
          grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
        });
        breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
      }
      const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
      if (viewNameRoot) {
        viewNameRoot.textContent = '';
        viewNameRoot.appendChild(breadcrumbs.root);
      }
    }
  }

  try {
    treeNode.items.length = 0;
    const aaSequenceNode = treeNode.item('AA Sequences');
    aaSequenceNode.onSelected.subscribe(async () => {
      addBenchlingView('AA Sequences', 'BenchlingLink:getAASequences');
    });

    const dnaSequencesNode = treeNode.item('DNA Sequences');
    dnaSequencesNode.onSelected.subscribe(async () => {
      addBenchlingView('DNA Sequences', 'BenchlingLink:getDNASequences');
    });

    const assayResultsNode = treeNode.item('Assay Results');
    assayResultsNode.onSelected.subscribe(async () => {
      addBenchlingView('Assay Results', 'BenchlingLink:getAssayResults');
    });

    const assayRunsNode = treeNode.item('Assay Runs');
    assayRunsNode.onSelected.subscribe(async () => {
      addBenchlingView('Assay Runs', 'BenchlingLink:getAssayRuns');
    });

    const moleculesNode = treeNode.item('Molecules');
    moleculesNode.onSelected.subscribe(async () => {
      addBenchlingView('Molecules', 'BenchlingLink:getMolecules');
    });

    const platesNode = treeNode.item('Plates');
    platesNode.onSelected.subscribe(async () => {
      addBenchlingView('Plates', 'BenchlingLink:getPlates');
    });

    const mixturesNode = treeNode.item('Mixtures');
    mixturesNode.onSelected.subscribe(async () => {
      addBenchlingView('Mixtures', 'BenchlingLink:getMixtures');
    });

    const dnaOligosNode = treeNode.item('DNA Oligos');
    dnaOligosNode.onSelected.subscribe(async () => {
      addBenchlingView('DNA Oligos', 'BenchlingLink:getDnaOligos');
    });

    const projectsNode = treeNode.item('Projects');
    projectsNode.onSelected.subscribe(async () => {
      addBenchlingView('Projects', 'BenchlingLink:getProjects');
    });
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
//input: string names_anyOf {nullable:true}
//input: string names_anyOf_caseSensitive {nullable:true}
//input: string creatorIds {nullable:true}
//input: string authorIds_anyOf {nullable:true}
//input: string returning {nullable:true}
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
    names_anyOf,
    names_anyOf_caseSensitive,
    creatorIds,
    authorIds_anyOf,
    returning,
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
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
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
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
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
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
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
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
}

//name: Get Molecules
//input: string sort {nullable:true}
//input: string createdAt {nullable:true}
//input: string modifiedAt {nullable:true}
//input: string name {nullable:true}
//input: string nameIncludes {nullable:true}
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
//input: string authorIds_anyOf {nullable:true}
//input: string chemicalSubstructure_mol {nullable:true}
//input: string chemicalSubstructure_smiles {nullable:true}
//output: dataframe df
export async function getMolecules(
  sort?: string,
  createdAt?: string,
  modifiedAt?: string,
  name?: string,
  nameIncludes?: string,
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
  authorIds_anyOf?: string,
  chemicalSubstructure_mol?: string,
  chemicalSubstructure_smiles?: string,
): Promise<DG.DataFrame> {
  const params = {
    sort,
    createdAt,
    modifiedAt,
    name,
    nameIncludes,
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
    authorIds_anyOf,
    chemicalSubstructure_mol,
    chemicalSubstructure_smiles,
  };
  return await queryMolecules(params);
}

//name: Create Molecule
//input: string name
//input: string smiles
//input: string formula {nullable:true}
//output: dataframe df
export async function createMolecule(
  name: string,
  smiles: string,
  formula?: string,
): Promise<DG.DataFrame> {
  const body: any = { name, smiles };
  if (formula) body.formula = formula;
  const result = await postMolecule(body);
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
}

//name: Get Projects
//input: string sort {nullable:true}
//input: string archiveReason {nullable:true}
//input: string ids {nullable:true}
//input: string name {nullable:true}
//output: dataframe df
export async function getProjects(
  sort?: string,
  archiveReason?: string,
  ids?: string,
  name?: string,
): Promise<DG.DataFrame> {
  const params = {
    sort,
    archiveReason,
    ids,
    name,
  };
  const projects = await (await import('./projectsApi')).queryProjects(params);
  return projects;
}

//name: Get Plates
//input: string sort {nullable:true}
//input: string schemaId {nullable:true}
//input: string schemaFields {nullable:true}
//input: string createdAt {nullable:true}
//input: string modifiedAt {nullable:true}
//input: string name {nullable:true}
//input: string nameIncludes {nullable:true}
//input: int emptyPositions {nullable:true}
//input: int emptyPositions_gte {nullable:true}
//input: int emptyPositions_gt {nullable:true}
//input: int emptyPositions_lte {nullable:true}
//input: int emptyPositions_lt {nullable:true}
//input: int emptyContainers {nullable:true}
//input: int emptyContainers_gte {nullable:true}
//input: int emptyContainers_gt {nullable:true}
//input: int emptyContainers_lte {nullable:true}
//input: int emptyContainers_lt {nullable:true}
//input: string ancestorStorageId {nullable:true}
//input: string storageContentsId {nullable:true}
//input: string storageContentsIds {nullable:true}
//input: string archiveReason {nullable:true}
//input: string ids {nullable:true}
//input: string barcodes {nullable:true}
//input: string names_anyOf {nullable:true}
//input: string names_anyOf_caseSensitive {nullable:true}
//input: string returning {nullable:true}
//input: string creatorIds {nullable:true}
//output: dataframe df
export async function getPlates(
  sort?: string,
  schemaId?: string,
  schemaFields?: string,
  createdAt?: string,
  modifiedAt?: string,
  name?: string,
  nameIncludes?: string,
  emptyPositions?: number,
  emptyPositions_gte?: number,
  emptyPositions_gt?: number,
  emptyPositions_lte?: number,
  emptyPositions_lt?: number,
  emptyContainers?: number,
  emptyContainers_gte?: number,
  emptyContainers_gt?: number,
  emptyContainers_lte?: number,
  emptyContainers_lt?: number,
  ancestorStorageId?: string,
  storageContentsId?: string,
  storageContentsIds?: string,
  archiveReason?: string,
  ids?: string,
  barcodes?: string,
  names_anyOf?: string,
  names_anyOf_caseSensitive?: string,
  returning?: string,
  creatorIds?: string,
): Promise<DG.DataFrame> {
  const params = {
    sort,
    schemaId,
    schemaFields,
    createdAt,
    modifiedAt,
    name,
    nameIncludes,
    emptyPositions,
    emptyPositions_gte,
    emptyPositions_gt,
    emptyPositions_lte,
    emptyPositions_lt,
    emptyContainers,
    emptyContainers_gte,
    emptyContainers_gt,
    emptyContainers_lte,
    emptyContainers_lt,
    ancestorStorageId,
    storageContentsId,
    storageContentsIds,
    archiveReason,
    ids,
    barcodes,
    names_anyOf,
    names_anyOf_caseSensitive,
    returning,
    creatorIds,
  };
  return await queryPlates(params);
}

//name: Create Plate
//input: string name
//input: string schemaId
//input: string barcode {nullable:true}
//input: string containerSchemaId {nullable:true}
//input: string fields {nullable:true}
//input: string parentStorageId {nullable:true}
//input: string projectId {nullable:true}
//input: string wells {nullable:true}
//output: dataframe df
export async function createPlate(
  name: string,
  schemaId: string,
  barcode?: string,
  containerSchemaId?: string,
  fields?: string,
  parentStorageId?: string,
  projectId?: string,
  wells?: string,
): Promise<DG.DataFrame> {
  const body: any = { name, schemaId };
  if (barcode) body.barcode = barcode;
  if (containerSchemaId) body.containerSchemaId = containerSchemaId;
  if (fields) body.fields = JSON.parse(fields);
  if (parentStorageId) body.parentStorageId = parentStorageId;
  if (projectId) body.projectId = projectId;
  if (wells) body.wells = JSON.parse(wells);
  const result = await postPlate(body);
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
}

//name: Get Mixtures
//input: string sort {nullable:true}
//input: string createdAt {nullable:true}
//input: string modifiedAt {nullable:true}
//input: string name {nullable:true}
//input: string nameIncludes {nullable:true}
//input: string folderId {nullable:true}
//input: string mentionedIn {nullable:true}
//input: string projectId {nullable:true}
//input: string registryId {nullable:true}
//input: string schemaId {nullable:true}
//input: string schemaFields {nullable:true}
//input: string archiveReason {nullable:true}
//input: string mentions {nullable:true}
//input: string ids {nullable:true}
//input: string names_anyOf {nullable:true}
//input: string names_anyOf_caseSensitive {nullable:true}
//input: string entityRegistryIds_anyOf {nullable:true}
//input: string ingredientComponentEntityIds {nullable:true}
//input: string ingredientComponentEntityIds_anyOf {nullable:true}
//input: string authorIds_anyOf {nullable:true}
//output: dataframe df
export async function getMixtures(
  sort?: string,
  createdAt?: string,
  modifiedAt?: string,
  name?: string,
  nameIncludes?: string,
  folderId?: string,
  mentionedIn?: string,
  projectId?: string,
  registryId?: string,
  schemaId?: string,
  schemaFields?: string,
  archiveReason?: string,
  mentions?: string,
  ids?: string,
  names_anyOf?: string,
  names_anyOf_caseSensitive?: string,
  entityRegistryIds_anyOf?: string,
  ingredientComponentEntityIds?: string,
  ingredientComponentEntityIds_anyOf?: string,
  authorIds_anyOf?: string,
): Promise<DG.DataFrame> {
  const params = {
    sort,
    createdAt,
    modifiedAt,
    name,
    nameIncludes,
    folderId,
    mentionedIn,
    projectId,
    registryId,
    schemaId,
    schemaFields,
    archiveReason,
    mentions,
    ids,
    names_anyOf,
    names_anyOf_caseSensitive,
    entityRegistryIds_anyOf,
    ingredientComponentEntityIds,
    ingredientComponentEntityIds_anyOf,
    authorIds_anyOf,
  };
  return await queryMixtures(params);
}

//name: Create Mixture
//input: string name
//input: string ingredients
//input: string schemaId
//input: string units
//input: string aliases {nullable:true}
//input: string amount {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//input: string entityRegistryId {nullable:true}
//input: string fields {nullable:true}
//input: string folderId {nullable:true}
//output: dataframe df
export async function createMixture(
  name: string,
  ingredients: string,
  schemaId: string,
  units: string,
  aliases?: string,
  amount?: string,
  authorIds?: string,
  customFields?: string,
  entityRegistryId?: string,
  fields?: string,
  folderId?: string,
): Promise<DG.DataFrame> {
  const body: any = { name, schemaId, units };
  body.ingredients = JSON.parse(ingredients);
  if (aliases) body.aliases = JSON.parse(aliases);
  if (amount) body.amount = amount;
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  if (entityRegistryId) body.entityRegistryId = entityRegistryId;
  if (fields) body.fields = JSON.parse(fields);
  if (folderId) body.folderId = folderId;
  const result = await postMixture(body);
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
}

//name: Get DNA Oligos
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
//input: string names_anyOf {nullable:true}
//input: string names_anyOf_caseSensitive {nullable:true}
//input: string creatorIds {nullable:true}
//input: string authorIds_anyOf {nullable:true}
//input: string returning {nullable:true}
//input: string customNotationId {nullable:true}
//output: dataframe df
export async function getDnaOligos(
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
  entityRegistryIds_anyOf?: string,
  names_anyOf?: string,
  names_anyOf_caseSensitive?: string,
  creatorIds?: string,
  authorIds_anyOf?: string,
  returning?: string,
  customNotationId?: string,
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
    names_anyOf,
    names_anyOf_caseSensitive,
    creatorIds,
    authorIds_anyOf,
    returning,
    customNotationId,
  };
  return await queryDnaOligos(params);
}

//name: Create DNA Oligo
//input: string name
//input: string bases
//input: string aliases {nullable:true}
//input: string annotations {nullable:true}
//input: string authorIds {nullable:true}
//input: string customFields {nullable:true}
//input: string fields {nullable:true}
//input: string folderId {nullable:true}
//input: string schemaId {nullable:true}
//input: string helm {nullable:true}
//output: dataframe df
export async function createDnaOligo(
  name: string,
  bases: string,
  aliases?: string,
  annotations?: string,
  authorIds?: string,
  customFields?: string,
  fields?: string,
  folderId?: string,
  schemaId?: string,
  helm?: string,
): Promise<DG.DataFrame> {
  const body: any = { name, bases };
  if (aliases) body.aliases = JSON.parse(aliases);
  if (annotations) body.annotations = JSON.parse(annotations);
  if (authorIds) body.authorIds = JSON.parse(authorIds);
  if (customFields) body.customFields = JSON.parse(customFields);
  if (fields) body.fields = JSON.parse(fields);
  if (folderId) body.folderId = folderId;
  if (schemaId) body.schemaId = schemaId;
  if (helm) body.helm = helm;
  const result = await postDnaOligo(body);
  return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
}