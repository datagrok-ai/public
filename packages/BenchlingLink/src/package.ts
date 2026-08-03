/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { queryAASequences, postAASequence } from './aaSequencesApi';
import { queryDNASequences, postDNASequence } from './dnaSequencesApi';
import { queryAssayResults, queryAssayRuns, postAssayResults, postAssayRuns } from './assayApi';
import { queryMolecules, postMolecule } from './moleculesApi';
import { dataFrameFromObjects } from './utils';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { queryPlates, postPlate } from './platesApi';
import { queryMixtures, postMixture } from './mixturesApi';
import { queryDnaOligos, postDnaOligo } from './dnaOligosApi';
import { queryProjects } from './projectsApi';

export * from './package.g';
export const _package = new DG.Package();

// JSON-parse a string body field only when caller supplied it. Used by every createX wrapper
// because InputForm delivers array/object inputs as JSON strings.
function assignJson(body: any, key: string, raw: string | undefined): void {
  if (raw)
    body[key] = JSON.parse(raw);
}

export class PackageFunctions {
  @grok.decorators.app({
    browsePath: 'Chem',
    name: 'Benchling'
  })
  static async benchlingLinkApp(): Promise<DG.ViewBase> {

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

  @grok.decorators.appTreeBrowser({app: 'Benchling'})
  static async benchlingLinkAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
    const addFuncNode = (label: string, funcName: string): void => {
      const node = treeNode.item(label);
      node.value = DG.Func.byName(funcName);
      node.onSelected.subscribe(async () => {
         const objHandler = DG.ObjectHandler.forEntity(node.value);
        if (objHandler)
          grok.shell.preview = await objHandler.renderPreview(node.value);
      });
    };

    try {
      treeNode.items.length = 0;
      addFuncNode('AA Sequences', 'BenchlingLink:getAASequences');
      addFuncNode('DNA Sequences', 'BenchlingLink:getDNASequences');
      addFuncNode('Assay Results', 'BenchlingLink:getAssayResults');
      addFuncNode('Assay Runs', 'BenchlingLink:getAssayRuns');
      addFuncNode('Molecules', 'BenchlingLink:getMolecules');
      addFuncNode('Plates', 'BenchlingLink:getPlates');
      addFuncNode('Mixtures', 'BenchlingLink:getMixtures');
      addFuncNode('DNA Oligos', 'BenchlingLink:getDnaOligos');
      addFuncNode('Projects', 'BenchlingLink:getProjects');
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
    }
  }

  @grok.decorators.func({name: 'Get AA Sequences'})
  static async getAASequences(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) aminoAcids?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string
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
    return await queryAASequences(params);
  }

  @grok.decorators.func({name: 'Get DNA Sequences'})
  static async getDNASequences(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) bases?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string
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
    return await queryDNASequences(params);
  }

  @grok.decorators.func({name: 'Get Assay Results'})
  static async getAssayResults(
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_lt?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_gt?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_lte?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_gte?: string,
    @grok.decorators.param({options: {nullable: true}}) minCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) maxCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) entityIds?: string,
    @grok.decorators.param({options: {nullable: true}}) storageIds?: string,
    @grok.decorators.param({options: {nullable: true}}) assayRunIds?: string,
    @grok.decorators.param({options: {nullable: true}}) automationOutputProcessorId?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_lt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_gt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_lte?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_gte?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string
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
    return await queryAssayResults(params);
  }

  @grok.decorators.func({name: 'Get Assay Runs'})
  static async getAssayRuns(
    schemaId: string,
    @grok.decorators.param({options: {nullable: true}}) minCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) maxCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) ids?: string
  ): Promise<DG.DataFrame> {
    const params = {
      schemaId,
      minCreatedTime,
      maxCreatedTime,
      ids,
    };
    return await queryAssayRuns(params);
  }

  @grok.decorators.func({name: 'Create AA Sequence'})
  static async createAASequence(
    name: string,
    aminoAcids: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, aminoAcids };
    assignJson(body, 'aliases', aliases);
    assignJson(body, 'annotations', annotations);
    assignJson(body, 'authorIds', authorIds);
    assignJson(body, 'customFields', customFields);
    assignJson(body, 'fields', fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    const result = await postAASequence(body);
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Create DNA Sequence'})
  static async createDNASequence(
    name: string,
    bases: string,
    isCircular: boolean,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, bases, isCircular };
    assignJson(body, 'aliases', aliases);
    assignJson(body, 'annotations', annotations);
    assignJson(body, 'authorIds', authorIds);
    assignJson(body, 'customFields', customFields);
    assignJson(body, 'fields', fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    const result = await postDNASequence(body);
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Create Assay Result'})
  static async createAssayResult(
    schemaId: string,
    fields: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) fieldValidation?: string,
    @grok.decorators.param({options: {nullable: true}}) id?: string,
  ): Promise<DG.DataFrame> {
    const item: any = { schemaId, fields: JSON.parse(fields) };
    if (projectId) item.projectId = projectId;
    assignJson(item, 'fieldValidation', fieldValidation);
    if (id) item.id = id;
    const result = await postAssayResults({ assayResults: [item] });
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Create Assay Run'})
  static async createAssayRun(
    schemaId: string,
    fields: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) validationComment?: string,
    @grok.decorators.param({options: {nullable: true}}) validationStatus?: string,
    @grok.decorators.param({options: {nullable: true}}) id?: string,
  ): Promise<DG.DataFrame> {
    const item: any = { schemaId, fields: JSON.parse(fields) };
    if (projectId) item.projectId = projectId;
    if (validationComment) item.validationComment = validationComment;
    if (validationStatus) item.validationStatus = validationStatus;
    if (id) item.id = id;
    const result = await postAssayRuns({ assayRuns: [item] });
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Get Molecules'})
  static async getMolecules(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) chemicalSubstructure_mol?: string,
    @grok.decorators.param({options: {nullable: true}}) chemicalSubstructure_smiles?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      folderId,
      mentionedIn: mentionedIn ? mentionedIn.split(',').map((s) => s.trim()).filter(Boolean) : undefined,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions: mentions ? mentions.split(',').map((s) => s.trim()).filter(Boolean) : undefined,
      ids,
      entityRegistryIds_anyOf,
      names_anyOf,
      authorIds_anyOf,
      chemicalSubstructure_mol,
      chemicalSubstructure_smiles,
    };
    return await queryMolecules(params);
  }

  @grok.decorators.func({name: 'Create Molecule'})
  static async createMolecule(
    name: string,
    value: string,
    @grok.decorators.param({options: {nullable: true}}) structureFormat?: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
  ): Promise<DG.DataFrame> {
    const format = structureFormat ?? 'smiles';
    if (format !== 'smiles' && format !== 'molfile')
      throw new Error(`createMolecule: unknown structureFormat '${format}'. Expected 'smiles' or 'molfile'.`);
    const body: any = { name, chemicalStructure: { structureFormat: format, value } };
    assignJson(body, 'aliases', aliases);
    assignJson(body, 'authorIds', authorIds);
    assignJson(body, 'customFields', customFields);
    assignJson(body, 'fields', fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    const result = await postMolecule(body);
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Get Projects'})
  static async getProjects(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      archiveReason,
      ids,
      name,
    };
    return await queryProjects(params);
  }

  @grok.decorators.func({name: 'Get Plates'})
  static async getPlates(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_gte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_gt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_lte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_lt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_gte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_gt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_lte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_lt?: number,
    @grok.decorators.param({options: {nullable: true}}) ancestorStorageId?: string,
    @grok.decorators.param({options: {nullable: true}}) storageContentsId?: string,
    @grok.decorators.param({options: {nullable: true}}) storageContentsIds?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) barcodes?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
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

  @grok.decorators.func({name: 'Create Plate'})
  static async createPlate(
    name: string,
    schemaId: string,
    @grok.decorators.param({options: {nullable: true}}) barcode?: string,
    @grok.decorators.param({options: {nullable: true}}) containerSchemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) parentStorageId?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) wells?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, schemaId };
    if (barcode) body.barcode = barcode;
    if (containerSchemaId) body.containerSchemaId = containerSchemaId;
    assignJson(body, 'fields', fields);
    if (parentStorageId) body.parentStorageId = parentStorageId;
    if (projectId) body.projectId = projectId;
    assignJson(body, 'wells', wells);
    const result = await postPlate(body);
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Get Mixtures'})
  static async getMixtures(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) ingredientComponentEntityIds?: string,
    @grok.decorators.param({options: {nullable: true}}) ingredientComponentEntityIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
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

  @grok.decorators.func({name: 'Create Mixture'})
  static async createMixture(
    name: string,
    ingredients: string,
    schemaId: string,
    units: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) amount?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryId?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, schemaId, units, ingredients: JSON.parse(ingredients) };
    assignJson(body, 'aliases', aliases);
    if (amount) body.amount = amount;
    assignJson(body, 'authorIds', authorIds);
    assignJson(body, 'customFields', customFields);
    if (entityRegistryId) body.entityRegistryId = entityRegistryId;
    assignJson(body, 'fields', fields);
    if (folderId) body.folderId = folderId;
    const result = await postMixture(body);
    return dataFrameFromObjects([result]);
  }

  @grok.decorators.func({name: 'Get DNA Oligos'})
  static async getDnaOligos(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) bases?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string,
    @grok.decorators.param({options: {nullable: true}}) customNotationId?: string,
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

  @grok.decorators.func({name: 'Create DNA Oligo'})
  static async createDnaOligo(
    name: string,
    bases: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) helm?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, bases };
    assignJson(body, 'aliases', aliases);
    assignJson(body, 'annotations', annotations);
    assignJson(body, 'authorIds', authorIds);
    assignJson(body, 'customFields', customFields);
    assignJson(body, 'fields', fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    if (helm) body.helm = helm;
    const result = await postDnaOligo(body);
    return dataFrameFromObjects([result]);
  }
}
