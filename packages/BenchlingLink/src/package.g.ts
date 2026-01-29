import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Benchling
//output: view result
//meta.role: app
//meta.browsePath: Chem
export async function benchlingLinkApp() : Promise<any> {
  return await PackageFunctions.benchlingLinkApp();
}

//input: dynamic treeNode 
export async function benchlingLinkAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.benchlingLinkAppTreeBrowser(treeNode);
}

//name: Get AA Sequences
//input: string sort { nullable: true }
//input: string createdAt { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string name { nullable: true }
//input: string nameIncludes { nullable: true }
//input: string aminoAcids { nullable: true }
//input: string folderId { nullable: true }
//input: string mentionedIn { nullable: true }
//input: string projectId { nullable: true }
//input: string registryId { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string archiveReason { nullable: true }
//input: string mentions { nullable: true }
//input: string ids { nullable: true }
//input: string entityRegistryIds_anyOf { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string names_anyOf_caseSensitive { nullable: true }
//input: string creatorIds { nullable: true }
//input: string authorIds_anyOf { nullable: true }
//input: string returning { nullable: true }
//output: dataframe result
export async function getAASequences(sort?: string, createdAt?: string, modifiedAt?: string, name?: string, nameIncludes?: string, aminoAcids?: string, folderId?: string, mentionedIn?: string, projectId?: string, registryId?: string, schemaId?: string, schemaFields?: string, archiveReason?: string, mentions?: string, ids?: string, entityRegistryIds_anyOf?: string, names_anyOf?: string, names_anyOf_caseSensitive?: string, creatorIds?: string, authorIds_anyOf?: string, returning?: string) : Promise<any> {
  return await PackageFunctions.getAASequences(sort, createdAt, modifiedAt, name, nameIncludes, aminoAcids, folderId, mentionedIn, projectId, registryId, schemaId, schemaFields, archiveReason, mentions, ids, entityRegistryIds_anyOf, names_anyOf, names_anyOf_caseSensitive, creatorIds, authorIds_anyOf, returning);
}

//name: Get DNA Sequences
//input: string sort { nullable: true }
//input: string createdAt { nullable: true }
//input: string name { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string nameIncludes { nullable: true }
//input: string bases { nullable: true }
//input: string folderId { nullable: true }
//input: string mentionedIn { nullable: true }
//input: string projectId { nullable: true }
//input: string registryId { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string archiveReason { nullable: true }
//input: string mentions { nullable: true }
//input: string ids { nullable: true }
//input: string entityRegistryIds_anyOf { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string names_anyOf_caseSensitive { nullable: true }
//input: string creatorIds { nullable: true }
//input: string authorIds_anyOf { nullable: true }
//input: string returning { nullable: true }
//output: dataframe result
export async function getDNASequences(sort?: string, createdAt?: string, name?: string, modifiedAt?: string, nameIncludes?: string, bases?: string, folderId?: string, mentionedIn?: string, projectId?: string, registryId?: string, schemaId?: string, schemaFields?: string, archiveReason?: string, mentions?: string, ids?: string, entityRegistryIds_anyOf?: string, names_anyOf?: string, names_anyOf_caseSensitive?: string, creatorIds?: string, authorIds_anyOf?: string, returning?: string) : Promise<any> {
  return await PackageFunctions.getDNASequences(sort, createdAt, name, modifiedAt, nameIncludes, bases, folderId, mentionedIn, projectId, registryId, schemaId, schemaFields, archiveReason, mentions, ids, entityRegistryIds_anyOf, names_anyOf, names_anyOf_caseSensitive, creatorIds, authorIds_anyOf, returning);
}

//name: Get Assay Results
//input: string schemaId { nullable: true }
//input: string createdAt_lt { nullable: true }
//input: string createdAt_gt { nullable: true }
//input: string createdAt_lte { nullable: true }
//input: string createdAt_gte { nullable: true }
//input: double minCreatedTime { nullable: true }
//input: double maxCreatedTime { nullable: true }
//input: string sort { nullable: true }
//input: string entityIds { nullable: true }
//input: string storageIds { nullable: true }
//input: string assayRunIds { nullable: true }
//input: string automationOutputProcessorId { nullable: true }
//input: string ids { nullable: true }
//input: string modifiedAt_lt { nullable: true }
//input: string modifiedAt_gt { nullable: true }
//input: string modifiedAt_lte { nullable: true }
//input: string modifiedAt_gte { nullable: true }
//input: string archiveReason { nullable: true }
//output: dataframe result
export async function getAssayResults(schemaId?: string, createdAt_lt?: string, createdAt_gt?: string, createdAt_lte?: string, createdAt_gte?: string, minCreatedTime?: number, maxCreatedTime?: number, sort?: string, entityIds?: string, storageIds?: string, assayRunIds?: string, automationOutputProcessorId?: string, ids?: string, modifiedAt_lt?: string, modifiedAt_gt?: string, modifiedAt_lte?: string, modifiedAt_gte?: string, archiveReason?: string) : Promise<any> {
  return await PackageFunctions.getAssayResults(schemaId, createdAt_lt, createdAt_gt, createdAt_lte, createdAt_gte, minCreatedTime, maxCreatedTime, sort, entityIds, storageIds, assayRunIds, automationOutputProcessorId, ids, modifiedAt_lt, modifiedAt_gt, modifiedAt_lte, modifiedAt_gte, archiveReason);
}

//name: Get Assay Runs
//input: string schemaId { nullable: true }
//input: double minCreatedTime { nullable: true }
//input: double maxCreatedTime { nullable: true }
//input: string ids { nullable: true }
//output: dataframe result
export async function getAssayRuns(schemaId?: string, minCreatedTime?: number, maxCreatedTime?: number, ids?: string) : Promise<any> {
  return await PackageFunctions.getAssayRuns(schemaId, minCreatedTime, maxCreatedTime, ids);
}

//name: Create AA Sequence
//input: string name 
//input: string aminoAcids 
//input: string aliases { nullable: true }
//input: string annotations { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//input: string fields { nullable: true }
//input: string folderId { nullable: true }
//input: string schemaId { nullable: true }
//output: dataframe result
export async function createAASequence(name: string, aminoAcids: string, aliases?: string, annotations?: string, authorIds?: string, customFields?: string, fields?: string, folderId?: string, schemaId?: string) : Promise<any> {
  return await PackageFunctions.createAASequence(name, aminoAcids, aliases, annotations, authorIds, customFields, fields, folderId, schemaId);
}

//name: Create DNA Sequence
//input: string name 
//input: string bases 
//input: string aliases { nullable: true }
//input: string annotations { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//input: string fields { nullable: true }
//input: string folderId { nullable: true }
//input: string schemaId { nullable: true }
//output: dataframe result
export async function createDNASequence(name: string, bases: string, aliases?: string, annotations?: string, authorIds?: string, customFields?: string, fields?: string, folderId?: string, schemaId?: string) : Promise<any> {
  return await PackageFunctions.createDNASequence(name, bases, aliases, annotations, authorIds, customFields, fields, folderId, schemaId);
}

//name: Create Assay Result
//input: string schemaId 
//input: string fields { nullable: true }
//input: string entityIds { nullable: true }
//input: string storageIds { nullable: true }
//input: string assayRunId { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//output: dataframe result
export async function createAssayResult(schemaId: string, fields?: string, entityIds?: string, storageIds?: string, assayRunId?: string, authorIds?: string, customFields?: string) : Promise<any> {
  return await PackageFunctions.createAssayResult(schemaId, fields, entityIds, storageIds, assayRunId, authorIds, customFields);
}

//name: Create Assay Run
//input: string schemaId 
//input: string fields { nullable: true }
//input: string name { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//output: dataframe result
export async function createAssayRun(schemaId: string, fields?: string, name?: string, authorIds?: string, customFields?: string) : Promise<any> {
  return await PackageFunctions.createAssayRun(schemaId, fields, name, authorIds, customFields);
}

//name: Get Molecules
//input: string sort { nullable: true }
//input: string createdAt { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string name { nullable: true }
//input: string nameIncludes { nullable: true }
//input: string folderId { nullable: true }
//input: string mentionedIn { nullable: true }
//input: string projectId { nullable: true }
//input: string registryId { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string archiveReason { nullable: true }
//input: string mentions { nullable: true }
//input: string ids { nullable: true }
//input: string entityRegistryIds_anyOf { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string authorIds_anyOf { nullable: true }
//input: string chemicalSubstructure_mol { nullable: true }
//input: string chemicalSubstructure_smiles { nullable: true }
//output: dataframe result
export async function getMolecules(sort?: string, createdAt?: string, modifiedAt?: string, name?: string, nameIncludes?: string, folderId?: string, mentionedIn?: string, projectId?: string, registryId?: string, schemaId?: string, schemaFields?: string, archiveReason?: string, mentions?: string, ids?: string, entityRegistryIds_anyOf?: string, names_anyOf?: string, authorIds_anyOf?: string, chemicalSubstructure_mol?: string, chemicalSubstructure_smiles?: string) : Promise<any> {
  return await PackageFunctions.getMolecules(sort, createdAt, modifiedAt, name, nameIncludes, folderId, mentionedIn, projectId, registryId, schemaId, schemaFields, archiveReason, mentions, ids, entityRegistryIds_anyOf, names_anyOf, authorIds_anyOf, chemicalSubstructure_mol, chemicalSubstructure_smiles);
}

//name: Create Molecule
//input: string name 
//input: string smiles 
//input: string formula { nullable: true }
//output: dataframe result
export async function createMolecule(name: string, smiles: string, formula?: string) : Promise<any> {
  return await PackageFunctions.createMolecule(name, smiles, formula);
}

//name: Get Projects
//input: string sort { nullable: true }
//input: string archiveReason { nullable: true }
//input: string ids { nullable: true }
//input: string name { nullable: true }
//output: dataframe result
export async function getProjects(sort?: string, archiveReason?: string, ids?: string, name?: string) : Promise<any> {
  return await PackageFunctions.getProjects(sort, archiveReason, ids, name);
}

//name: Get Plates
//input: string sort { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string createdAt { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string name { nullable: true }
//input: string nameIncludes { nullable: true }
//input: double emptyPositions { nullable: true }
//input: double emptyPositions_gte { nullable: true }
//input: double emptyPositions_gt { nullable: true }
//input: double emptyPositions_lte { nullable: true }
//input: double emptyPositions_lt { nullable: true }
//input: double emptyContainers { nullable: true }
//input: double emptyContainers_gte { nullable: true }
//input: double emptyContainers_gt { nullable: true }
//input: double emptyContainers_lte { nullable: true }
//input: double emptyContainers_lt { nullable: true }
//input: string ancestorStorageId { nullable: true }
//input: string storageContentsId { nullable: true }
//input: string storageContentsIds { nullable: true }
//input: string archiveReason { nullable: true }
//input: string ids { nullable: true }
//input: string barcodes { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string names_anyOf_caseSensitive { nullable: true }
//input: string returning { nullable: true }
//input: string creatorIds { nullable: true }
//output: dataframe result
export async function getPlates(sort?: string, schemaId?: string, schemaFields?: string, createdAt?: string, modifiedAt?: string, name?: string, nameIncludes?: string, emptyPositions?: number, emptyPositions_gte?: number, emptyPositions_gt?: number, emptyPositions_lte?: number, emptyPositions_lt?: number, emptyContainers?: number, emptyContainers_gte?: number, emptyContainers_gt?: number, emptyContainers_lte?: number, emptyContainers_lt?: number, ancestorStorageId?: string, storageContentsId?: string, storageContentsIds?: string, archiveReason?: string, ids?: string, barcodes?: string, names_anyOf?: string, names_anyOf_caseSensitive?: string, returning?: string, creatorIds?: string) : Promise<any> {
  return await PackageFunctions.getPlates(sort, schemaId, schemaFields, createdAt, modifiedAt, name, nameIncludes, emptyPositions, emptyPositions_gte, emptyPositions_gt, emptyPositions_lte, emptyPositions_lt, emptyContainers, emptyContainers_gte, emptyContainers_gt, emptyContainers_lte, emptyContainers_lt, ancestorStorageId, storageContentsId, storageContentsIds, archiveReason, ids, barcodes, names_anyOf, names_anyOf_caseSensitive, returning, creatorIds);
}

//name: Create Plate
//input: string name 
//input: string schemaId 
//input: string barcode { nullable: true }
//input: string containerSchemaId { nullable: true }
//input: string fields { nullable: true }
//input: string parentStorageId { nullable: true }
//input: string projectId { nullable: true }
//input: string wells { nullable: true }
//output: dataframe result
export async function createPlate(name: string, schemaId: string, barcode?: string, containerSchemaId?: string, fields?: string, parentStorageId?: string, projectId?: string, wells?: string) : Promise<any> {
  return await PackageFunctions.createPlate(name, schemaId, barcode, containerSchemaId, fields, parentStorageId, projectId, wells);
}

//name: Get Mixtures
//input: string sort { nullable: true }
//input: string createdAt { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string name { nullable: true }
//input: string nameIncludes { nullable: true }
//input: string folderId { nullable: true }
//input: string mentionedIn { nullable: true }
//input: string projectId { nullable: true }
//input: string registryId { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string archiveReason { nullable: true }
//input: string mentions { nullable: true }
//input: string ids { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string names_anyOf_caseSensitive { nullable: true }
//input: string entityRegistryIds_anyOf { nullable: true }
//input: string ingredientComponentEntityIds { nullable: true }
//input: string ingredientComponentEntityIds_anyOf { nullable: true }
//input: string authorIds_anyOf { nullable: true }
//output: dataframe result
export async function getMixtures(sort?: string, createdAt?: string, modifiedAt?: string, name?: string, nameIncludes?: string, folderId?: string, mentionedIn?: string, projectId?: string, registryId?: string, schemaId?: string, schemaFields?: string, archiveReason?: string, mentions?: string, ids?: string, names_anyOf?: string, names_anyOf_caseSensitive?: string, entityRegistryIds_anyOf?: string, ingredientComponentEntityIds?: string, ingredientComponentEntityIds_anyOf?: string, authorIds_anyOf?: string) : Promise<any> {
  return await PackageFunctions.getMixtures(sort, createdAt, modifiedAt, name, nameIncludes, folderId, mentionedIn, projectId, registryId, schemaId, schemaFields, archiveReason, mentions, ids, names_anyOf, names_anyOf_caseSensitive, entityRegistryIds_anyOf, ingredientComponentEntityIds, ingredientComponentEntityIds_anyOf, authorIds_anyOf);
}

//name: Create Mixture
//input: string name 
//input: string ingredients 
//input: string schemaId 
//input: string units 
//input: string aliases { nullable: true }
//input: string amount { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//input: string entityRegistryId { nullable: true }
//input: string fields { nullable: true }
//input: string folderId { nullable: true }
//output: dataframe result
export async function createMixture(name: string, ingredients: string, schemaId: string, units: string, aliases?: string, amount?: string, authorIds?: string, customFields?: string, entityRegistryId?: string, fields?: string, folderId?: string) : Promise<any> {
  return await PackageFunctions.createMixture(name, ingredients, schemaId, units, aliases, amount, authorIds, customFields, entityRegistryId, fields, folderId);
}

//name: Get DNA Oligos
//input: string sort { nullable: true }
//input: string createdAt { nullable: true }
//input: string modifiedAt { nullable: true }
//input: string name { nullable: true }
//input: string nameIncludes { nullable: true }
//input: string bases { nullable: true }
//input: string folderId { nullable: true }
//input: string mentionedIn { nullable: true }
//input: string projectId { nullable: true }
//input: string registryId { nullable: true }
//input: string schemaId { nullable: true }
//input: string schemaFields { nullable: true }
//input: string archiveReason { nullable: true }
//input: string mentions { nullable: true }
//input: string ids { nullable: true }
//input: string entityRegistryIds_anyOf { nullable: true }
//input: string names_anyOf { nullable: true }
//input: string names_anyOf_caseSensitive { nullable: true }
//input: string creatorIds { nullable: true }
//input: string authorIds_anyOf { nullable: true }
//input: string returning { nullable: true }
//input: string customNotationId { nullable: true }
//output: dataframe result
export async function getDnaOligos(sort?: string, createdAt?: string, modifiedAt?: string, name?: string, nameIncludes?: string, bases?: string, folderId?: string, mentionedIn?: string, projectId?: string, registryId?: string, schemaId?: string, schemaFields?: string, archiveReason?: string, mentions?: string, ids?: string, entityRegistryIds_anyOf?: string, names_anyOf?: string, names_anyOf_caseSensitive?: string, creatorIds?: string, authorIds_anyOf?: string, returning?: string, customNotationId?: string) : Promise<any> {
  return await PackageFunctions.getDnaOligos(sort, createdAt, modifiedAt, name, nameIncludes, bases, folderId, mentionedIn, projectId, registryId, schemaId, schemaFields, archiveReason, mentions, ids, entityRegistryIds_anyOf, names_anyOf, names_anyOf_caseSensitive, creatorIds, authorIds_anyOf, returning, customNotationId);
}

//name: Create DNA Oligo
//input: string name 
//input: string bases 
//input: string aliases { nullable: true }
//input: string annotations { nullable: true }
//input: string authorIds { nullable: true }
//input: string customFields { nullable: true }
//input: string fields { nullable: true }
//input: string folderId { nullable: true }
//input: string schemaId { nullable: true }
//input: string helm { nullable: true }
//output: dataframe result
export async function createDnaOligo(name: string, bases: string, aliases?: string, annotations?: string, authorIds?: string, customFields?: string, fields?: string, folderId?: string, schemaId?: string, helm?: string) : Promise<any> {
  return await PackageFunctions.createDnaOligo(name, bases, aliases, annotations, authorIds, customFields, fields, folderId, schemaId, helm);
}
