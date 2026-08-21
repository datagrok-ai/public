/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {AlignmentHandler, buildTreeBrowser, routeView} from './app/catalog-app';
import {showPublishDialog} from './app/publish-dialog';
import {ensureProgram, ProgramSpec, setupSchemaSecurity} from './domain/security';
import {
  approvePublication, discover, getPublicationHistory, publishWorkflowRun as publishRun,
  PublishRequest, rejectPublication, updateCuration,
} from './service/publication-service';

export * from './package.g';

export const _package = new DG.Package();

//tags: autostart
//meta.autostartImmediate: true
export function _initArtifactAlignment(): void {
  DG.ObjectHandler.register(new AlignmentHandler());
}

//name: Artifact Catalog
//description: Program-aligned catalog of published Compute2 workflow runs
//tags: app
//meta.icon: images/catalog.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function artifactCatalogApp(path?: string): Promise<DG.ViewBase> {
  return await routeView(path);
}

//name: artifactCatalogTreeBrowser
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: Artifact Catalog
export async function artifactCatalogTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  await buildTreeBrowser(treeNode);
}

//name: publishWorkflowRun
//description: Publishes a saved Compute2 workflow run into a program as a frozen copy
//input: object request
//input: funccall sourceCall {optional: true}
//output: object result
export async function publishWorkflowRun(request: PublishRequest, sourceCall?: DG.FuncCall): Promise<object> {
  return publishRun(sourceCall != null ? {...request, sourceCall} : request);
}

//name: publishWorkflowRunDialog
//description: Opens the "Publish to program" dialog for a Compute2 run — a saved run id or a live FuncCall
//input: string sourceMetaCallId {optional: true}
//input: funccall sourceCall {optional: true}
//input: string defaultName {optional: true}
export async function publishWorkflowRunDialog(
  sourceMetaCallId?: string, sourceCall?: DG.FuncCall, defaultName?: string): Promise<void> {
  const source = sourceCall ?? sourceMetaCallId;
  if (source == null)
    throw new Error('Either sourceCall or sourceMetaCallId is required');
  await showPublishDialog(source, defaultName);
}

//name: approveArtifactPublication
//input: string rowId
export async function approveArtifactPublication(rowId: string): Promise<void> {
  await approvePublication(rowId);
}

//name: rejectArtifactPublication
//input: string rowId
//input: string reason
export async function rejectArtifactPublication(rowId: string, reason: string): Promise<void> {
  await rejectPublication(rowId, reason);
}

//name: updateArtifactCuration
//input: string rowId
//input: object curation
export async function updateArtifactCuration(rowId: string, curation: {tags?: string[], path?: string}): Promise<void> {
  await updateCuration(rowId, curation);
}

//name: discoverArtifacts
//description: The discovery slice: approved rows only, one per publication
//input: string programId {optional: true}
//output: object result
export async function discoverArtifacts(programId?: string): Promise<object> {
  return discover(programId ? programId : undefined);
}

//name: getArtifactHistory
//input: string publicationId
//output: object result
export async function getArtifactHistory(publicationId: string): Promise<object> {
  return getPublicationHistory(publicationId);
}

//name: ensureArtifactProgram
//description: Idempotently provisions a program with its viewers/contributors/approvers groups and row grants
//input: object spec
//output: object result
export async function ensureArtifactProgram(spec: ProgramSpec): Promise<object> {
  const info = await ensureProgram(spec);
  return {id: info.id, code: info.code, viewers: info.viewers.id,
    contributors: info.contributors.id, approvers: info.approvers.id};
}

//name: setupArtifactSecurity
//description: One-time idempotent security setup: umbrella groups, registry visibility, per-column write gating
export async function setupArtifactSecurity(): Promise<void> {
  await setupSchemaSecurity();
}

//name: AaTestStep
//description: Test fixture: a workflow step with a dataframe output
//input: int a = 3
//output: dataframe result
export function aaTestStep(a: number): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('int', 'x', Array.from({length: a}, (_, i) => i)),
  ]);
}

//name: AaTestProvider
//description: Test fixture: stands in for a pipeline provider function
//output: object result
export function aaTestProvider(): object {
  return {};
}
