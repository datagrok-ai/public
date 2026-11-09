/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {routeNcaApp} from './views/nca-studies-router';
import {buildNcaTree} from './views/nca-tree-browser';
import {computePipeline} from './pipeline/pipeline-compute';
import {detectPipelineMapping} from './pipeline/pipeline-mapping';
import {summarizePipeline} from './pipeline/pipeline-summary';
import {saveToStudies} from './pipeline/pipeline-save';
import {seedDemoStudies, productionSeedDeps} from './pipeline/demo-seed';
import {DEMO_STUDIES_DIR} from './adapter/demo-studies';
import {loadSampleByName} from './views/nca-launch';
import {detectPowerGridSupport} from './views/powergrid-support';
import {registerNcaRepositoryFactory, EmsNcaRepository} from './adapter';

export const _package = new DG.Package();

// The session's storage (nca-ems-storage-migration): the EMS schema `nca`. A factory, not an
// instance — the repository is built on first use, after the platform is up; tests swap in an
// in-memory one with `setNcaRepository`.
registerNcaRepositoryFactory(() => new EmsNcaRepository());

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Non-Compartmental Analysis
//description: Non-Compartmental Analysis for preclinical PK studies
//tags: app
//meta.browsePath: Compute
//input: string path {meta.url: true; optional: true}
//output: view result
export async function ncaApp(path?: string): Promise<DG.ViewBase> {
  // The PowerGrid check runs at app start (user-directed 2026-10-19) and is cached for the session;
  // its answer is the default of the Settings `Violated column` switch.
  detectPowerGridSupport();
  // metaUrl routing (AD-1): `path` is the trailing segment(s) after /apps/Nca.
  // The router parses it into landing / project-overview / analysis and handles
  // the path-aware hub-dedup (AC-A5): empty path → return the open hub;
  // non-empty path → route first to the project/analysis (no hub interception).
  return routeNcaApp(path ?? '');
}

//input: dynamic treeNode
//meta.app: Non-Compartmental Analysis
export async function ncaAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  // Browse-tree builder, paired to the `Non-Compartmental Analysis` app by the
  // `//meta.app` binding above (datagrok-app-tree-pattern.md §"The two halves").
  // Studies → groups, analyses → leaves, reactive on NCA_STUDIES_CHANGED_EVENT.
  detectPowerGridSupport(); // an analysis opened from the tree needs the same default
  await buildNcaTree(treeNode);
}

// ---- Pipeline functions (headless, Flow-ready) ----

//name: ncaDetectMapping
//description: Detects the NCA column mapping of a table. Returns JSON with mapping, confidence and issues.
//meta.includeInFlow: true
//input: dataframe table {nullable: false}
//output: string mapping
export function ncaDetectMapping(table: DG.DataFrame): string {
  return JSON.stringify(detectPipelineMapping(table));
}

//name: ncaCompute
//description: Headless NCA: per-profile verdict and reasons; QC flags are site-tunable screening heuristics
//meta.includeInFlow: true
//input: dataframe table {nullable: false}
//input: string mapping = "" {optional: true; description: Mapping JSON or ncaDetectMapping output. Empty = auto-detect}
//input: string rules = "" {optional: true; description: NcaRules JSON plus pipeline policy keys. Empty = app defaults}
//output: dataframe results
export function ncaCompute(table: DG.DataFrame, mapping?: string, rules?: string): DG.DataFrame {
  return computePipeline(table, mapping ?? '', rules ?? '');
}

//name: ncaSummarize
//description: Per-arm summary of ncaCompute results - N, mean, CV%, profiles to review, F% vs the same-sex IV arm.
//meta.includeInFlow: true
//input: dataframe results {nullable: false}
//output: dataframe summary
export function ncaSummarize(results: DG.DataFrame): DG.DataFrame {
  return summarizePipeline(results);
}

//name: ncaSaveToStudy
//description: Saves a pipeline run into NCA Studio. A re-run updates untouched analyses, saves beside reviewed ones
//meta.includeInFlow: true
//input: dataframe table {nullable: false; description: The observations table given to ncaCompute}
//input: dataframe results {nullable: false; description: The output of ncaCompute}
//input: string share = "" {optional: true; description: JSON with view and edit group names for a NEW study}
//input: string owner = "" {optional: true; description: Study owner key - set it stably for scheduled runs}
//output: dataframe saved
export async function ncaSaveToStudy(table: DG.DataFrame, results: DG.DataFrame,
  share?: string, owner?: string): Promise<DG.DataFrame> {
  return saveToStudies(table, results, share ?? '', owner ?? '');
}

/* eslint-disable max-len -- a function header line cannot wrap; the seeding requirement is stated in full here (nca-hub-demo-studies F7) */
//name: ncaSeedDemoStudies
//description: Makes sure the five simulated NCA demo studies are installed and current, shared View to All users. Any user can run it: the first caller installs the demo and owns it. Re-run after a sci-comp update: the owner or an administrator refreshes out-of-date demo studies; anyone else is told who can. See the package README, Demo data.
//output: dataframe report
export async function ncaSeedDemoStudies(): Promise<DG.DataFrame> {
  return seedDemoStudies(productionSeedDeps((file) => loadSampleByName(`${DEMO_STUDIES_DIR}/${file}`)));
}
/* eslint-enable max-len */
