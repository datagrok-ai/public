/**
 * Views barrel — the public surface of the view layer.
 */

export {NcaAppView} from './nca-app-view';

export {
  loadSampleByName, resolveAnalysisName, launchAnalysis, mountAnalysisIntoView,
  ensureDemoTables,
} from './nca-launch';

// Bundled demo-table catalogue (leaf module — imports nothing; see its header).
export {DEMO_TABLE_FILES, DEMO_ID_TAG, PROVENANCE_TAG} from './nca-demo-tables';
export type {DemoTableEntry} from './nca-demo-tables';

// Hub onboarding help panel (nca-hub-onboarding).
export {
  openHubHelp, restoreHubHelpSync, buildHubHelpElement, HUB_HELP_MD,
  INFO_ICON_TOKEN, ADD_ICON_TOKEN, SHARE_ICON_TOKEN,
  DEMO_BALLOON_MD, buildDemoBalloonElement,
} from './nca-hub-help';

export {mountMappingPane} from './nca-mapping-pane';
export type {MountMappingPaneOptions, MappingPaneHandle, ComputeCallbackCtx} from './nca-mapping-pane';

export {configureResultsDataFrame, configureResultsGrid, openResultsView} from './results-view';

// NCA Studies workspace (nca-projects-workspace).
export {
  parseStudyPath, updateNcaStudyPath, routeNcaApp,
} from './nca-studies-router';
export type {RouteState} from './nca-studies-router';
export {messageView, openAnalysis} from './nca-analysis-open';
export {renderStudiesListInto, openStudyOverview, createStudyFlow} from './nca-studies-list';
export {renderStudyOverview} from './nca-study-overview';
export {openShareDialog} from './nca-share-dialog';
export {openAddToStudyDialog} from './nca-add-to-study';
export type {AddToStudyInput} from './nca-add-to-study';
