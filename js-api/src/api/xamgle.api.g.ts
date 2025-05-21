/// this file was generated automatically from xamgle classes declarations
import { toDart } from "../wrappers";
let api = (typeof window !== 'undefined' ? window : global.window) as any;

export interface SettingsInterface {
  /// Automatically save workspace locally
  autoSaveWorkspace: boolean;

  /// Automatically detect column semantic types
  autoDetectSemanticTypes: boolean;

  /// Show '?' link on docked viewers
  dockShowHelpLink: boolean;

  /// Show '⚙' settings icon on docked viewers
  dockShowSettingsLink: boolean;

  /// Show hamburger link on docked viewers
  dockShowMenuLink: boolean;

  /// Show notification when a file is imported
  notifyOnFileImport: boolean;

  allowAsyncFileImport: boolean;

  /// Controls whether tables pane automatically shows when more
  /// than one table is open
  autoShowTablesPane: boolean;

  /// Controls whether columns pane automatically shows when a table is open
  autoShowColumnsPane: boolean;

  /// Show viewer settings in the property dialog when user clicks on a viewer.
  showViewerSettingsOnClick: boolean;

  /// Show viewer settings in the property dialog when a viewer is added.
  showViewerSettingsOnAddition: boolean;

  /// Show user icon on top.
  showUserIcon: boolean;

  /// Always show filters icons, not only on hover
  showFiltersIconsConstantly: boolean;

  /// Determines when to show range sliders
  showRangeSlidersOnViewers: string;

  /// Auto-apply existing layout after selected rows are extracted
  applyLayoutWhenExtractingRows: boolean;

  /// Persist history of actions along with tables and columns
  dataHistory: boolean;

  isServer: boolean;

  showRecentlyOpenedViewsInHistory: boolean;

  warnOnUnsavedChanges: boolean;

  hiddenMenus: Array<string>;

  showCurrentRowInProperties: boolean;

  showFilteredRowsInProperties: boolean;

  showSelectedRowsInProperties: boolean;

  showSelectedColumnsInProperties: boolean;

  showCurrentColumnInProperties: boolean;

  showMenu: boolean;

  showTables: boolean;

  showColumns: boolean;

  showProperties: boolean;

  showToolbox: boolean;

  showBrowse: boolean;

  showProjects: boolean;

  showFavorites: boolean;

  showStatusBar: boolean;

  showVariables: boolean;

  showConsole: boolean;

  showHelp: boolean;

  enableBetaViewers: boolean;

  saveProjectWithViewLayout: boolean;

  allowWidgetsAsColumns: boolean;

  allowEventScripts: boolean;

  /// Displays Index Files check box in connection edit dialog. If it is disabled then check
  /// box will not be shown to the user and indexing of files will not be configurable.
  enableConnectionIndexFiles: boolean;

  // When true, use WebGPU for rendering and advanced computations if possible.
  enableScatterPlotWebGPUAcceleration: boolean;

  dateFormat: string;

  integerNumberFormat: string;

  floatingPointNumberFormat: string;

  hiddenPanels: Array<string>;

  panelOrder: string;

  /// Will load default settings from server on platform start
  loadDefaultsOnStart: boolean;

  /// If not empty, will be used to load package.js instead of the Datagrok backend
  webpackDevUrl: string;

  clientFuncCacheEnabled: boolean;

  clientFilesCacheEnabled: boolean;

  dataFrameBatchSize: number;

  /// Whenever an error occurs, automatically create report with extended logs and send it to the server. Use UsageAnalysis to browse them.
  /// This includes console logs, server logs, data connectivity logs, Docker logs, etc.
  autoReportErrors: boolean;


}
