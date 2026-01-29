export interface FuncRoleDescription {
  role: string;
  description: string;
  header: string;
  signature?: string;
}

export const FUNC_TYPES = {
  /** An application that gets shown in the app store.
    * Signature: app() */
  APP: 'app',

  /** Context-specific widget that appears on the context panel
    * Signature: panel(x: any): Widget */
  PANEL: 'panel',

  /** Gets invoked when the containing package is initialized
    * Signature: init() */
  INIT: 'init',

  /** Gets invoked at platform startup. Use it wisely as the whole package will get initialized.
    * Signature: autostart() */
  AUTOSTART: 'autostart',

  /** Semantic type detector for a column. Gets invoked when a new dataframe is imported into the platform.
   *  Implementation should either set column.semType directly, or return the semantic type that will get assigned.
   *  Signature: semTypeDetector(Column): string */
  SEM_TYPE_DETECTOR: 'semTypeDetector',

  /** Creates a viewer (or editor) for a file with the specified extension.
   *  The extension is derived from the `fileViewer-<extension>` tag.
   *  Used in the file system browser.
   *  Signature: fileViewer(FileInfo): View */
  FILE_VIEWER: 'fileViewer',

  /** Exports a file. Gets added to the "export" menu at startup.
   *  Signature: fileExporter() */
  FILE_EXPORTER: 'fileExporter',

  /** Handles custom file formats.
   * The `meta.ext` parameter should contain a comma-separated list of extensions.
   * Signature: fileImporter(x: string | TypedArray): DataFrame[] */
  FILE_IMPORTER: 'file-handler',

  /** Creates a cell renderer that is used for rendering cells for specific semantic types.
   *  Semantic type is derived from the `cellRenderer-<semType>` tag.
   *  Signature: cellRenderer(): GridCellRenderer */
  CELL_RENDERER: 'cellRenderer',

  /** Edits package settings.
   *  Signature: packageSettingsEditor(): Widget */
  PACKAGE_SETTINGS_EDITOR: 'packageSettingsEditor',

  /** Makes a widget appear on the welcome screen
   *  Signature: dashboard(): DG.Widget */
  DASHBOARD: 'dashboard',

  /**
   * Function analysis. Examples: sensitivity analysis, parameter editor
   * Func => View */
  FUNCTION_ANALYSIS: 'functionAnalysis',

  /** Converts values. Has one input and one output */
  CONVERTER: 'converter',

  WIDGET: 'widget',
  WIDGETS: 'widgets',
  EDITOR: 'editor',
  TRANSFORM: 'Transform',
  FILTER: 'filter',
  VIEWER: 'viewer',
  VALUE_EDITOR: 'valueEditor',
  CELL_EDITOR: 'cellEditor',
  UNIT_CONVERTER: 'unitConverter',
  MOLECULE_SKETCHER: 'moleculeSketcher',
  TOOLTIP: 'tooltip',
  FOLDER_VIEWER: 'folderViewer',
  SCRIPT_HANDLER: 'scriptHandler',

  HIT_TRIAGE_FUNCTION: 'HitTriageFunction',
  HIT_TRIAGE_DATA_SOURCE: 'HitTriageDataSource',
  HIT_TRIAGE_SUBMIT_FUNCTION: 'HitTriageSubmitFunction',
  HIT_DESIGNER_FUNCTION: 'HitDesignerFunction',

  DIM_RED_PREPROCESS: 'dim-red-preprocessing-function',
  DIM_RED_POSTPROCESS: 'dim-red-postprocessing-function',

  MONOMER_LIB_PROVIDER: 'monomer-lib-provider',

  SEARCH_PROVIDER: 'searchProvider',
  NOTATION_REFINER: 'notationRefiner',
};

export const functionRoles: FuncRoleDescription[] = [
  {
    role: FUNC_TYPES.APP,
    description: 'An application that gets shown in the app store.',
    header: 'tags',
    signature: 'app(): void | View',
  },
  {
    role: FUNC_TYPES.PANEL,
    description: 'Context-specific widget that appears on the context panel.',
    header: 'tags',
    signature: 'panel(...args): Widget | Viewer | graphics | void',
  },
  {
    role: FUNC_TYPES.INIT,
    description: 'Gets invoked when the containing package is initialized.',
    header: 'tags',
    signature: 'init(): void',
  },
  {
    role: FUNC_TYPES.AUTOSTART,
    description: 'Gets invoked at platform startup. Use it wisely as the whole package will get initialized.',
    header: 'tags',
    signature: 'autostart(): void',
  },
  {
    role: FUNC_TYPES.SEM_TYPE_DETECTOR,
    description: 'Semantic type detector for a column. Gets invoked when a new dataframe is imported into the platform.\n   *  Implementation should either set column.semType directly, or return the semantic type that will get assigned.',
    header: 'tags',
    signature: 'semTypeDetector(col: Column): string',
  },
  {
    role: FUNC_TYPES.FILE_VIEWER,
    header: 'tags',
    description: 'Creates a viewer (or editor) for a file with the specified extension.\n   *  The extension is derived from the `fileViewer-[extension]` tag.\n   *  Used in the file system browser.',
    signature: 'fileViewer(file: FileInfo): View',
  },
  {
    role: FUNC_TYPES.FILE_EXPORTER,
    header: 'tags',
    description: 'Exports a file. Gets added to the "export" menu at startup.',
    signature: 'fileExporter(): void',
  },
  {
    role: FUNC_TYPES.FILE_IMPORTER,
    header: 'tags',
    description: 'Handles custom file formats.\n   * The `meta.ext` parameter should contain a comma-separated list of extensions',
    signature: 'fileImporter(x: string | TypedArray): DataFrame[]',
  },
  {
    role: FUNC_TYPES.CELL_RENDERER,
    header: 'tags',
    description: 'Creates a cell renderer that is used for rendering cells for specific semantic types.\n   *  Semantic type is derived from the `cellRenderer-[semType]` tag.',
    signature: 'cellRenderer(): GridCellRenderer',
  },
  {
    role: FUNC_TYPES.PACKAGE_SETTINGS_EDITOR,
    header: 'tags',
    description: 'Edits package settings.',
    signature: 'packageSettingsEditor(): Widget',
  },
  {
    role: FUNC_TYPES.DASHBOARD,
    description: 'Makes a widget appear on the welcome screen.',
    header: 'tags',
    signature: 'dashboard(): Widget',
  },
  {
    role: FUNC_TYPES.FUNCTION_ANALYSIS,
    description: 'Function analysis that gets added to the function view. Examples: sensitivity analysis, parameter editor',
    header: 'tags',
    signature: 'functionAnalysis(x: Function): View',
  },
  {
    role: FUNC_TYPES.CONVERTER,
    description: 'Converts values. Has one input and one output',
    header: 'role',
    signature: 'converter(x: any): any',
  },
  {
    role: FUNC_TYPES.WIDGET,
    description: 'Creates a custom widget.',
    header: 'tags',
    signature: 'widget(...args): Widget',
  },
  {
    role: FUNC_TYPES.WIDGETS,
    description: 'Creates a custom widget.',
    header: 'tags',
    signature: 'widgets(...args): Widget',
  },
  {
    role: FUNC_TYPES.EDITOR,
    description: 'Creates a custom editor for a function call',
    header: 'tags',
    signature: 'editor(call: FuncCall): Widget | View | void',
  },
  {
    role: FUNC_TYPES.TRANSFORM,
    description: '',
    header: 'tags',
    signature: 'transform(table: DataFrame, ...args): any',
  },
  {
    role: FUNC_TYPES.FILTER,
    description: 'Creates a custom table filter',
    header: 'tags',
    signature: 'filter(): Filter',
  },
  {
    role: FUNC_TYPES.VIEWER,
    description: 'Creates a custom viewer',
    header: 'tags',
    signature: 'viewer(): Viewer',
  },
  {
    role: FUNC_TYPES.VALUE_EDITOR,
    description: 'Custom editor for specific input types',
    header: 'tags',
    signature: 'valueEditor(...args): any',
  },
  {
    role: FUNC_TYPES.CELL_EDITOR,
    description: 'CCustom editor for grid cells',
    header: 'tags',
    signature: 'cellEditor(cell: GridCell): void',
  },
  {
    role: FUNC_TYPES.UNIT_CONVERTER,
    description: 'Converts between measurement units',
    header: 'tags',
    signature: 'unitConverter(value: string, source: string, target: string): string',
  },
  {
    role: FUNC_TYPES.MOLECULE_SKETCHER,
    description: 'Creates a molecule sketcher widget',
    header: 'tags',
    signature: 'moleculeSketcher(): Widget',
  },
  {
    role: FUNC_TYPES.TOOLTIP,
    description: 'Provides a custom tooltip for columns',
    header: 'tags',
    signature: 'tooltip(col: Column): Widget',
  },
  {
    role: FUNC_TYPES.FOLDER_VIEWER,
    description: 'Provides a custom folder content preview.',
    header: 'tags',
    signature: 'folderViewer(folder: File, files: list<file>): Widget',
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_FUNCTION,
    description: 'Compute function for Hit Triage campaigns.',
    header: 'tags',
    signature: 'hitTriageFunction(table: DataFrame, moleculeCol: Column): DataFrame',
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_DATA_SOURCE,
    description: 'Provides a datasource for Hit Triage campaigns. Must return a dataframe containing molecules.',
    header: 'tags',
    signature: 'hitTriageDataSource(...args): DataFrame',
  },
  {
    role: FUNC_TYPES.HIT_TRIAGE_SUBMIT_FUNCTION,
    description: 'Processes or saves the computed dataset.',
    header: 'tags',
    signature: 'hitTriageSubmitFunction(df: DG.DataFrame, moleculesCol: string): void',
  },
  {
    role: FUNC_TYPES.HIT_DESIGNER_FUNCTION,
    description: 'Compute function for Hit Design campaigns.',
    header: 'tags',
    signature: 'hitDesignerFunction(molecule: string, ...args): DataFrame',
  },
  {
    role: FUNC_TYPES.DIM_RED_PREPROCESS,
    description: 'Preprocessing function for dimensionality reduction.',
    header: 'tags',
    signature: 'preprocess(col: Column, metric: string, ...args): any',
  },
  {
    role: FUNC_TYPES.DIM_RED_POSTPROCESS,
    description: 'Postprocessing function for dimensionality reduction.',
    header: 'tags',
    signature: 'postprocess(xCol: Column, yCol: Column, ...args): void',
  },
  {
    role: FUNC_TYPES.MONOMER_LIB_PROVIDER,
    description: 'Provides a monomer library provider.',
    header: 'tags',
    signature: 'getMonomerLibProvider(): IMonomerLibProvider',
  },
  {
    role: FUNC_TYPES.SEARCH_PROVIDER,
    description: 'Marks a function to be used as a search provider in the global search.',
    header: 'tags',
    signature: 'searchProvider(): SearchProvider',
  },
  {
    role: FUNC_TYPES.NOTATION_REFINER,
    description: 'Refines the biological sequence notation based on company specific rules',
    header: 'tags',
    signature: 'notationRefiner(column: Column, stats: any, separator: string): bool',
  },
];
