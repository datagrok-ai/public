export namespace decorators {
  /** A function that registers [DG.JsViewer](https://datagrok.ai/js-api/dg/classes/JsViewer).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param icon - path to an icon file from the package root.
   * @param toolbox - set to true to add the viewer icon to the toolbox.
   * @param trellisable - a flag that indicates whether this viewer can be an inner viewer of [DG.VIEWER.TRELLIS_PLOT]
   * @param viewerPath - viewer path in the top menu, should include the viewer name (Add | JavaScript Viewers | \<ViewerPath\>).
   *   The default path is "Add > JavaScript Viewers > \<Package Name\> > \<Friendly Viewer Name\>".
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/develop-custom-viewer}
   * 
   * Usage examples:
   * 
   * ```ts
   * @viewer()
   * class TestViewer {
   *   constructor() {
   *     console.log('Viewer constructed');
   *   }
   * }
   * 
   * @viewer({
   *   name: 'Test Viewer',
   *   description: 'Creates a Test Viewer instance',
   *   icon: 'images/icon.png',
   *   toolbox: true,
   * })
   * class TestViewer {
   *   constructor() {
   *     console.log('Viewer constructed');
   *   }
   * }
   * ```
   */
  export function viewer(options?: {
    name?: string,
    description?: string,
    icon?: string,
    toolbox?: boolean,
    trellisable?: boolean,
    viewerPath?: string
  }) {
    return function (constructor: Function) { };
  }

  /** A function that registers [DG.Filter](http://datagrok.ai/js-api/dg/classes/Filter).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param semType - semantic type name (e.g., "Molecule"). See [DG.SEMTYPE](https://datagrok.ai/js-api/enums/dg.SEMTYPE)
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/custom-filters}
   */
  export function filter(options?: {
    name?: string,
    description?: string,
    semType?: string,
  }) {
    return function (constructor: Function) { };
  }

  /** A function that registers [DG.GridCellRender](https://datagrok.ai/js-api/grok/namespaces/decorators/functions/cellRenderer).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param cellType - cell type name (e.g., "html", "image").
   * @param columnTags - a string of column tags required for a match.
   *   Key-value pairs look like this: "quality=Macromolecule, units=separator".
   * @param virtual - a flag to enable rendering in virtual columns.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/custom-cell-renderers}
   */
  export function cellRenderer(options?: {
    name?: string,
    description?: string,
    cellType: string,
    columnTags?: string,
    virtual?: boolean,
  }) {
    return function (constructor: Function) { };
  }

  interface InputOptions {
    semType?: string;
    category?: string;
    optional?: boolean;
    editor?: string;
    nullable?: boolean;
    separators?: string[];
    choices?: string[] | string;
    format?: string;
    min?: string;
    max?: string;
    caption?: string;
    description?: string;
    initialValue?: string;
    viewer?: string;
    units?: string;
    type?: string;
    optionsType?: string;
    step?: string;
    "meta.url"?: boolean;
    metaUrl?: boolean;
  }

  interface Input extends InputOptions{
    name?: string;
    type?: string;
    options?: InputOptions
  }

  interface Output {
    name: string;
    type: string;
    options?: Record<string, string>;
  }

  interface Meta {
    cache?: string;
    cacheInvalidateOn?: string;
    ['cache.invalidateOn']?: string;
    browsePath?: string;
    icon?: string;
    demoPath?: string;
    demoSkip?: string;
    demoWait?: string;
    path?: string;
    vectorFunc?: string;
    ext?: string;
    cellType?: string;
    columnTags?: string;
    supportedSemTypes?: string;
    supportedTypes?: string;
    supportedDistanceFunctions?: string;
    supportedUnits?: string;
    action?: string;
    fileViewerCheck?: string;
    fileViewer?: string;
    keywords?: string;
    role?: string;
    mlname?: string;
    mlrole?: string;
    inputRegex?: string;
    runOnOpen?: string;
    runOnInput?: string;
    features?: string;
    toolbox?: string;
    gridChart?: string;
    virtual?: string;
    order?: string;
    autostartImmediate?: string;
    ['scriptHandler.language']?: string;
    scriptHandlerLanguage?: string;
    ['scriptHandler.extensions']?: string;
    scriptHandlerExtensions?: string;
    ['scriptHandler.commentStart']?: string;
    scriptHandlerCommentStart?: string;
    ['scriptHandler.templateScript']?: string;
    scriptHandlerTemplateScript?: string;
    ['scriptHandler.codeEditorMode']?: string;
    scriptHandlerCodeEditorMode?: string;
    ['scriptHandler.vectorizationFunction']?: string;
    scriptHandlerVectorizationFunction?: string;
    url?: string;
    propertyType?: string;
    semType?: string;
    /** Runs the function server-side as a celery task in a docker-based Node worker.
     * `true` — an auto-created worker container; a string — the package's `dockerfiles/<name>` container. */
    queue?: boolean | string;
    /** Runs the function server-side in the default worker queue — an alias of `queue: true`. */
    server?: boolean;
  }

  interface FunctionOptions {
    name?: string,
    friendlyName?: string,
    tags?: string[],
    description?: string,
    meta?: Meta | Record<string, string>,
    outputs?: Output[],
    result?: Output,
    sidebar?: string;
    editor?: string;
    cache?: string;
    ['cache.invalidateOn']?: string;
    cacheInvalidateOn?: string;
    ['top-menu']?: string;
    topMenu?: string;
    condition?: string;
    helpUrl?: string;
    ['help-url']?: string;
    connection?: string;
  }

  interface AppOptions extends FunctionOptions {
    browsePath?: string,
    icon?: string,
    url?: string
  }

  interface DashboardOptions extends FunctionOptions {
    test?: string,
    order?: string,
  }

  interface ModelOptions extends FunctionOptions {
    icon?: string,
    features?: string,
    runOnInput?: string,
    runOnOpen?: string
  }

  interface AppTreeBrowserOptions extends FunctionOptions {
    app?: string;
  }

  interface CellRendererOptions extends FunctionOptions {
    cellType?: string,
    columnTags?: string
  }

  interface DashboardOptions extends FunctionOptions {
    order?: string
  }

  interface FileViewerOptions extends FunctionOptions {
    fileViewer: string;
    fileViewerCheck?: string;
  }

  interface FileHandlerOptions extends FunctionOptions {
    ext: string;
    fileViewerCheck?: string;
  }

  interface DemoOptions extends FunctionOptions {
    path?: string;
    demoPath?: string;
    demoSkip?: string;
    demoWait?: string;
    test?: { test: string, wait: string, timeout?: string, skip?: string }
  }

  export function func(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function app(config: AppOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function autostart(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function init(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function editor(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function panel(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function dashboard(config?: DashboardOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function folderViewer(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function semTypeDetector(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function packageSettingsEditor(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  // export function cellRenderer(config: CellRendererOptions) {
  //   return function (
  //     target: any,
  //     propertyKey: string,
  //     descriptor: PropertyDescriptor
  //   ) { };
  // }

  export function functionAnalysis(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function converter(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function fileViewer(config?: FileViewerOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function fileExporter(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function fileHandler(config: FileHandlerOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function demo(config: DemoOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function treeBrowser(config?: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function model(config: ModelOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function appTreeBrowser(config?: AppTreeBrowserOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function param(options: Input) {
    return function (
      target: Object,
      propertyKey: string | symbol,
      parameterIndex: number
    ): void { };
  }

  /** Options of {@link entity}: the table a class maps to in the package's domain schema
   * (`databases/<schema>/schema.json`), with every table key the manifest knows. */
  export interface EntityOptions {
    /** The domain schema — the `databases/<schema>/` directory the manifest lands in. */
    schema: string;
    /** Table name; defaults to the class name in snake_case (`PlateWell` → `plate_well`). */
    table?: string;
    securityMode?: 'table' | 'master' | 'row';
    promotion?: 'lazy' | 'eager';
    defaultRowVisibility?: 'table' | 'none';
    delegate?: string;
    softDelete?: boolean;
    audit?: boolean;
    extensible?: boolean;
    idempotency?: boolean;
    ginIndex?: boolean;
    businessKey?: string[];
    friendlyName?: string;
    description?: string;
    singularName?: string;
    pluralName?: string;
    filters?: {column: string, type?: 'categories' | 'histogram' | 'range' | 'text' | 'bool', bins?: number, label?: string}[];
    relations?: {[name: string]: {via: string, target: string, viaSelf?: string, viaTarget?: string, allowCreate?: boolean}};
    schemas?: string[];
    /** Table CHECK constraints: name → `{check}` over this table's own relational columns
     * (`effective_to IS NULL OR effective_from < effective_to`). */
    constraints?: {[name: string]: {check: string}};
    /** Declared grants on the table's entity: group → permissions. */
    grants?: {[group: string]: ('view' | 'edit' | 'delete')[]};
  }

  export type EntityColumnType = 'string' | 'int' | 'float' | 'bool' | 'datetime' | 'string_list' | 'ref' | 'user' |
    'group' | 'file' | 'json';

  /** Options of {@link column}: one manifest column. */
  export interface ColumnOptions {
    /** Manifest column name; defaults to the property name. */
    name?: string;
    /** Column type; inferred from the property type when omitted (`string`, `number` → `float`,
     * `boolean` → `bool`, `Date` → `datetime`, `string[]` → `string_list`, `object` → `json`). */
    type?: EntityColumnType;
    /** Target table of a reference: its name, or `() => EntityClass` for an `@entity` class of
     * the same schema. Makes the column a `ref` unless the type is `user` or `group`. */
    ref?: string | (() => Function);
    onDelete?: 'cascade' | 'restrict' | 'setnull';
    required?: boolean;
    unique?: boolean;
    isName?: boolean;
    /** Write-once: the first non-null value wins; a later write may only repeat it. */
    immutable?: boolean;
    min?: number;
    max?: number;
    choices?: string[];
    semType?: string;
    friendlyName?: string;
    description?: string;
    default?: any;
    format?: string;
  }

  /** Maps a class to a table of the package's domain schema. A runtime no-op: `grok schema
   * generate` reads the decorated classes statically and writes
   * `databases/<schema>/schema.json` from them (property order = column order); `grok schema
   * seal|diff|migrate` then derive the snapshot and the migration scripts as usual.
   *
   * ```ts
   * @grok.decorators.entity({schema: 'lab', securityMode: 'row', businessKey: ['barcode']})
   * export class Plate {
   *   @grok.decorators.column({required: true, unique: true, isName: true})
   *   barcode!: string;
   *
   *   @grok.decorators.column({type: 'int', min: 1})
   *   rows!: number;
   *
   *   @grok.decorators.column({ref: () => Reader, onDelete: 'setnull'})
   *   reader_id?: string;
   * }
   * ```
   */
  export function entity(options: EntityOptions) {
    return function (constructor: Function) { };
  }

  /** Maps a property to a column of its {@link entity} table; a runtime no-op read by
   * `grok schema generate`. Undecorated properties are not columns. */
  export function column(options?: ColumnOptions) {
    return function (target: any, propertyKey: string | symbol): void { };
  }
}
