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
    choices?: string[];
    format?: string;
    min?: string;
    max?: string;
    caption?: string;
    description?: string;
    initialValue?: string;
    viewer?: string;
    units?: string;
  }

  interface Input {
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
    ['scriptHandler.extensions']?: string;
    ['scriptHandler.commentStart']?: string;
    ['scriptHandler.templateScript']?: string;
    ['scriptHandler.codeEditorMode']?: string;
    ['scriptHandler.vectorizationFunction']?: string;
    url?: string;
    propertyType?: string;
    semType?: string;
  }

  interface FunctionOptions {
    name?: string,
    tags?: string[],
    description?: string,
    meta?: Meta | Record<string, string>,
    outputs?: Output[],
    sidebar?: string;
    editor?: string;
    cache?: string;
    ['cache.invalidateOn']?: string;
  }

  interface AppOptions extends FunctionOptions{
    browsePath?: string,
    icon?: string, 
    url?: string
  }

  interface ModelOptions extends FunctionOptions{
    icon?: string,
    features?: string,
    runOnInput?: string,
    runOnOpen?: string
  }

  interface CellRendererOptions extends FunctionOptions{
    cellType?: string,
    columnTags?: string
  }

  interface DashboardOptions extends FunctionOptions{
    order?: string
  }

  interface FileViewerOptions extends FunctionOptions{
    fileViewer: string;
    fileViewerCheck?: string;
  }
  
  interface FileHandlerOptions extends FunctionOptions{
    ext: string;
    fileViewerCheck?: string;
  }
  
  interface DemoOptions extends FunctionOptions{
    path?: string;
    demoPath?: string;
    demoSkip?: string;
    demoWait?: string;
    test?: { test: string, wait: string, timeout?: string, skip?: string }
  }

  export function func(config: FunctionOptions) {
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

  export function autostart(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function init(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function editor(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function panel(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function folderViewer(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function semTypeDetector(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function packageSettingsEditor(config: FunctionOptions) {
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

  export function dashboard(config: DashboardOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function functionAnalysis(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function converter(config: FunctionOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function fileViewer(config: FileViewerOptions) {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  }

  export function fileExporter(config: FunctionOptions) {
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

  export function treeBrowser(config: FunctionOptions) {
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

  export function param(options: Input) {
    return function (
      target: Object,
      propertyKey: string | symbol,
      parameterIndex: number
    ): void { };
  }
}
