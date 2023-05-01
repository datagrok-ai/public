export namespace decorators {

  /** A function that registers [DG.JsViewer](https://datagrok.ai/js-api/classes/dg.JsViewer).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param icon - path to an icon file from the package root.
   * @param toolbox - set to true to add the viewer icon to the toolbox.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/develop-custom-viewer}
   * 
   * Usage examples:
   * 
   * ```ts
   * @grokViewer()
   * class TestViewer {
   *   constructor() {
   *     console.log('Viewer constructed');
   *   }
   * }
   * 
   * @grokViewer({
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
  export function grokViewer(options?: {
    name?: string,
    description?: string,
    icon?: string,
    toolbox?: boolean,
  }) {
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.GridCellRender](https://datagrok.ai/js-api/classes/dg.GridCellRenderer).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param cellType - cell type name (e.g., "html", "image").
   * @param columnTags - a string of column tags required for a match.
   *   Key-value pairs look like this: "quality=Macromolecule, units=separator".
   * @param virtual - flag ??
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/custom-cell-renderers}
   */
  export function grokCellRenderer(options?: {
    name?: string,
    description?: string,
    cellType: string,
    columnTags?: string,
    virtual?: boolean,
  }) {
    return function(constructor: Function) {};
  }
  
  /** A function that registers [DG.FUNC_TYPES.FILE_EXPORTER](https://datagrok.ai/js-api/modules/dg#func_types).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param extension - exported file extension.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/file-exporters}
   */
  export function grokFileExporter(options?: {
    name?: string,
    description?: string,
    extension: string,
  }) {
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.FUNC_TYPES.FILE_IMPORTER](https://datagrok.ai/js-api/modules/dg#func_types).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param extensions - supported file format or a list of formats.
   * @param inputType - file content type: "string" | "list". Default is "string".
   * Specify only if a list of bytes is expected as an input.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/file-handlers}
   */
  export function grokFileImporter(options?: {
    name?: string,
    description?: string,
    extensions: string | string[],
    inputType?: string,
  }) {
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.FUNC_TYPES.FILE_VIEWER](https://datagrok.ai/js-api/modules/dg#func_types).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param extensions - supported file format or a list of formats.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/create-custom-file-viewers}
   */
  export function grokFileViewer(options?: {
    name?: string,
    description?: string,
    extensions: string | string[],
  }) {
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.FUNC_TYPES.PACKAGE_SETTINGS_EDITOR](https://datagrok.ai/js-api/modules/dg#func_types).
   * @param name - function name in UI.
   * @param description - function description in UI.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/custom-package-settings-editors}
   */
    export function grokSettingsEditor(options?: {
      name?: string,
      description?: string,
    }) {
      return function(constructor: Function) {};
    }
}
