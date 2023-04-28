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
}
