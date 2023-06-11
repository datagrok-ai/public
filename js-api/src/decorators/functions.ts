export namespace decorators {

  /** A function that registers [DG.JsViewer](https://datagrok.ai/js-api/classes/dg.JsViewer).
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
    viewerPath?: string,
  }) {
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.Filter](https://datagrok.ai/js-api/classes/dg.Filter).
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
    return function(constructor: Function) {};
  }

  /** A function that registers [DG.GridCellRender](https://datagrok.ai/js-api/classes/dg.GridCellRenderer).
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
    return function(constructor: Function) {};
  }
}
