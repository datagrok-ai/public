export namespace decorators {

  /** A function that registers [DG.JsViewer].
   * @param name - function name in UI.
   * @param description - function description in UI.
   * @param icon - path to an icon file from the package root.
   * @param toolbox - set to true to add the viewer icon to the toolbox.
   * 
   * See also: {@link https://datagrok.ai/help/develop/how-to/develop-custom-viewer}
   * 
   * Usage examples:
   * 
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
   */
  export function grokViewer(options?: {
    name?: string,
    description?: string,
    icon?: string,
    toolbox?: boolean,
  }) {
    return function(constructor: Function) {};
  }

  export function grokCellRenderer(options?: {
    name?: string,
    description?: string,
  }) {}
}
