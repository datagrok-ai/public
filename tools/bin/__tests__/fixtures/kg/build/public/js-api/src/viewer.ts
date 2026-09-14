/** A visualization of a table. */
export class Viewer {
  get root(): HTMLElement {
    return document.createElement('div');
  }
}

export class JsViewer extends Viewer {
  onTableAttached(): void {
  }
}
