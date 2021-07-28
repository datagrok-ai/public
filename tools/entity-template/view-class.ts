import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// A sample class from the Notebooks package:
// https://github.com/datagrok-ai/public/tree/master/packages/Notebooks
// This class defines a new view for Jupyter Notebooks.
export class #{NAME} extends DG.ViewBase {
  TYPE: string;
  PATH: string;
  notebookId: any;

  constructor(params: object | null, path: string) {
    super(params, path);
    this.TYPE = 'Notebook';
    this.PATH = '/notebook';
  }

  // Override basic methods
  get type() {
    return this.TYPE;
  }

  get helpUrl() {
    return '/help/develop/jupyter-notebook.md';
  }

  get name() {
    return 'Notebook';
  }

  get path() {
    return `${this.PATH}/${this.notebookId}`;
  }

  // Icon
  getIcon() {
    let img = document.createElement('img');
    img.src = '/images/entities/jupyter.png';
    img.height = 18;
    img.width = 18;
    return img;
  }

  // View state serialization/deserialization
  saveStateMap() {
    return {'notebookId': this.notebookId};
  }

  loadStateMap(stateMap: { [x: string]: string }) {
    open(stateMap['notebookId']);
  }

  // URL path handler
  handlePath(path: string) {
    let id = path.replace(`${this.PATH}/`, '');
    open(id);
  }

  // URL path checker
  acceptsPath(path: string) {
    return path.startsWith(this.PATH);
  }
}
