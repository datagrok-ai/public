// A sample class from the Notebooks package:
// https://github.com/datagrok-ai/public/tree/master/packages/Notebooks
// This class defines a new view for Jupyter Notebooks.
export class #{NAME} extends DG.ViewBase {
    constructor(params, path) {
        super(params, path);
        this.TYPE = 'Notebook';
        this.PATH = '/notebook';
    }
    
    // Override basic metods
    get type() { return this.TYPE };
    get helpUrl() { return '/help/compute/jupyter-notebook.md'; }
    get name() { return 'Notebook' };
    get path() { return `${this.PATH}/${this.notebookId}` };
    
    // Icon
    getIcon() {
        let img = document.createElement('img');
        img.src = '/images/entities/jupyter.png';
        img.height = 18;
        img.width = 18;
        return img;
    };

    // View state serialization/deserialization
    saveStateMap() { return {'notebookId': this.notebookId }; }
    loadStateMap(stateMap) { open(stateMap['notebookId']); }
    
    // URL path handler
    handlePath(path) {
        let id = path.replace(`${this.PATH}/`, '');
        open(id);
    }
    
    // URL path checker
    acceptsPath(path) { return path.startsWith(this.PATH); }
}
