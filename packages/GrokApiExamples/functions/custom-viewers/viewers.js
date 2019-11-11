// This viewer does the following:
// * defines two properties, "question" and "answer". Properties are persistable and editable.
// * listens to changes in properties, attached table's selection and filter, and updates accordingly.
class JsDemoViewer extends JsViewer {

    constructor() {
        super();

        // Register properties and define fields initialized to properties' default values
        this.question = this.string('question', 'life');
        this.answer = this.int('answer', 42);

        // Properties that represent columns should be strings with the 'ColumnName' postfix
        this.valueColumnName = this.string('valueColumnName', null);
        this.registerCleanup(uit.handleResize(this.root, (w, h) => this.render()));
    }

    // override to handle property changes
    onPropertyChanged(prop) {
        gr.balloon.info(`${prop.name}: ${prop.get(this)}`);
    }

    onFrameAttached(dataFrameHandle) {
        this.dataFrame = new DataFrame(dataFrameHandle);
        this.dataFrame.selection.onChanged(() => this.render());
        this.dataFrame.filter.onChanged(() => this.render());

        this.render();
    }

    render() {

        this.root.innerHTML =
            `${this.dataFrame.toString()}<br>
            Column: ${this.valueColumnName}<br>
            Question: ${this.question}<br>
            Answer: ${this.answer}<br>
            Selected: ${this.dataFrame.selection.trueCount}<br>
            Filtered: ${this.dataFrame.filter.trueCount}<br>
            Size: ${this.root.clientWidth} x ${this.root.clientHeight}`;
    }
}

// Register viewer with the platform. This enables the following:
// * Add viewer from Add | JsDemoViewer, or from the toolbar 'viewers' popup
// * Persist viewer as part of the layout
// * Common viewer operations under the "Viewer" popup menu, such as cloning, embedding, etc
gr.registerViewer('JsDemoViewer', 'JavaScript-based viewer', () => new JsDemoViewer());

// Add viewer to a table view
demog = gr.testData('demog', 5000);
view = gr.addTableView(demog);
hist = view.addViewer('JsDemoViewer');