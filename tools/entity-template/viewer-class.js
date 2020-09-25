// A sample class from https://datagrok.ai/help/develop/js-api#custom-viewers
// This viewer does the following:
// * listens to changes of filter and selection in the attached table,
// * updates the number of filtered/selected rows accordingly.
export class #{NAME_TITLECASE}Viewer extends DG.JsViewer {
    onFrameAttached() {
        subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
        subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));

        this.render();
    }

    render() {
        this.root.innerHTML =
            `${this.dataFrame.toString()}<br>
            Selected: ${this.dataFrame.selection.trueCount}<br>
            Filtered: ${this.dataFrame.filter.trueCount}`;
    }
}
