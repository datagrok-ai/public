export class DockConfig {
    public constructor() {
        this.escClosesWindow = true;
        this.escClosesDialog = true;
        this.dialogRootElement = document.body;
        this.moveOnlyWithinDockConatiner = false;
        this.enableBrowserWindows = true;
    }

    escClosesWindow?: boolean;
    escClosesDialog?: boolean;
    dialogRootElement: HTMLElement;
    moveOnlyWithinDockConatiner?: boolean;
    enableBrowserWindows?: boolean;
    // When false, docking a panel into an existing tab strip does not switch focus
    // to it (the first panel of an empty strip is still activated). Default: true.
    activatePanelOnAdd?: boolean;
}
