import { FillDockContainer } from "./FillDockContainer";
import { TabHost } from "./TabHost";
import { DocumentTabPage } from "./DocumentTabPage";
import { DockManager } from "./DockManager";
import { TabHostDirection } from "./enums/TabHostDirection";
import { IState } from "./interfaces/IState";
import { IDockContainer } from "./interfaces/IDockContainer";

/**
 * The document manager is then central area of the dock layout hierarchy.
 * This is where more important panels are placed (e.g. the text editor in an IDE,
 * 3D view in a modelling package etc
 */

export class DocumentManagerContainer extends FillDockContainer {

    constructor(dockManager: DockManager) {
        super(dockManager, TabHostDirection.TOP);

        this.minimumAllowedChildNodes = 0;
        this.element.classList.add('document-manager');
        this.tabHost.createTabPage = this._createDocumentTabPage;
        this.tabHost.displayCloseButton = true;
    }

    private _createDocumentTabPage(tabHost: TabHost, container: IDockContainer) {
        return new DocumentTabPage(tabHost, container);
    }

    saveState(state: IState) {
        super.saveState(state);
        state.documentManager = true;
    }

    /** Returns the selected document tab */
    selectedTab() {
        return this.tabHost.activeTab;
    }
}
