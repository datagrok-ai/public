import { TabHost } from "./TabHost";
import { TabPage } from "./TabPage";
import { Utils } from "./Utils";
import { IDockContainer } from "./interfaces/IDockContainer";
import { PanelContainer } from "./PanelContainer";

/**
 * Specialized tab page that doesn't display the panel's frame when docked in a tab page
 */
export class DocumentTabPage extends TabPage {

    constructor(host: TabHost, container: IDockContainer) {
        super(host, container);

        // If the container is a panel, extract the content element and set it as the tab's content
        if (this.container.containerType === 'panel') {
            this.panel = container as PanelContainer;
            this.containerElement = this.panel.elementContentWrapper;

            // detach the container element from the panel's frame.
            // It will be reattached when this tab page is destroyed
            // This enables the panel's frame (title bar etc) to be hidden
            // inside the tab page
            Utils.removeNode(this.containerElement);
        }
    }

    destroy() {
        super.destroy();

        // Restore the panel content element back into the panel frame
        //Utils.removeNode(this.containerElement);
        this.panel.elementContentHost.appendChild(this.containerElement);
    }
}
