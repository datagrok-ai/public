import { DockManager } from "./DockManager";
import { DockModel } from "./DockModel";
import { DockNode } from "./DockNode";
import { PanelContainer } from "./PanelContainer";
import { HorizontalDockContainer } from "./HorizontalDockContainer";
import { VerticalDockContainer } from "./VerticalDockContainer";
import { DocumentManagerContainer } from "./DocumentManagerContainer";
import { FillDockContainer } from "./FillDockContainer";
import { Dialog } from "./Dialog";
import { Utils } from "./Utils";
import { IPanelInfo } from "./interfaces/IPanelInfo";
import { INodeInfo } from "./interfaces/INodeInfo";
import { IDockContainer } from "./interfaces/IDockContainer";

/**
 * Deserializes the dock layout hierarchy from JSON and creates a dock hierarhcy graph
 */
export class DockGraphDeserializer {

    dockManager: DockManager;
    documentManagerNode: DockNode;

    constructor(dockManager: DockManager) {
        this.dockManager = dockManager;
    }

    async deserialize(_json: string): Promise<DockModel> {
        let info = JSON.parse(_json);
        let model = new DockModel();
        model.rootNode = await this._buildGraph(info.graphInfo);
        model.dialogs = await this._buildDialogs(info.dialogsInfo);
        model.documentManagerNode = this.documentManagerNode;
        return model;
    }

    async _buildGraph(nodeInfo: INodeInfo) {
        let childrenInfo = nodeInfo.children;
        let children: DockNode[] = [];
        for (let childInfo of childrenInfo) {
            let childNode = await this._buildGraph(childInfo);
            if (childNode !== null) {
                children.push(childNode);
            }
        };

        // Build the container owned by this node
        let container = await this._createContainer(nodeInfo, children);
        if (container === null) {
            return null;
        }
        // Build the node for this container and attach it's children
        let node = new DockNode(container);
        if (container instanceof DocumentManagerContainer)
            this.documentManagerNode = node;
        node.children = children;
        for (let childNode of node.children.reverse()) {
            childNode.parent = node;
        };
        node.children.reverse();
        // node.container.setActiveChild(node.container);
        return node;
    }

    async _createContainer(nodeInfo: INodeInfo, children: DockNode[]) {
        let containerType = nodeInfo.containerType;
        let containerState = nodeInfo.state;
        let container;

        let childContainers: IDockContainer[] = [];
        for (let childNode of children) {
            childContainers.push(childNode.container);
        }

        if (containerType === 'panel') {
            container = await PanelContainer.loadFromState(containerState, this.dockManager);
            if (!container?.prepareForDocking)
                return null;
            container.prepareForDocking();
            Utils.removeNode(container.elementPanel);
        }
        else if (containerType === 'horizontal')
            container = new HorizontalDockContainer(this.dockManager, childContainers);
        else if (containerType === 'vertical')
            container = new VerticalDockContainer(this.dockManager, childContainers);
        else if (containerType === 'fill') {
            // Check if this is a document manager

            // TODO: Layout engine compares the string 'fill', so cannot create another subclass type
            // called document_manager and have to resort to this hack. use RTTI in layout engine
            let typeDocumentManager = containerState.documentManager;
            if (typeDocumentManager)
                container = new DocumentManagerContainer(this.dockManager);
            else
                container = new FillDockContainer(this.dockManager);
        }
        else
            throw new Error('Cannot create dock container of unknown type: ' + containerType);

        // Restore the state of the container
        container.loadState(containerState);

        // container.performLayout(childContainers);
        return container;
    }

    async _buildDialogs(dialogsInfo: IPanelInfo[]) {
        let dialogs: Dialog[] = [];
        for (let dialogInfo of dialogsInfo) {
            let containerType = dialogInfo.containerType;
            let containerState = dialogInfo.state;
            let container;
            if (containerType === 'panel') {
                container = await PanelContainer.loadFromState(containerState, this.dockManager);
                if (container.prepareForDocking) {
                    Utils.removeNode(container.elementPanel);
                    container.isDialog = true;
                    let dialog = new Dialog(container, this.dockManager);
                    if (dialogInfo.position.x > document.body.clientWidth ||
                        dialogInfo.position.y > document.body.clientHeight - 70) {
                        dialogInfo.position.x = 20;
                        dialogInfo.position.y = 70;
                    }
                    dialog.setPosition(dialogInfo.position.x, dialogInfo.position.y);
                    dialog.isHidden = dialogInfo.isHidden;
                    if (dialog.isHidden)
                        dialog.hide();
                    dialogs.push(dialog);
                }
            }

        }
        return dialogs;
    }
}
