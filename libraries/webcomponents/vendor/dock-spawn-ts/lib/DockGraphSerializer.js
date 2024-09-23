/**
 * The serializer saves / loads the state of the dock layout hierarchy
 */
export class DockGraphSerializer {
    serialize(model) {
        const graphInfo = this._buildGraphInfo(model.rootNode);
        const dialogs = this._buildDialogsInfo(model.dialogs.sort((x, y) => x.elementDialog.style.zIndex - y.elementDialog.style.zIndex));
        return JSON.stringify({ graphInfo: graphInfo, dialogsInfo: dialogs });
    }
    _buildGraphInfo(node) {
        const nodeState = {};
        node.container.saveState(nodeState);
        const childrenInfo = [];
        node.children.forEach((childNode) => {
            childrenInfo.push(this._buildGraphInfo(childNode));
        });
        const nodeInfo = {
            containerType: node.container.containerType,
            state: nodeState,
            children: childrenInfo,
        };
        return nodeInfo;
    }
    _buildDialogsInfo(dialogs) {
        const dialogsInfo = [];
        dialogs.forEach((dialog) => {
            const panelState = {};
            const panelContainer = dialog.panel;
            panelContainer.saveState(panelState);
            const panelInfo = {
                containerType: panelContainer.containerType,
                state: panelState,
                position: dialog.getPosition(),
                isHidden: dialog.isHidden,
            };
            dialogsInfo.push(panelInfo);
        });
        return dialogsInfo;
    }
}
