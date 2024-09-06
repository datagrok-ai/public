import {DockModel} from './DockModel.js';
import {DockNode} from './DockNode.js';
import {Dialog} from './Dialog.js';
import {IPanelInfo} from './interfaces/IPanelInfo.js';
import {INodeInfo} from './interfaces/INodeInfo.js';
import {IState} from './interfaces/IState.js';

/**
 * The serializer saves / loads the state of the dock layout hierarchy
 */
export class DockGraphSerializer {
  serialize(model: DockModel) {
    const graphInfo = this._buildGraphInfo(model.rootNode);
    const dialogs = this._buildDialogsInfo(model.dialogs.sort((x, y)=><number><any>x.elementDialog.style.zIndex-<number><any>y.elementDialog.style.zIndex));
    return JSON.stringify({graphInfo: graphInfo, dialogsInfo: dialogs});
  }

  _buildGraphInfo(node: DockNode): INodeInfo {
    const nodeState: IState = {};
    node.container.saveState(nodeState);

    const childrenInfo: INodeInfo[] = [];
    node.children.forEach((childNode) => {
      childrenInfo.push(this._buildGraphInfo(childNode));
    });

    const nodeInfo: INodeInfo = {
      containerType: node.container.containerType,
      state: nodeState,
      children: childrenInfo,
    };
    return nodeInfo;
  }

  _buildDialogsInfo(dialogs: Dialog[]): IPanelInfo[] {
    const dialogsInfo: IPanelInfo[] = [];
    dialogs.forEach((dialog) => {
      const panelState: IState = {};
      const panelContainer = dialog.panel;
      panelContainer.saveState(panelState);

      const panelInfo: IPanelInfo = {
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
