/**
 * @typedef {string} DockType
 **/

/** @enum {DockType} */
export const DOCK_TYPE = {
    LEFT: "left",
    RIGHT: "right",
    TOP: "top",
    BOTTOM: "bottom",
    FILL: "fill",
};


export class DockNode {
    constructor(d) { this.d = d; }

    /** @returns {DockContainer} */
    get container() { return new DockContainer(grok_DockNode_Get_Container(this.d)); }

    detachFromParent() { return grok_DockNode_DetachFromParent(this.d); }
};


export class DockContainer {
    constructor(d) { this.d = d; }

    get element() { return grok_DockContainer_Get_Element(this.d); }
}


/** Window docking manager */
export class DockManager {
    constructor(d) { this.d = d; }

    get element() { return grok_DockManager_Get_Element(this.d); }

    /**
     * The document view is then central area of the dock layout hierarchy.
     * This is where more important panels are placed (e.g. the text editor in an IDE,
     * 3D view in a modelling package etc
     */
    get documentContainer() { return new DockContainer(grok_DockManaget_Get_DocumentView(this.d)); }

    /**
     * Docks the element relative to the reference node.
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
    dock(element, dockType, refNode, title = '', ratio = 0.5) {
        return new DockNode(grok_DockManager_Dock(this.d, refNode == null ? null : refNode.d, element, dockType, title, ratio));
    }

    /**
     * Docks the element relative to the reference node.
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
    dockDialog(element, dockType, refNode, title = '') {
        return new DockNode(grok_DockManager_DockDialog(this.d, refNode == null ? null : refNode.d, element, dockType, title));
    }
}