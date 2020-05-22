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

};


export class DockContainer {
    constructor(d) { this.d = d; }
}


/** Window docking manager */
export class DockManager {
    constructor(d) { this.d = d; }

    /**
     * Docks the element relative to the reference node.
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
    dock(element, dockType, refNode, title = '', ratio = 0.5) {
        return new DockNode(grok_DockManager_Dock(this.d, refNode.d, element, dockType, title, ratio));
    }

    /**
     * Docks the element relative to the reference node.
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
    dockDialog(element, dockType, refNode, title = '') {
        return new DockNode(grok_DockManager_DockDialog(this.d, refNode.d, element, dockType, title));
    }
}