import {DockType} from "./const";

/**
 * Represents a dockable window.
 * See also {@link DockContainer}, {@link DockNode}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockNode {

    /** @returns {DockContainer} */
    get container(): DockContainer

    /** Detaches this node from parent. */
    detachFromParent(): void

    /** Removes a child node.
     * @param {DockNode} childNode */
    removeChild(childNode: DockNode): void
}


/**
 * Represents a dockable window.
 * See also {@link DockContainer}, {@link DockNode}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockContainer {

    /** Container element.
     * @returns {HTMLDivElement} */
    get containerElement(): HTMLDivElement;

    /** Destroys and detaches the container. */
    destroy(): void;

    /** Undocks a panel and converts it into a floating dialog window
     *  It is assumed that only leaf nodes (panels) can be undocked */
    float(): void

    /** Removes a dock container from the dock layout hierarcy
     *  @returns {DockNode} - the node that was removed from the dock tree */
    //remove() { return new DockNode(grok_DockContainer_Remove(this.d)); }
}


/**
 * Window docking manager.
 *
 * Dock manager manages all the dock panels in a hierarchy, similar to visual studio.
 * It owns a Html Div element inside which all panels are docked
 * Initially the document manager takes up the central space and acts as the root node
 *
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * Learn more: {@link https://github.com/coderespawn/dock-spawn} for details.
 */
export class DockManager {

    get element(): HTMLDivElement

    get rootNode(): DockNode

    /**
     * The document view is then central area of the dock layout hierarchy.
     * This is where more important panels are placed (e.g. the text editor in an IDE,
     * 3D view in a modeling package etc
     */
    get documentContainer(): DockContainer

    /**
     * Docks the element relative to the reference node.
     * @param {HTMLElement} element - Element to dock
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {DockNode|null} refNode - reference node
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
    dock(element: HTMLHtmlElement, dockType?: DockType, refNode?: DockNode | null, title?: string, ratio?: number): DockNode

    // /**
    //  * Docks the element relative to the reference node.
    //  * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
    //  * @param {number} ratio - Ratio of the area to take (relative to the reference node).
    //  * @param {string=} title - Name of the resulting column. Default value is agg(colName).
    //  * @returns {DockNode}
    //  * */
    // dockDialog(element, dockType, refNode, title = '') {
    //     return new DockNode(grok_DockManager_DockDialog(this.d, refNode == null ? null : refNode.d, element, dockType, title));
    // }
}