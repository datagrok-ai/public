import { toJs } from './wrappers';
import { _toIterable } from './utils';
let api = window;
/**
 * Dock node.
 * See also {@link DockContainer}, {@link DockManager}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockNode {
  constructor(d) {
    this.d = d;
  }
  /** @returns {DockContainer} */
  get container() {
    return new DockContainer(api.grok_DockNode_Get_Container(this.d));
  }
  /** Detaches this node from parent. */
  detachFromParent() {
    return api.grok_DockNode_DetachFromParent(this.d);
  }
  /** Removes a child node.
     * @param {DockNode} childNode */
  removeChild(childNode) {
    return api.grok_DockNode_RemoveChild(this.d, childNode);
  }
  /** @returns {DockNode} */
  get parent() {
    return toJs(api.grok_DockNode_Parent(this.d));
  }
  /** @returns {Iterable.<DockNode>} */
  get children() {
    return _toIterable(api.grok_DockNode_Children(this.d));
  }
}
/**
 * Represents a dockable window.
 * See also {@link DockContainer}, {@link DockNode}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockContainer {
  constructor(d) {
    this.d = d;
  }
  /** Container element.
     * @returns {HTMLDivElement} */
  get containerElement() {
    return api.grok_DockContainer_Get_ContainerElement(this.d);
  }
  /** Destroys and detaches the container. */
  destroy() {
    api.grok_DockContainer_Destroy(this.d);
  }
  /** Undocks a panel and converts it into a floating dialog window
     *  It is assumed that only leaf nodes (panels) can be undocked */
  float() {
    api.grok_DockContainer_Float(this.d);
  }
}
/**
 * Window docking manager.
 *
 * Dock manager manages all the dock panels in a hierarchy, similar to visual studio.
 * It owns an HTML Div element inside which all panels are docked.
 * Initially, the dock manager takes up the central space and acts as the root node.
 *
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * Learn more: {@link https://github.com/coderespawn/dock-spawn} for details.
 */
export class DockManager {
  constructor(d) {
    this.d = d;
  }
  get element() {
    return api.grok_DockManager_Get_Element(this.d);
  }
  get rootNode() {
    return toJs(api.grok_DockManager_Get_RootNode(this.d));
  }
  /**
     * The document view is then central area of the dock layout hierarchy.
     * This is where more important panels are placed (e.g. the text editor in an IDE,
     * 3D view in a modeling package, etc.)
     */
  get documentContainer() {
    return new DockContainer(api.grok_DockManaget_Get_DocumentContainer(this.d));
  }
  /**
     * Docks the element relative to the reference node.
     * @param {HTMLElement | Viewer} element - Element to dock
     * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
     * @param {DockNode|null} refNode - reference node
     * @param {number} ratio - Ratio of the area to take (relative to the reference node).
     * @param {string=} title - Name of the resulting column. Default value is agg(colName).
     * @returns {DockNode}
     * */
  dock(element, dockType = DG.DOCK_TYPE.LEFT, refNode = null, title, ratio = 0.5) {
    return new DockNode(api.grok_DockManager_Dock(this.d, refNode === null ? null : refNode.d, element, dockType, title, ratio));
  }
  /**
     * Undocks the element.
     * @param {HTMLElement | DockNode} object - Element to undock
     * */
  close(object) {
    // @ts-ignore
    if (object.d === undefined)
      api.grok_DockManager_UndockByElement(this.d, object);
    else
    // @ts-ignore
      api.grok_DockManager_UndockNode(this.d, object.d);
  }
  /**
     * Finds the node of an element.
     * @param {HTMLElement} element - Element to find the node for.
     * @returns {DockNode}
     * */
  findNode(element) {
    return toJs(api.grok_DockManager_FindNode(this.d, element));
  }
}
