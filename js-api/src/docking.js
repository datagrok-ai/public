import {toJs} from "./wrappers";
import {_toIterable} from "./utils";
import {DockType} from "./const";


/**
 * Represents a dockable window.
 * See also {@link DockContainer}, {@link DockNode}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockNode {
  constructor(d) {
    this.d = d;
  }

  /** @returns {DockContainer} */
  get container() {
    return new DockContainer(grok_DockNode_Get_Container(this.d));
  }

  /** Detaches this node from parent. */
  detachFromParent() {
    return grok_DockNode_DetachFromParent(this.d);
  }

  /** Removes a child node.
   * @param {DockNode} childNode */
  removeChild(childNode) {
    return grok_DockNode_RemoveChild(this.d, childNode);
  }

  /** @returns {DockNode} */
  get parent() {
    return toJs(grok_DockNode_Parent(this.d));
  }

  /** @returns {Iterable.<DockNode>} */
  get children() {
    return _toIterable(grok_DockNode_Children(this.d));
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
    return grok_DockContainer_Get_ContainerElement(this.d);
  }

  /** Destroys and detaches the container. */
  destroy() {
    grok_DockContainer_Destroy(this.d);
  }

  /** Undocks a panel and converts it into a floating dialog window
   *  It is assumed that only leaf nodes (panels) can be undocked */
  float() {
    grok_DockContainer_Float(this.d);
  }

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
  constructor(d) {
    this.d = d;
  }

  get element() {
    return grok_DockManager_Get_Element(this.d);
  }

  get rootNode() {
    return toJs(grok_DockManager_Get_RootNode(this.d));
  }

  /**
   * The document view is then central area of the dock layout hierarchy.
   * This is where more important panels are placed (e.g. the text editor in an IDE,
   * 3D view in a modeling package etc
   */
  get documentContainer() {
    return new DockContainer(grok_DockManaget_Get_DocumentContainer(this.d));
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
    return new DockNode(grok_DockManager_Dock(this.d, refNode === null ? null : refNode.d, element, dockType, title, ratio));
  }

  /**
   * Undocks the element.
   * @param {HTMLElement | DockNode} object - Element to undock
   * */
  close(object) {
    if (object.d === undefined)
      grok_DockManager_UndockByElement(this.d, object);
    else
      grok_DockManager_UndockNode(this.d, object.d);
  }

  /**
   * Finds node of the element.
   * @param {HTMLElement} element - Element to find
   * @returns {DockNode}
   * */
  findNode(element) {
    return toJs(grok_DockManager_FindNode(this.d, element));
  }

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