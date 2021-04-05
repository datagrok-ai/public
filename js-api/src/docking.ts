import {toJs} from './wrappers';
// @ts-ignore
import {_toIterable} from './utils';
import {DockType} from './const';
import { Viewer } from './viewer';


declare let DG: any;
let api = <any>window;

/**
 * Dock node.
 * See also {@link DockContainer}, {@link DockManager}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockNode {
  d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** @returns {DockContainer} */
  get container(): DockContainer {
    return new DockContainer(api.grok_DockNode_Get_Container(this.d));
  }

  /** Detaches this node from parent. */
  detachFromParent(): void {
    return api.grok_DockNode_DetachFromParent(this.d);
  }

  /** Removes a child node.
   * @param {DockNode} childNode */
  removeChild(childNode: DockNode): void {
    return api.grok_DockNode_RemoveChild(this.d, childNode);
  }

  /** @returns {DockNode} */
  get parent(): DockNode {
    return toJs(api.grok_DockNode_Parent(this.d));
  }

  /** @returns {Iterable.<DockNode>} */
  get children(): Iterable<DockNode> {
    return _toIterable(api.grok_DockNode_Children(this.d));
  }
}


/**
 * Represents a dockable window.
 * See also {@link DockContainer}, {@link DockNode}.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
 * */
export class DockContainer {
  d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Container element.
   * @returns {HTMLDivElement} */
  get containerElement(): HTMLDivElement {
    return api.grok_DockContainer_Get_ContainerElement(this.d);
  }

  /** Destroys and detaches the container. */
  destroy(): void {
    api.grok_DockContainer_Destroy(this.d);
  }

  /** Undocks a panel and converts it into a floating dialog window
   *  It is assumed that only leaf nodes (panels) can be undocked */
  float(): void {
    api.grok_DockContainer_Float(this.d);
  }

  /** Removes a dock container from the dock layout hierarchy
   *  @returns {DockNode} - the node that was removed from the dock tree */
  //remove() { return new DockNode(api.grok_DockContainer_Remove(this.d)); }
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
  d: any;

  constructor(d: any) {
    this.d = d;
  }

  get element(): HTMLDivElement {
    return api.grok_DockManager_Get_Element(this.d);
  }

  get rootNode(): DockNode {
    return toJs(api.grok_DockManager_Get_RootNode(this.d));
  }

  /**
   * The document view is then central area of the dock layout hierarchy.
   * This is where more important panels are placed (e.g. the text editor in an IDE,
   * 3D view in a modeling package, etc.)
   */
  get documentContainer(): DockContainer {
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
  dock(element: HTMLElement | Viewer, dockType: DockType = DG.DOCK_TYPE.LEFT, refNode: DockNode | null = null, title?: string, ratio: number = 0.5): DockNode {
    return new DockNode(api.grok_DockManager_Dock(this.d, refNode === null ? null : refNode.d, element, dockType, title, ratio));
  }

  /**
   * Undocks the element.
   * @param {HTMLElement | DockNode} object - Element to undock
   * */
  close(object: HTMLElement | DockNode): void {
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
  findNode(element: HTMLElement): DockNode {
    return toJs(api.grok_DockManager_FindNode(this.d, element));
  }

  // /**
  //  * Docks the element relative to the reference node.
  //  * @param {DockType} dockType - Dock type (left | right | top | bottom | fill).
  //  * @param {number} ratio - Ratio of the area to take (relative to the reference node).
  //  * @param {string=} title - Name of the resulting column. Default value is agg(colName).
  //  * @returns {DockNode}
  //  * */
  // dockDialog(element, dockType, refNode, title = '') {
  //     return new DockNode(api.grok_DockManager_DockDialog(this.d, refNode == null ? null : refNode.d, element, dockType, title));
  // }
}