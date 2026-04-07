/**
 * Tree view and tag editor classes: TagEditor, TagElement, TreeViewNode, TreeViewGroup.
 * @module widgets/tree
 */

import {toDart, toJs} from "../wrappers";
import {Observable} from "rxjs";
import {__obs, _sub, StreamSubscription} from "../events";
import {IDartApi} from "../api/grok_api.g";
import {HttpDataSource} from "../dapi";

import '../../css/tags-input.css';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class TagEditor {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(): TagEditor {
    return toJs(api.grok_TagEditor());
  }

  get root(): HTMLElement {
    return api.grok_TagEditor_Get_Root(this.dart);
  }

  get tags(): TagElement[] {
    return api.grok_TagEditor_Get_Tags(this.dart);
  }

  addTag(tag: TagElement | string, notify: boolean = true) {
    return api.grok_TagEditor_AddTag(this.dart, tag, notify);
  }

  removeTag(tag: TagElement | string): void {
    api.grok_TagEditor_RemoveTag(this.dart, tag);
  }

  clearTags(): void {
    api.grok_TagEditor_ClearTags(this.dart);
  }

  set acceptsDragDrop(predicate: (...params: any[]) => boolean) {
    api.grok_TagEditor_Set_AcceptsDragDrop(this.dart, (x: any) => predicate(toJs(x, false)));
  };

  set doDrop(action: Function) {
    api.grok_TagEditor_Set_DoDrop(this.dart, (x: any) => action(toJs(x, false)));
  }

  onChanged(callback: Function): StreamSubscription {
    return _sub(api.grok_TagEditor_OnChanged(this.dart, callback));
  }
}


export class TagElement {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  get tag(): TagElement | string {
    return api.grok_TagElement_Get_Tag(this.dart);
  };

  set tag(x: TagElement | string) {
    api.grok_TagElement_Set_Tag(this.dart, x);
  };
}


/** Tree view node.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * */
export class TreeViewNode<T = any> {
  dart: any;

  /** @constructs {TreeView} from the Dart object */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Visual root */
  get root(): HTMLElement {
    return api.grok_TreeViewNode_Root(this.dart);
  }

  /** Top-most node. */
  get rootNode(): TreeViewGroup {
    let x: TreeViewNode = this;
    while (x.parent)
      x = x.parent;
    return x as TreeViewGroup;
  }

  /** Node's parent */
  get parent(): TreeViewNode { return toJs(api.grok_TreeViewNode_Parent(this.dart)); }

  /** Caption label */
  get captionLabel(): HTMLElement { return api.grok_TreeViewNode_CaptionLabel(this.dart); }

  /** Check box element  */
  get checkBox(): HTMLElement | null { return api.grok_TreeViewNode_CheckBox(this.dart); }

  /** Returns `true` if checked */
  get checked(): boolean { return api.grok_TreeViewNode_Get_Checked(this.dart); }
  set checked(checked: boolean) { api.grok_TreeViewNode_Set_Checked(this.dart, checked); }

  /** Node text */
  get text(): string { return api.grok_TreeViewNode_Get_Text(this.dart); }
  set text(value: string) {api.grok_TreeViewNode_Set_Text(this.dart, value); }

  /** Node icon */
  get icon(): Element { return api.grok_TreeViewNode_Get_Icon(this.dart); }
  set icon(value: Element) {api.grok_TreeViewNode_Set_Icon(this.dart, value); }

  /** Auxiliary information associated with the node. */
  get tag(): any { return api.grok_TreeViewNode_Get_Tag(this.dart); }
  set tag(t : any) { api.grok_TreeViewNode_Set_Tag(this.dart, t); }

  /** Node value. Normally, when you click on the node, the context panel shows this object. */
  get value(): T { return api.grok_TreeViewNode_Get_Value(this.dart); };
  set value(v: T) { api.grok_TreeViewNode_Set_Value(this.dart, v)};

  /** Enables checkbox */
  enableCheckBox(checked: boolean = false): void {
    api.grok_TreeViewNode_EnableCheckBox(this.dart, checked);
  }

  /** Occurs when the selected node is changed. */
  get onSelected(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-current', this.dart); }

  /** Removes the node and its children from the parent */
  remove(): void {
    api.grok_TreeViewNode_Remove(this.dart);
  }
}

export class TreeViewGroup extends TreeViewNode {
  /** Creates new tree */
  static tree(): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_Tree());
  }

  static fromItemCategories(items: any[], props: string[], options?: {
    removeEmpty: boolean, itemToElement?: (item:any) => Element, itemToString?: (item: any) => string, itemToValue?: (item: any) => any
  }): TreeViewGroup {

    function init(node: TreeViewGroup, path: string[]) {

      //
      const pathItems = items
          .filter((item) => props.every((p: string) => !(options?.removeEmpty ?? false) || item[p] != null))
          .filter((item) => path.every((p, i) => item[props[i]] == path[i]));

      // leafs
      if (path.length == props.length) {
        let itemToValue = options!.itemToValue ?? function (i) {return i;};
        let itemToString = options!.itemToElement ?? options!.itemToString;
        for (let item of pathItems)
          node.item(itemToString!(item), itemToValue(item));
      }
      else {
        let categories: Set<string> = new Set<string>();
        for (let item of pathItems)
          categories.add(item[props[path.length]]);

        for (let category of categories)
          init(node.group(category), [...path, category]);

        // lazy loading - use it for big collections and async queries
        // for (let category of categories)
        //   node.group(category)
        //     .onNodeExpanding
        //     .subscribe((expandingNode) => init(expandingNode, [...path, category]));
      }
    }

    let rootNode = TreeViewGroup.tree();
    init(rootNode, []);
    return rootNode;
  }

  /** Gets all node items */
  get items(): TreeViewNode[] {
    return api.grok_TreeViewNode_Items(this.dart).map((i: any) => toJs(i));
  }

  /** Gets the node's children */
  get children(): TreeViewNode[] {
    return api.grok_TreeViewNode_Children(this.dart).map((i: any) => toJs(i));
  }

  /** Removes all children (going down recursively) that satisfy the predicate */
  removeChildrenWhere(predicate: (node: TreeViewNode) => boolean): void {
    for (const child of this.children) {
      if (predicate(child)) {
        child.remove();
      }
      else if (child instanceof TreeViewGroup) {
        child.removeChildrenWhere(predicate);
      }
    }
  }

  /** Controls expanded state. */
  get expanded(): boolean { return api.grok_TreeViewNode_Get_Expanded(this.dart); }
  set expanded(isExpanded: boolean) { api.grok_TreeViewNode_Set_Expanded(this.dart, isExpanded); }

  /** Indicates whether check or uncheck is applied to a node only or to all node's children */
  get autoCheckChildren(): boolean { return api.grok_TreeViewNode_GetAutoCheckChildren(this.dart); }
  set autoCheckChildren(auto: boolean) { api.grok_TreeViewNode_SetAutoCheckChildren(this.dart, auto); }

  /** Currently selected node. */
  get currentItem(): TreeViewNode { return toJs(api.grok_TreeViewNode_Get_CurrentItem(this.dart)); }
  set currentItem(node: TreeViewNode) { api.grok_TreeViewNode_Set_CurrentItem(this.dart, toDart(node)); }

  /** Adds new group and returns it */
  group(text: string | Element, value: object | null = null, expanded: boolean = true, index: number | null = null): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_Group(this.dart, text, value, expanded, index));
  }

  /** Returns existing, or creates a new node group */
  getOrCreateGroup(text: string, value: object | null = null, expanded: boolean = true): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_GetOrCreateGroup(this.dart, text, value, expanded));
  }

  /** Adds new item to group */
  item(text: string | Element, value: object | null = null): TreeViewNode {
    return toJs(api.grok_TreeViewNode_Item(this.dart, text, value));
  }

  /** Adds new items to group */
  addItems(items: any[]): TreeViewNode {
    for (const i of items)
      this.item(i, i);
    return this;
  }


  get onNodeExpanding(): Observable<TreeViewGroup> { return __obs('d4-tree-view-node-expanding', this.dart); }
  get onNodeAdded(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-added', this.dart); }
  get onNodeCheckBoxToggled(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-checkbox-toggled', this.dart); }
  get onChildNodeExpandedChanged(): Observable<TreeViewGroup> { return __obs('d4-tree-view-child-node-expanded-changed', this.dart); }
  get onChildNodeExpanding(): Observable<TreeViewGroup> { return __obs('d4-tree-view-child-node-expanding', this.dart); }
  // get onChildNodeContextMenu(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-context-menu', this.dart); }
  get onNodeContextMenu(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-context-menu', this.dart); }
  get onSelectedNodeChanged(): Observable<TreeViewNode> { return __obs('d4-tree-view-selected-node-changed', this.dart); }
  get onNodeMouseEnter(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-mouse-enter', this.dart); }
  get onNodeMouseLeave(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-mouse-leave', this.dart); }
  get onNodeEnter(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-enter', this.dart); }

  async loadSources(source: HttpDataSource<any>): Promise<void> {
    return api.grok_TreeViewGroup_Load_Sources(this.dart, source.dart);
  }
}
