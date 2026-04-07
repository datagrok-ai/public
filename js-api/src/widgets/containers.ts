/**
 * Container widgets: Accordion, AccordionPane, TabControl, TabPane, ToolboxPage.
 * @module widgets/containers
 */

import {toJs} from "../wrappers";
import {Observable} from "rxjs";
import {__obs} from "../events";
import {IDartApi} from "../api/grok_api.g";
import {DartWidget} from "./base";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/**
 * Accordion control with collapsible/expandable panes.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/accordion}
 * @extends {DartWidget}
 * */
export class Accordion extends DartWidget {

  /** @constructs Accordion */
  constructor(dart: any) {
    super(dart);
  }

  /** An object this accordion is associated with */
  get context(): any {
    return toJs(api.grok_Accordion_Get_Context(this.dart));
  }
  set context(x: any) {
    api.grok_Accordion_Set_Context(this.dart, x);
  }

  /** Creates a new instance of Accordion */
  static create(key: any = null): Accordion {
    return toJs(api.grok_Accordion(key));
  }

  /** @type {AccordionPane[]} */
  get panes(): AccordionPane[] {
    return api.grok_TabControlBase_Get_Panes(this.dart).map(toJs);
  }

  /** Header element on top of the accordion */
  get header(): HTMLElement { return api.grok_Accordion_Get_Header(this.dart); }
  set header(header) { api.grok_Accordion_Set_Header(this.dart, header); }

  /** Whether tab header should be hidden if there is only one tab */
  get autoHideTabHeader(): boolean { return api.grok_Accordion_Get_AutoHideTabHeader(this.dart); }
  set autoHideTabHeader(x) { api.grok_Accordion_Set_AutoHideTabHeader(this.dart, x); }

  /** Returns a pane with the specified name.
   * @param {string} name
   * @returns {AccordionPane} */
  getPane(name: string): AccordionPane {
    return toJs(api.grok_TabControlBase_GetPane(this.dart, name));
  }

  /** Adds a title element. */
  addTitle(element: HTMLElement): void {
    return api.grok_Accordion_AddTitle(this.dart, element);
  }

  /** Adds a pane */
  addPane(name: string, getContent: () => HTMLElement, expanded: boolean = false, before: AccordionPane | null = null,
    allowDragOut: boolean = true): AccordionPane {
    return toJs(api.grok_Accordion_AddPane(this.dart, name, getContent, expanded, before !== null ? before.dart : null, null, allowDragOut));
  }

  /** Adds a pane with the count indicator next to the title.
   * getCount() is executed immediately. */
  addCountPane(name: string, getContent: () => HTMLElement, getCount: () => number, expanded: boolean = false, before: AccordionPane | null = null,
    allowDragOut: boolean = true): AccordionPane {
    return toJs(api.grok_Accordion_AddPane(this.dart, name, getContent, expanded, before !== null ? before.dart : null, getCount, allowDragOut));
  }

  /** Removed the specified pane. */
  removePane(pane: AccordionPane) {
    api.grok_Accordion_RemovePane(this.dart, pane.dart);
  }

  /** Finalizes accordion construction */
  end() {
    api.grok_Accordion_End(this.dart);
  }
}


/** A pane in the {@link Accordion} control. */
export class AccordionPane extends DartWidget {
  declare dart: any;

  constructor(dart: any) {
    super(dart);
  }

  /** Expanded state
   * @type {boolean} */
  get expanded(): boolean {
    return api.grok_AccordionPane_Get_Expanded(this.dart);
  }

  set expanded(v: boolean) {
    api.grok_AccordionPane_Set_Expanded(this.dart, v);
  }

  /** @type {string} */
  get name(): string {
    return api.grok_AccordionPane_Get_Name(this.dart);
  }

  set name(name: string) {
    api.grok_AccordionPane_Set_Name(this.dart, name);
  }
}


/** Tab control that hosts panes inside. See also {@link TabPane} */
export class TabControl extends DartWidget {

  constructor(dart: any) {
    super(dart);
    this.dart = dart;
  }

  /** Creates a new TabControl */
  static create(vertical: boolean = false): TabControl {
    return toJs(api.grok_TabControl(vertical));
  }

  /** Visual root */
  get root(): HTMLDivElement {
    return api.grok_Widget_Get_Root(this.dart);
  }

  /** Header shown on top of the control */
  get header(): HTMLDivElement {
    return api.grok_TabControlBase_Get_Header(this.dart);
  }

  /** Panes currently present in the pane control.
   * Do not change the array, use {@link addPane} instead */
  get panes(): TabPane[] {
    return api.grok_TabControlBase_Get_Panes(this.dart).map(toJs);
  }

  /** Gets the pane with the specified name */
  getPane(name: string): TabPane {
    return toJs(api.grok_TabControlBase_GetPane(this.dart, name));
  }

  /** Adds a new pane with the specified name */
  addPane(name: string, getContent: () => HTMLElement, icon: any = null, options?: {allowClose: boolean}): TabPane {
    return toJs(api.grok_TabControlBase_AddPane(this.dart, name, getContent, icon, options?.allowClose ?? false));
  }

  /** Removes all panes */
  clear(): void {
    api.grok_TabControlBase_Clear(this.dart);
  }

  /** Currently visible pane */
  get currentPane(): TabPane { return api.grok_TabControlBase_Get_CurrentPane(this.dart); }
  set currentPane(v: TabPane) { api.grok_TabControlBase_Set_CurrentPane(this.dart, v.dart); }

  /** Occurs before the active pane is changed */
  get onBeforeTabChanged(): Observable<any> { return __obs('d4-tabcontrol-before-tab-changed', this.dart); }

  /** Occurs after the active pane is changed */
  get onTabChanged(): Observable<any> { return __obs('d4-tabcontrol-tab-changed', this.dart); }

  /** Occurs after a pane is added */
  get onTabAdded(): Observable<any> { return __obs('d4-tabcontrol-tab-added', this.dart); }

  /** Occurs after a pane is removed */
  get onTabRemoved(): Observable<any> { return __obs('d4-tabcontrol-tab-removed', this.dart); }
}


/** Represents a pane of either {@link TabControl} or {@link Accordion} */
export class TabPane extends DartWidget {

  /** Creates TabPane from the Dart handle */
  constructor(dart: any) {
    super(dart);
  }

  /** {@link TabControl} this pane belongs to */
  get parent(): TabControl {
    return toJs(api.grok_TabPane_Get_Parent(this.dart));
  }

  /** A control shown on top of the pane */
  get header(): HTMLDivElement {
    return api.grok_TabPane_Get_Header(this.dart);
  }

  /** Content */
  get content(): HTMLDivElement {
    return api.grok_TabPane_Get_Content(this.dart);
  }

  /** Whether the pane is expanded. Applicable to Accordion's panes only. */
  get expanded(): boolean { return api.grok_AccordionPane_Get_Expanded(this.dart); }
  set expanded(v: boolean) { api.grok_AccordionPane_Set_Expanded(this.dart, v); }

  /** Tab pane name */
  get name(): string { return api.grok_AccordionPane_Get_Name(this.dart); }
  set name(name: string) { api.grok_AccordionPane_Set_Name(this.dart, name); }
}


export class ToolboxPage {
  dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }

  get accordion(): Accordion {
    return toJs(api.grok_ToolboxPage_Get_Accordion(this.dart));
  }
}
