/**
 * Menu and Balloon classes.
 * @module widgets/menu
 */

import {toDart, toJs} from "../wrappers";
import {Observable} from "rxjs";
import {__obs, EventData} from "../events";
import {DataFrame} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import {
  IMenuItemsOptions,
  IMenuColorPaletteOptions,
  IMenuFontEditorOptions,
  IMenuSingleColumnSelectorOptions,
  IMenuMultiColumnSelectorOptions,
  IMenuHeaderOptions,
  IMenuItemOptions,
  IShowMenuOptions
} from "./types";

declare let ui: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/**
 * Menu (either top menu or popup menu).
 * Top menu sample: {@link https://public.datagrok.ai/js/samples/ui/menu}
 * Popup menu sample: {@link https://public.datagrok.ai/js/samples/ui/popup-menu}
 *
 * @example
 * DG.Menu.popup()
 *   .item('Show info', () => grok.shell.info('Info'))
 *   .separator()
 *   .items(['First', 'Second'], (s) => grok.shell.info(s))
 *   .show();
 * */
export class Menu {
  dart: any;
  _check: HTMLDivElement = ui.div('d4-menu-item-check');

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Creates a top menu. */
  static create(): Menu {
    return toJs(api.grok_Menu());
  }

  /** Creates a popup menu. */
  static popup(): Menu {
    return toJs(api.grok_Menu_Context());
  }

  /** Visual root */
  get root(): HTMLElement { return api.grok_Menu_Get_Root(this.dart); }

  /** Whether the menu closes when clicked. */
  get closeOnClick(): boolean { return api.grok_Menu_Get_CloseOnClick(this.dart); }
  set closeOnClick(value: boolean) { api.grok_Menu_Set_CloseOnClick(this.dart, value); }

  /** Finds a child menu item with the specified text. */
  find(text: string): Menu {
    return toJs(api.grok_Menu_Find(this.dart, text));
  }

  /** Executes the onClick function for that menu item.
   * Only works for items, not groups. */
  click(): void {
    api.grok_Menu_Click(this.dart);
  }

  /** Removes a child menu item with the specified text. */
  remove(text: string): void {
    api.grok_Menu_Remove(this.dart, text);
  }

  /** Removes all child menu items. */
  clear(): void {
    api.grok_Menu_Clear(this.dart);
  }

  /** Returns an existing menu group or adds a new group with the specified text. */
  group(text: string, order: number | null = null): Menu {
    return toJs(api.grok_Menu_Group(this.dart, text, order));
  }

  /** Ends a group of menu items and returns to the higher menu level.
   * @returns {Menu} `this` menu itself. */
  endGroup(): Menu {
    return toJs(api.grok_Menu_EndGroup(this.dart));
  }

  /** Adds a menu group with the specified text and handler. */
  item(text: string, onClick: () => void, order: number | null = null, options: IMenuItemOptions | null = null): Menu {
    return toJs(api.grok_Menu_Item(this.dart, text, onClick, order, options));
  }

  /** For each item in items, adds a menu group with the specified text and handler. */
  items<T = any>(items: T[], onClick: (item: T) => void, options: IMenuItemsOptions<T> | null = null): Menu {
    return toJs(api.grok_Menu_Items(this.dart, items, onClick, options?.isValid, options?.isChecked,
      options?.hasOwnProperty('toString') ? options.toString : null,
      options?.getTooltip, options?.onMouseEnter, options?.radioGroup));
  }

  /** Adds color palettes colors to menu.
   * @param colors - Array of arrays of color choices.
   * @param options - Optional params and functions, see {@link IMenuColorPaletteOptions}.
   * @returns {Menu} `this` menu itself. */
  colorPalette(colors: number[][], options?: IMenuColorPaletteOptions): Menu {
    return toJs(api.grok_Menu_ColorPalette(this.dart, colors, options?.getInitialValue, options?.onSelect,
      options?.onPreview, options?.asGroup, options?.visible, options?.categorical ?? false, options?.resetColorMs ?? 200, options?.closeOnClick ?? true));
  }

  /** Adds font editor to menu.
   * @param initial - Initial font to be set first or reset.
   * @param options - Optional params and functions, see {@link IMenuFontEditorOptions}.
   * @returns {Menu} `this` menu itself. */
  fontEditor(initial: string, options?: IMenuFontEditorOptions): Menu {
    return toJs(api.grok_Menu_FontEditor(this.dart, initial, options?.fontSizeMin, options?.fontSizeMax,
      options?.fontSizeStep ?? 1, options?.fontFamilies, options?.asGroup, options?.onChange));
  }

  /** Adds single-column selector to menu.
   * @param dataFrame - Data frame to be used for the selector,where column choices are taken from.
   * @param options - Optional params and functions, see {@link IMenuSingleColumnSelectorOptions}.
   * @returns {Menu} `this` menu itself. */
  singleColumnSelector(dataFrame: DataFrame, options?: IMenuSingleColumnSelectorOptions): Menu {
    return toJs(api.grok_Menu_SingleColumSelector(this.dart, dataFrame.dart, options?.initialValue,
      !options?.onChange ? null : (grid: any, c: any, currentRowChanged: boolean) => options?.onChange?.(toJs(grid), toJs(c), currentRowChanged),
      options?.asGroup, options?.nullable ?? false, options?.visible ?? true, options?.editable ?? false, options?.closeOnClick ?? false,
      options?.changeOnHover ?? true, (c: any) => options?.columnFilter?.(toJs(c)) ?? true));
  }

  /** Adds multi-column selector to menu.
   * @param dataFrame - Data frame to be used for the selector,where column choices are taken from.
   * @param options - Optional params and functions, see {@link IMenuMultiColumnSelectorOptions}.
   * @returns {Menu} `this` menu itself. */
  multiColumnSelector(dataFrame: DataFrame, options?: IMenuMultiColumnSelectorOptions): Menu {
    return toJs(api.grok_Menu_MultiColumSelector(this.dart, dataFrame.dart, options?.initialValue,
      !options?.onChange ? null : (grid: any) => options?.onChange?.(toJs(grid)),
      options?.asGroup, options?.visible ?? true, options?.editable ?? false, (c: any) => options?.columnFilter?.(toJs(c)) ?? true));
  }

  /** Adds a header title.
   * @param text - Header title text.
   * @param options - Optional params and functions, see {@link IMenuHeaderOptions}.
   * @returns {Menu} `this` menu itself. */
  header(text: string, options?: IMenuHeaderOptions): Menu {
    return toJs(api.grok_Menu_Header(this.dart, text, options?.onClick, options?.hasHoverEffect ?? false, options?.getDescription));
  }

  /** Adds a separator line.
   *  @returns {Menu} */
  separator(): Menu {
    return toJs(api.grok_Menu_Separator(this.dart));
  }

  /** Shows the menu.
   * @returns {Menu} `this` menu itself. */
  show(options?: IShowMenuOptions): Menu {
    return toJs(api.grok_Menu_Show(this.dart, options?.element, options?.causedBy, options?.x, options?.y, options?.nextToElement));
  }

  hide(): void {
    api.grok_Menu_Hide(this.dart);
  }

  /** Binds the menu to the specified {@link options.element} */
  bind(element: HTMLElement): Menu {
    element.oncontextmenu = (ev) => {
      ev.preventDefault();
      this.show({causedBy: ev});
    }
    return this;
  }

  get onContextMenuItemClick() {
    return __obs('d4-menu-item-click', this.dart);
  }

  get onClose(): Observable<EventData> { return api.grok_Menu_OnClose(this.dart); }

  toString(): string {
    return api.grok_MenuItem_ToString(this.dart);
  }
}


/** Balloon-style visual notifications. */
export class Balloon {

  /** Shows information message (green background) */
   info(s: string | HTMLElement): void {
    api.grok_Balloon(s, 'info', toDart({}));
  }

  /** Shows information message (red background) */
  error(s: string | HTMLElement): void {
    api.grok_Balloon(s, 'error', toDart({}));
  }

  warning(s: string | HTMLElement): void {
    api.grok_Balloon(s, 'warning', toDart({}));
  }

  /** Closes all balloons currently shown */
  static closeAll(): void {
    api.grok_Balloon_CloseAll();
  }
}
