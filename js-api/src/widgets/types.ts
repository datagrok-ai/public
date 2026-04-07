/**
 * Type definitions for widgets module.
 * @module widgets/types
 */

import type {ColumnGrid} from "../grid";
import type {Column} from "../dataframe";
import type {typeaheadConfig, Dictionary} from 'typeahead-standalone/dist/types';

export type RangeSliderStyle = 'barbell' | 'lines' | 'thin_barbell';

export type SliderOptions = {
  style?: RangeSliderStyle
}

export type ICodeEditorOptions = {
  root?: HTMLDivElement;
}

export type TypeAheadConfig = Omit<typeaheadConfig<Dictionary>, 'input' | 'className'>;

export type MarkdownConfig = {
  value?: string;
};

export type CodeConfig = {
  script?: string;
  mode?: string;
  placeholder?: string;
};

/** See {@link Menu.items} */
export interface IMenuItemsOptions<T = any> {

  /** Whether a check box appears before the item */
  isChecked?: (item: T) => boolean;

  /** If result is not null, the item is grayed out and the result is shown in the tooltip */
  isValid?: (item: T) => string | null;

  /** Text to be shown on the menu item */
  toString?: (item: T) => string;

  /** Tooltip */
  getTooltip?: (item: T) => string;

  /** Gets invoked when the mouse enters the item */
  onMouseEnter?: (item: T) => void;

  /** Identifies a group of items where only one can be checked at a time. */
  radioGroup?: string;
}

/** See {@link Menu.colorPalette} */
export interface IMenuColorPaletteOptions {

  /** Returns value put into the selector as initial or to reset. */
  getInitialValue?: () => number[];

  /** Called when color is selected by click.
   * @param {number[]} list - A sequence of selected color palette.
  */
  onSelect?: (list: number[]) => void;

  /** Called when color is hovered or reset.
   * @param {number[]} list - A sequence of selected color palette.
  */
  onPreview?: (list: number[]) => void;

  /** Either the current item is a separate subgroup by defined `string` name or inside main menu. */
  asGroup?: string | null;

  /** Whether the current item is visible. */
  visible?: boolean | null;

  /** Either the palette is displayed as a gradient (if `false`) or as a separate colors sequence (if `true`). */
  categorical?: boolean;

  /** Whether hover effect preview is allowed. */
  allowPreview?: boolean;

  /** Delay when color value reset to default after leaving hovered item. */
  resetColorMs?: number;

  /** Whether to close the menu after color is selected. */
  closeOnClick?: boolean;
}

/** See {@link Menu.fontEditor} */
export interface IMenuFontEditorOptions {
  /** Minimum font size. */
  fontSizeMin: number;

  /** Maximum font size. */
  fontSizeMax: number;

  /** Step of font size. */
  fontSizeStep: number;

  /** List of available font families. */
  fontFamilies: Iterable<String>;

  /** A name of the group in the menu. */
  asGroup: string;

  /** Called when font is changed */
  onChange: (value: string) => void;
}

/** See {@link IMenuSingleColumnSelectorOptions} and {@link IMenuMultiColumnSelectorOptions} */
export interface IMenuColumnSelectorOptions<T> {

  /** Value put into the selector as initial or to reset. */
  initialValue?: T;

  /** Called when selector value is changed */
  onChange?: (...args: any[]) => void;

  /** Either the current item is a separate subgroup by defined `string` name or inside main menu. */
  asGroup?: string | null,

  /** Whether the current item is visible. */
  visible?: boolean;

  /** Whether the current item can be changed or its visibility can be toggled. */
  editable?: boolean;

  /** Filters set selector columns to be displayed. */
  columnFilter?: (c: Column) => boolean;
}

/** See {@link Menu.singleColumnSelector} */
export interface IMenuSingleColumnSelectorOptions extends IMenuColumnSelectorOptions<string> {

  /** Called when selector value is changed.
   * @param {ColumnGrid} g - selector column grid instance.
   * @param {Column} c - A column to be selected.
   * @param {boolean} currentRowChanged - Either the current row is changed by click or just selected on hover.
  */
  onChange?: (grid: ColumnGrid, column: Column, currentRowChanged: boolean) => void;

  /** Whether selector contains empty value to indicate none selected. */
  nullable?: boolean;

  /** Whether the current item is closed after click. */
  closeOnClick?: boolean;

  /** Whether a hovered value is selected. */
  changeOnHover?: boolean;
}

/** See {@link Menu.multiColumnSelector} */
export interface IMenuMultiColumnSelectorOptions extends IMenuColumnSelectorOptions<string[]> {

  /** Called when selector value is changed.
   * @param {ColumnGrid} g - selector column grid instance.
  */
  onChange?: (grid: ColumnGrid) => void;
}

/** See {@link Menu.header} */
export interface IMenuHeaderOptions {

  /** Called when header is clicked. */
  onClick?: (item: any) => void;

  /** Whether hover effect is applied. */
  hasHoverEffect?: boolean;

  /** Tooltip to be shown on the menu item. */
  getDescription?: () => string;
}

export interface IMenuItemOptions {
  /** Identifies a group of items where only one can be checked at a time. */
  radioGroup?: string;

  /** Position in the menu */
  order?: number;

  /** Shortcut to be shown on the item. NOTE: it does not handle the keypress, just shows the shortcut*/
  shortcut?: string;

  /** Whether the menu is visible; if false, the menu is not added. Might be handy in for-loops and fluent API. */
  visible?: boolean;

  /** A function that gets called each time an item is shown.
   * Should return null if the item is enabled, otherwise the reason why it's disabled.
   * The reason for being disabled is shown in a tooltip. */
  isEnabled?: () => (string | null);

  /** For items preceded by checkboxes, indicates if the item is checked. */
  check?: boolean;

  /** Tooltip to be shown on the menu item */
  description?: string;

  /** Gets invoked when the mouse enters the item */
  onMouseEnter?: () => void;

  /** Gets invoked when the mouse leaves the item */
  onMouseLeave?: () => void;
}


export interface IShowMenuOptions {
  element?: HTMLElement,
  causedBy?: MouseEvent,
  x?: number,
  y?: number,
  nextToElement?: boolean
}

export type fileShares = 'S3';
