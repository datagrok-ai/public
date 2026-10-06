// Matches VALIDATOR_DEBOUNCE_TIME in reactive-tree-driver, so disabled styling
// settles in step with RTD validation/lock cycling.
export const DISABLED_STYLE_DEBOUNCE_TIME = 250;

export interface RibbonItemBase {
  onClick: () => void;
  tooltip?: string;
  disabled?: boolean;
  disabledReason?: string | (() => string | null);
  /** 'none' disables without any visual change */
  disabledStyle?: 'default' | 'none';
  /** Debounce (ms) for the disabled style toggle, default DISABLED_STYLE_DEBOUNCE_TIME */
  disabledStyleDebounce?: number;
  /** How the disabled reason is surfaced: on hover, as a shell warning on click, or both */
  disabledReasonMode?: 'tooltip' | 'popup' | 'both' | 'none';
}

export interface RibbonPanelItem extends RibbonItemBase {
  icon: string | {path: string, width?: number, height?: number};
  active?: boolean;
}

export interface RibbonMenuItem extends RibbonItemBase {
  text: string;
  icon?: string;
  check?: boolean;
}

/** The imperative surface the ribbon service writes to; DG-backed in production, fake in tests */
export interface RibbonHost {
  getPanels(): HTMLElement[][];
  setPanels(panels: HTMLElement[][]): void;
  readonly closing: boolean;
  rebuildGroup(name: string, rank: number, els: HTMLElement[],
    onClick: (el: HTMLElement) => void,
    opts: {isValid: (el: HTMLElement) => string | null, isChecked?: (el: HTMLElement) => boolean}): void;
  removeGroup(name: string): void;
  warn(message: string): void;
}

export type DebounceFactory = <A extends unknown[]>(
  fn: (...args: A) => void, wait: () => number) => (...args: A) => void;
