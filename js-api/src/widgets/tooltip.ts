/**
 * Tooltip class and the platform's single tooltip instance.
 * @module widgets/tooltip
 */

import * as rxjs from 'rxjs';

import {IDartApi} from '../api/grok_api.g';
import {IndexPredicate} from '../const';
import {DataFrame} from '../dataframe';
import {__obs} from '../events';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** Options for tooltip display. Future-proof — add fields here as tooltip features grow. */
export interface ITooltipOptions {
  /** Delay in milliseconds before the tooltip appears. Defaults to 0 (immediate). */
  delay?: number;
}

/** Represents a tooltip. */
export class Tooltip {
  private _pendingTimer: ReturnType<typeof setTimeout> | null = null;

  private _cancelPending(): void {
    if (this._pendingTimer !== null) {
      clearTimeout(this._pendingTimer);
      this._pendingTimer = null;
    }
  }

  /** Hides the tooltip. Also cancels any pending {@link showDelayed} call. */
  hide(): void {
    this._cancelPending();
    api.grok_Tooltip_Hide();
  }

  /** Associated the specified visual element with the corresponding item.
   * Example: {@link https://public.datagrok.ai/js/samples/ui/tooltips/tooltips}
  */
  bind(element: HTMLElement, tooltip?: string | null | (() => string | HTMLElement | null), tooltipPosition?: 'left' | 'right' | 'top' | 'bottom' | undefined | null): HTMLElement {
    if (tooltip != null)
      api.grok_Tooltip_SetOn(element, tooltip, tooltipPosition);
    return element;
  }

  /** Shows the tooltip at the specified position.
   *
   * Any pending delayed show is cancelled first. Passing `null`/`undefined` for `content`
   * hides the tooltip and cancels any pending show — ideal for hover handlers where a single
   * call replaces the show+hide+debounce quartet:
   *
   * ```ts
   * onMouseMove(e) {
   *   const t = getTooltipFor(e);   // string | HTMLElement | null
   *   ui.tooltip.show(t, e.x + 16, e.y + 16, {delay: 200});
   * }
   * ```
   */
  show(content: HTMLElement | string | null | undefined, x: number, y: number,
       options?: ITooltipOptions): void {
    this._cancelPending();
    if (content == null) {
      api.grok_Tooltip_Hide();
      return;
    }
    const delay = options?.delay ?? 0;
    if (delay <= 0) {
      api.grok_Tooltip_Show(content, x, y);
      return;
    }
    this._pendingTimer = setTimeout(() => {
      this._pendingTimer = null;
      api.grok_Tooltip_Show(content, x, y);
    }, delay);
  }

  showRowGroup(dataFrame: DataFrame, indexPredicate: IndexPredicate, x: number, y: number): void {
    api.grok_Tooltip_ShowRowGroup(dataFrame.dart, indexPredicate, x, y);
  }

  /** Returns a tooltip element. */
  get root(): HTMLElement {
    return api.grok_Tooltip_Get_Root();
  }

  get isVisible(): boolean { return api.grok_Tooltip_Get_IsVisible(); }

  get onTooltipRequest(): rxjs.Observable<any> { return __obs('d4-tooltip-request'); }
  get onTooltipShown(): rxjs.Observable<any> { return __obs('d4-tooltip-shown'); }
  get onTooltipClosed(): rxjs.Observable<any> { return __obs('d4-tooltip-closed'); }
}

export const tooltip = new Tooltip();
