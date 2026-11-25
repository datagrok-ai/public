import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {validateBandFormula, validateLineFormula} from './formulas';
import {EditorItem} from './types';

/** What a form needs from the editor it belongs to. */
export interface FormContext {
  dataFrame: DG.DataFrame,
  /** Single-axis viewers (box plot, histogram, bar chart) restrict the column picker
   *  to their value column. `undefined` means "no restriction" (scatter / line / density). */
  singleAxisColumnName?: string,
  /** Per-host default for the Arrange dropdown when an item's `zIndex` is unset. */
  defaultItemZIndex: number,
  /** Shows the changed item in the table and on the preview; false if the preview failed to render it. */
  onItemChanged: () => boolean,
  onFormulaValidation: (isValid: boolean) => void,
  /** Runs [change] so that the axes it moves on the preview are not synced back to the item. */
  applyChange: (change: () => void) => void,
}

/** A color as the viewers draw it: layouts saved before the alpha may hold an ARGB number, and no alpha is opaque. */
export const savedColor = (color: number | string): number =>
  DG.Color.opaqueIfNoAlpha(typeof color === 'string' ? DG.Color.fromHtml(color) : color);

/** Inputs of one item, grouped into accordion panes; the parts shared by the lines and the regions. */
export abstract class ItemForm<T extends EditorItem> {
  public abstract readonly root: HTMLElement;

  constructor(protected item: T, protected ctx: FormContext) {}

  protected createPane(): HTMLElement {
    return ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});
  }

  protected createAccordion(panes: {[name: string]: HTMLElement}): HTMLElement {
    const combinedPanels = ui.accordion();
    for (const [name, pane] of Object.entries(panes))
      combinedPanels.addPane(name, () => pane, true);
    return ui.div([combinedPanels.root]);
  }

  // Pairs ui-input-roots on a single row when there's enough width; wraps to
  // separate rows when the container is narrower than ~2× the basis.
  protected inlineInputs(...children: HTMLElement[]): HTMLElement {
    const row = ui.divH(children);
    row.style.flexWrap = 'wrap';
    row.style.gap = '4px 8px';
    for (const child of children) {
      // flex-shrink: 0 so wrap is decided by the basis, not by intrinsic min-content
      // (color picker shrinks far below 180px, while a bool input doesn't, which would
      // otherwise leave the color row side-by-side at narrow widths and the bool row stacked).
      child.style.flex = '1 0 180px';
      child.style.minWidth = '0';
    }
    return row;
  }

  // A handler that writes the value of an input to the [key] field of the item
  protected setField<K extends keyof T>(key: K): (value: T[K]) => void {
    return (value) => {
      this.item[key] = value;
      this.ctx.onItemChanged();
    };
  }

  protected inputDescription(): HTMLElement {
    const ibDescription = ui.input.textArea('Description', {value: this.item.description ?? '',
      onValueChanged: this.setField('description')});
    ibDescription.input.setAttribute('style', 'height: 40px; font-family: inherit; font-size: inherit;');
    return ibDescription.root;
  }

  // [fallbackZIndex] applies while `zIndex` is unset
  protected inputArrange(fallbackZIndex: number): HTMLElement {
    const effectiveZIndex = this.item.zIndex ?? fallbackZIndex;

    const ibArrange = ui.input.choice('Arrange', {
      value: effectiveZIndex >= 0 ? 'above markers' : 'below markers',
      items: ['above markers', 'below markers'],
      onValueChanged: (value) => {
        this.item.zIndex = value === 'above markers' ? 100 : -100;
        this.ctx.onItemChanged();
      }});

    ibArrange.input.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibArrange.root]);
  }

  // [resultOk] is false when the preview failed to render the formula
  protected setFormulaValidationResult(resultOk: boolean, value: string, ibFormula: DG.InputBase<string>,
    isBand: boolean = false): void {
    const structureError = isBand ?
      validateBandFormula(value) :
      validateLineFormula(value, this.ctx.singleAxisColumnName);
    // Expression syntax validation comes last
    const validationTooltip = structureError || (resultOk ? '' : 'Invalid formula syntax');
    resultOk = resultOk && validationTooltip === '';
    ibFormula.setTooltip(validationTooltip);
    (ibFormula.input as HTMLInputElement).classList.toggle('d4-forced-invalid', !resultOk);
    this.ctx.onFormulaValidation(resultOk);
  }
}
