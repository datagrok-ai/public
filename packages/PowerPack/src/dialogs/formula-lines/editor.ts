import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AnnotationRegionForm} from './annotation-region-form';
import {FormulaLineForm} from './formula-line-form';
import {ItemList} from './item-list';
import {FormContext} from './item-form';
import {EditorItem} from './types';

/**
 * Editor Helper for Formula Lines (form with corresponding inputs).
 */
export class Editor {
  /** True while a change of the item moves the axes of the preview, which must not be synced back to the item. */
  public applyingChange: boolean = false;

  private form: HTMLElement = ui.form([]);
  private formulaLineForm?: FormulaLineForm;

  public get root(): HTMLElement {return this.form;}

  /** The "Column" input of the constant line being edited. */
  public get columnInput(): DG.InputBase<DG.Column | null> | undefined {return this.formulaLineForm?.columnInput;}

  constructor(
    public items: ItemList,
    private dataFrame: DG.DataFrame,
    private onItemChangedAction: (rowIdx: number) => boolean,
    private onFormulaValidation: (item: EditorItem, isValid: boolean) => void,
    private singleAxisColumnName: string | undefined,
    private defaultItemZIndex: number,
  ) {}

  /** Creates and fills editor for the item of the [rowIdx] row */
  public update(rowIdx: number) {
    const item = this.items.at(rowIdx);
    this.formulaLineForm = undefined;

    let newForm: HTMLElement;
    if (!item)
      newForm = this.emptyForm();
    else if (this.items.isFormulaLine(rowIdx)) {
      this.formulaLineForm = new FormulaLineForm(item, this.formContext(rowIdx, item));
      newForm = this.formulaLineForm.root;
    } else
      newForm = new AnnotationRegionForm(item, this.formContext(rowIdx, item)).root;

    this.form.replaceWith(newForm);
    this.form = newForm;
  }

  /** Runs [change] with {@link applyingChange} set. */
  public applyChange(change: () => void): void {
    this.applyingChange = true;
    try {
      change();
    } finally {
      this.applyingChange = false;
    }
  }

  private formContext(rowIdx: number, item: EditorItem): FormContext {
    return {
      dataFrame: this.dataFrame,
      singleAxisColumnName: this.singleAxisColumnName,
      defaultItemZIndex: this.defaultItemZIndex,
      onItemChanged: () => this.onItemChangedAction(rowIdx),
      onFormulaValidation: (isValid) => this.onFormulaValidation(item, isValid),
      applyChange: (change) => this.applyChange(change),
    };
  }

  private emptyForm(): HTMLElement {
    return ui.div(['No formula line or annotation region selected, add one to edit.'], {classes: 'ui-form', style: {
      padding: '12px',
      textAlign: 'center',
      whiteSpace: 'normal',
      wordBreak: 'normal',
    }});
  }
}
