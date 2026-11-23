import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getTitleFromFormula, isAxisToken, wrapCol} from './formulas';
import {FormContext, ItemForm, savedColor} from './item-form';
import {ITEM_CAPTION, ITEM_TYPE} from './types';

/** Inputs of a formula line or a band. */
export class FormulaLineForm extends ItemForm<DG.FormulaLine> {
  public readonly root: HTMLElement;

  /** The "Column" input of a constant line. */
  public columnInput?: DG.InputBase<DG.Column | null>;

  // Kept for the other inputs, because the title may depend on the formula
  private ibTitle?: DG.InputBase;

  constructor(item: DG.FormulaLine, ctx: FormContext) {
    super(item, ctx);

    let caption = ITEM_CAPTION.BAND;
    let [itemY, expression] = ['', ''];
    if (item.type !== ITEM_TYPE.BAND) {
      // getMeta parses the formula and throws on invalid input — keep the form alive so
      // the user can fix an in-progress edit instead of crashing on row switch.
      try {
        const itemMeta = DG.FormulaLinesHelper.getMeta(item);
        [itemY, expression] = [itemMeta.funcName!, itemMeta.expression!];
        caption = itemMeta.argName ? ITEM_CAPTION.LINE : ITEM_CAPTION.CONST_LINE;
      } catch {
        caption = ITEM_CAPTION.LINE;
      }
    }

    const mainPane = this.createPane();
    const formatPane = this.createPane();
    const tooltipPane = this.createPane();

    // An axis token has no column to put in the Column picker, and picking one there would
    // rewrite `Y = 150` into `${someColumn} = 150`; edit those as text instead.
    mainPane.append(caption === ITEM_CAPTION.CONST_LINE && !isAxisToken(itemY) ?
      this.inputConstant(itemY, expression) :
      this.inputFormula());
    if (caption === ITEM_CAPTION.BAND)
      mainPane.append(this.inputColumn2());

    formatPane.append(this.inputColor());
    if (caption !== ITEM_CAPTION.BAND)
      formatPane.append(this.inputStyle());
    formatPane.append(this.inputRange());
    // Lines default to above markers; bands fall back to the per-host default.
    formatPane.append(this.inputArrange(item.type === ITEM_TYPE.BAND ? ctx.defaultItemZIndex : 100));

    /** "Show on plot" and "Show on tooltip" share a row */
    tooltipPane.append(this.inputTitle());
    tooltipPane.append(this.inlineInputs(
      this.inputBool('Show on plot', 'showOnPlot'), this.inputBool('Show on tooltip', 'showOnTooltip')));
    tooltipPane.append(this.inputDescription());

    this.root = this.createAccordion({[caption]: mainPane, 'Format': formatPane, 'Tooltip': tooltipPane});
  }

  private setTitleIfEmpty(oldFormula: string, newFormula: string) {
    const placeholder = (this.ibTitle!.input as HTMLInputElement).placeholder;
    if ([oldFormula, getTitleFromFormula(oldFormula)].includes(placeholder) || newFormula === this.ibTitle!.value)
      this.ibTitle!.value = newFormula;
  }

  private inputFormula(): HTMLElement {
    const item = this.item;
    const isBand = item.type === ITEM_TYPE.BAND;

    const ibFormula = ui.input.textArea('', {value: item.formula ?? '',
      onValueChanged: (value) => {
        const oldFormula = item.formula!;
        item.formula = value;
        const resultOk = this.ctx.onItemChanged();
        this.setFormulaValidationResult(resultOk, value, ibFormula, isBand);
        this.setTitleIfEmpty(oldFormula, item.formula);
      }});

    const elFormula = ibFormula.input as HTMLInputElement;
    elFormula.placeholder = 'Formula';
    elFormula.setAttribute('style', 'height: 60px;');

    ui.tools.initFormulaAccelerators(ibFormula, this.ctx.dataFrame);

    // Flag pre-existing invalid bands without waiting for the user to edit.
    if (isBand)
      this.setFormulaValidationResult(true, item.formula ?? '', ibFormula, true);

    return ibFormula.root;
  }

  private inputConstant(colName: string, value: string): HTMLElement {
    const item = this.item;

    const setFormula = (column: DG.Column, value: string) => {
      const oldFormula = item.formula!;
      item.formula = `${wrapCol(column.name)} = ${value}`;
      this.ctx.onItemChanged();
      this.setTitleIfEmpty(oldFormula, item.formula);
    };

    const ibColumn = this.inputColumn('Column', colName, (column) => setFormula(column, ibValue.value));
    this.columnInput = ibColumn;

    const ibValue = ui.input.string('Value', {value: value,
      onValueChanged: (value) => setFormula(ibColumn.value as DG.Column, value)});
    ibValue.nullable = false;

    return ui.div([ibColumn.root, ibValue.root], 'ui-form');
  }

  private inputColumn2(): HTMLElement {
    return ui.divH([this.inputColumn('Adjacent column', this.item.column2, (column) => this.ctx.applyChange(() => {
      this.item.column2 = column.name;
      this.ctx.onItemChanged();
    })).root]);
  }

  // A single-axis viewer offers only its value column
  private inputColumn(caption: string, colName: string | undefined,
    onValueChanged: (column: DG.Column) => void): DG.InputBase<DG.Column | null> {
    const axisName = this.ctx.singleAxisColumnName;
    //@ts-ignore
    return ui.input.column(caption, {
      nullable: false,
      table: this.ctx.dataFrame,
      filter: axisName ? (c: DG.Column) => c.name === axisName : undefined,
      value: colName ? this.ctx.dataFrame.col(colName) ?? undefined : undefined,
      onValueChanged: onValueChanged,
    });
  }

  // The alpha of the color is the opacity. A color saved before the alpha is read
  // the way the viewers draw it, and written back in the new format on the first edit.
  private inputColor(): HTMLElement {
    const item = this.item;
    const defaultColor = '#000000';
    const c = savedColor(DG.FormulaLinesHelper.setDefaults(item).color ?? defaultColor);

    const ibColor = ui.input.color('Color', {
      value: DG.Color.toHtml(c), useAlphaChannel: true,
      onValueChanged: (value) => {
        delete item.opacity;
        // "#rrggbb" is the format saved before the alpha, which draws a band at 30%
        item.color = value.length === 7 ? `${value}ff` : value;
        this.ctx.onItemChanged();
      }});
    (ibColor.input as HTMLInputElement).placeholder = defaultColor;
    return ibColor.root;
  }

  private inputStyle(): HTMLElement {
    const ibStyle = ui.input.choice('Style', {value: this.item.style ?? 'solid',
      items: ['solid', 'dotted', 'dashed', 'longdash', 'dotdash'], onValueChanged: this.setField('style')});

    const ibWidth = ui.input.int('', {value: this.item.width ?? 1, onValueChanged: this.setField('width')});
    ibWidth.addPostfix('px');

    const elWidth = ibWidth.input as HTMLInputElement;
    elWidth.placeholder = '1';
    elWidth.setAttribute('style', 'width: 61px; padding-right: 24px;');

    return ui.divH([ibStyle.root, ibWidth.root]);
  }

  private inputRange(): HTMLElement {
    return ui.divH([this.inputRangeBound('min', 'Range'), this.inputRangeBound('max', '')]);
  }

  private inputRangeBound(bound: 'min' | 'max', caption: string): HTMLElement {
    const item = this.item;

    const ibBound = ui.input.string(caption, {value: `${item[bound] ?? ''}`,
      onValueChanged: (value) => {
        item[bound] = value.length === 0 ? undefined : Number(value);
        this.ctx.onItemChanged();
      }});

    const elBound = ibBound.input as HTMLInputElement;
    elBound.placeholder = bound;
    elBound.setAttribute('style', 'width: 98px;');

    return ibBound.root;
  }

  private inputTitle(): HTMLElement {
    const item = this.item;
    if (!item.title || item.title === item.formula)
      item.title = getTitleFromFormula(item.formula!);

    /** A title equal to the formula is shown as the placeholder, so that it keeps following the formula. */
    const formTitleValue = (value: string): string => {
      const potentialValue = getTitleFromFormula(item.formula!);
      if (value === '' || value === item.formula || value === potentialValue) {
        elTitle.placeholder = potentialValue;
        elTitle.value = '';
        return potentialValue;
      }
      elTitle.placeholder = '';
      return value;
    };

    this.ibTitle = ui.input.string('Title', {value: item.title ?? '',
      onValueChanged: (value) => {
        item.title = formTitleValue(value);
        this.ctx.onItemChanged();
      }});

    const elTitle = this.ibTitle.input as HTMLInputElement;
    elTitle.setAttribute('style', 'width: 204px; max-width: none;');
    formTitleValue(elTitle.value);

    return ui.divH([this.ibTitle.root]);
  }

  private inputBool(caption: string, key: 'showOnPlot' | 'showOnTooltip'): HTMLElement {
    return ui.input.bool(caption, {value: this.item[key] ?? true, onValueChanged: this.setField(key)}).root;
  }
}
