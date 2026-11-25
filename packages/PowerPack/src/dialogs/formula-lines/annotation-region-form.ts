import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getAreaColumn, setAreaColumn} from './formulas';
import {FormContext, ItemForm, savedColor} from './item-form';
import {ITEM_TYPE} from './types';

/** Inputs of an annotation region: an area drawn on the viewer, or the space between two formula lines. */
export class AnnotationRegionForm extends ItemForm<DG.AnnotationRegion> {
  public readonly root: HTMLElement;

  constructor(item: DG.AnnotationRegion, ctx: FormContext) {
    super(item, ctx);

    const isArea = item.type === ITEM_TYPE.AREA_REGION_ANNOTATION;
    const mainPane = this.createPane();
    const descriptionPane = this.createPane();

    if (isArea) {
      mainPane.append(this.inputAreaColumn('x'));
      mainPane.append(this.inputAreaColumn('y'));
      mainPane.append(this.areaPointsInput());
    } else {
      mainPane.append(this.inputAnnotationFormula('formula1'));
      mainPane.append(this.inputAnnotationFormula('formula2'));
    }

    mainPane.append(this.inputColor());
    mainPane.append(this.inlineInputs(
      this.areaInputColor('Outline Color', 'outlineColor', DG.Color.toHtml(DG.Color.gray)),
      this.inputLineWidth(),
    ));
    mainPane.append(this.inputArrange(ctx.defaultItemZIndex));

    descriptionPane.append(this.inputAnnotationHeader());
    descriptionPane.append(this.inputDescription());
    descriptionPane.append(this.areaInputColor('Color', 'headerColor'));

    this.root = this.createAccordion({[isArea ? 'Area' : 'Formula']: mainPane, 'Description': descriptionPane});
  }

  private inputAreaColumn(axis: 'x' | 'y'): HTMLElement {
    const item = this.item as DG.AreaAnnotationRegion;

    const ibColumn = ui.input.column(axis.toUpperCase() + ' column', {
      nullable: false,
      table: this.ctx.dataFrame,
      value: this.ctx.dataFrame.col(getAreaColumn(item, axis)) ?? undefined,
      onValueChanged: (value) => {
        setAreaColumn(item, axis, value?.name);
        this.ctx.onItemChanged();
      }});

    return ui.divH([ibColumn.root]);
  }

  private areaPointsInput(): HTMLElement {
    const item = this.item as DG.AreaAnnotationRegion;
    const isPoint = (p: unknown) => Array.isArray(p) && p.length === 2 &&
      typeof p[0] === 'number' && typeof p[1] === 'number';
    const value = item.area ? JSON.stringify(item.area).replaceAll(',', ', ') : '';

    const textArea = ui.input.textArea('Points', {
      value: value.substring(1, value.length - 1),
      onValueChanged: (value) => {
        let error = 'Points should be a list of [x, y] number pairs, e.g., [1, 2], [3, 4]';
        try {
          const parsed: [number, number][] = JSON.parse(`[${value}]`);
          if (parsed.length < 3)
            error = 'Area must have at least 3 points';
          else if (parsed.every(isPoint)) {
            item.area = parsed;
            this.ctx.onItemChanged();
            error = '';
          }
        } catch {}

        textArea.setTooltip(error);
        textArea.input.classList.toggle('d4-forced-invalid', error !== '');
        this.ctx.onFormulaValidation(error === '');
      }});
    textArea.input.setAttribute('style', 'height: 75px; font-family: inherit; font-size: inherit;');

    return textArea.root;
  }

  private inputAnnotationFormula(key: 'formula1' | 'formula2'): HTMLElement {
    const item = this.item as DG.FormulaAnnotationRegion;

    const ibFormula = ui.input.textArea(key === 'formula1' ? 'Formula 1' : 'Formula 2', {
      value: item[key] ? item[key] : '',
      onValueChanged: (value) => {
        item[key] = value;
        const resultOk = this.ctx.onItemChanged();
        this.setFormulaValidationResult(resultOk, value, ibFormula);
      }});

    const elFormula = ibFormula.input as HTMLInputElement;
    elFormula.setAttribute('style', 'height: 40px; font-family: inherit; font-size: inherit;');

    ui.tools.initFormulaAccelerators(ibFormula, this.ctx.dataFrame);

    return ibFormula.root;
  }

  // The alpha of the fill is the opacity. A color saved before the alpha is read the way
  // the viewers draw it, and written back in the new format on the first edit.
  private inputColor(): HTMLElement {
    const item = this.item;
    const defaultColor = DG.Color.toHtml(DG.Color.gray);

    const ibColor = ui.input.color('Region Color', {
      value: DG.Color.toHtml(DG.AnnotationRegionsHelper.getFillColor(item)), useAlphaChannel: true,
      onValueChanged: (value) => {
        delete item.opacity;
        delete item.fillColor;
        item.fill = DG.Color.fromHtml(value);
        this.ctx.onItemChanged();
      }});
    (ibColor.input as HTMLInputElement).placeholder = defaultColor;
    return ibColor.root;
  }

  private areaInputColor(header: string, key: 'outlineColor' | 'headerColor',
    defaultColor: string = '#000000'): HTMLElement {
    const item = this.item;
    const saved: number | string | undefined = item[key];
    const ibColor = ui.input.color(header, {
      value: saved != null ? DG.Color.toHtml(savedColor(saved)) : defaultColor, useAlphaChannel: true,
      onValueChanged: (value) => {
        item[key] = DG.Color.fromHtml(value);
        this.ctx.onItemChanged();
      }});
    (ibColor.input as HTMLInputElement).placeholder = defaultColor;
    return ibColor.root;
  }

  private inputLineWidth(): HTMLElement {
    const ibWidth = ui.input.int('Width', {value: this.item.outlineWidth ?? 1,
      onValueChanged: this.setField('outlineWidth')});
    ibWidth.addPostfix('px');

    const elWidth = ibWidth.input as HTMLInputElement;
    elWidth.placeholder = '1';
    elWidth.setAttribute('style', 'padding-right: 24px;');

    return ibWidth.root;
  }

  private inputAnnotationHeader(): HTMLElement {
    const ibHeader = ui.input.string('Title', {value: this.item.header ?? '', onValueChanged: this.setField('header')});
    ibHeader.input.setAttribute('style', 'width: 204px; max-width: none;');
    return ui.divH([ibHeader.root]);
  }
}
