/* Do not change these import lines to match external modules in webpack configuration */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

import {RadioButtonFilter} from './filters/radio-button-filter';
import {MultiValueFilter} from './filters/multi-value-filter';
import {TableSummaryViewer, TimeWidget} from './widgets';

export const _package = new DG.Package();

//name: Single Choice
//description: A filter that lets you select exactly one category
//tags: filter
//output: filter result
export function radioButtonFilter() {
  return new RadioButtonFilter();
}

//name: Multi Choice
//description: A filter that works with columns of multi-value cells (such as lists of identifiers)
//tags: filter
//output: filter result
export function multiValueFilter() {
  return new MultiValueFilter();
}

//name: TimeWidget
//description: Shows current time
//output: widget result
export function timeWidget() {
  return new TimeWidget();
}

//name: TableSummary
//tags: viewer
//output: viewer result
export function tableSummary() {
  return new TableSummaryViewer();
}

//name: inputDemo
export function inputDemo() {
  const medication = {name: 'Aspirin', quantity: '20 mg'};
  const fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});

  ui.dialog()
    .add(ui.input.form(medication, [DG.Property.js('name', DG.TYPE.STRING), fooProp]))
    .show();
}

//name: fooInput
//tags: valueEditor
//meta.propertyType: string
//meta.semType: foo
//output: object
export function fooInput(): DG.InputBase {
  return new FooInput();
}

export class FooInput extends DG.JsInputBase<string> {
  get inputType(): string {
    return 'foo';
  }
  get dataType(): string {
    return DG.TYPE.STRING;
  }

  valueEditor: HTMLInputElement = ui.element('input');
  unitsEditor: HTMLInputElement = ui.element('input');
  editor = ui.divH([this.valueEditor, this.unitsEditor]);

  constructor() {
    super();
    $(this.valueEditor).on('input', (e) => {this.fireInput(); this.fireChanged();});
    $(this.unitsEditor).on('input', (e) => {this.fireInput(); this.fireChanged();});
  }

  getInput(): HTMLElement {
    return this.editor;
  }

  getValue(): string {
    return this.valueEditor.value + ' ' + this.unitsEditor.value;
  }

  getStringValue(): string {return this.value;}
  setStringValue(value: string) {this.value = value;}

  setValue(value: string): void {
    if (value == this.getValue())
      return;

    const values = value.split(' ');
    this.valueEditor.value = values[0];
    this.unitsEditor.value = values.length > 1 ? values[1] : '';
    this.fireChanged();
  }
}
