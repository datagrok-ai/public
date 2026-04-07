/* Do not change these import lines to match external modules in webpack configuration */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';

import {RadioButtonFilter} from './filters/radio-button-filter';
import {MultiValueFilter} from './filters/multi-value-filter';
import {TableSummaryViewer, TimeWidget} from './widgets';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func({
    name: 'Single Choice',
    description: 'A filter that lets you select exactly one category',
    outputs: [{type: 'filter', name: 'result'}],
    meta: {role: 'filter'},
  })
  static radioButtonFilter() {
    return new RadioButtonFilter();
  }

  @grok.decorators.func({
    name: 'Multi Choice',
    description: 'A filter that works with columns of multi-value cells (such as lists of identifiers)',
    outputs: [{type: 'filter', name: 'result'}],
    meta: {role: 'filter'},
  })
  static multiValueFilter() {
    return new MultiValueFilter();
  }

  @grok.decorators.func({
    name: 'TimeWidget',
    description: 'Shows current time',
  })
  static timeWidget() : DG.Widget {
    return new TimeWidget();
  }


  @grok.decorators.func({
    name: 'TableSummary',
    meta: {showInGallery: 'false', role: 'viewer'},
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static tableSummary() {
    return new TableSummaryViewer();
  }


  @grok.decorators.func()
  static inputDemo() {
    const medication = {name: 'Aspirin', quantity: '20 mg'};
    const fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});

    ui.dialog()
      .add(ui.input.form(medication, [DG.Property.js('name', DG.TYPE.STRING), fooProp]))
      .show();
  }


  @grok.decorators.func({
    meta: {
      propertyType: 'string',
      semType: 'foo',
      role: 'valueEditor',
    },
    outputs: [{type: 'object', name: 'result'}],
  })
  static fooInput(): DG.InputBase {
    return new FooInput();
  }
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
