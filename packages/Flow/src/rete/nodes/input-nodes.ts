/** Input nodes — emit `//input:` annotation lines in the generated script.
 *
 * Each input node has exactly one output socket carrying the declared DG type.
 * Per-type qualifiers (nullable, min/max, semType, choices, …) live in
 * `node.properties` and are rendered into the annotation line by the emitter.
 *
 * Every input node carries an inline value editor on its body (the `value`
 * control → `InputValueControl`), mirrored in the side panel: a configured
 * value is fed straight into the prepared run, so neither Run nor autorun
 * needs the parameter dialog (see `utils/input-values.ts`). */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';
import {InputValueControl} from './input-value-control';

const COLOR_INPUT = categoricalColor(CAT.green);

abstract class InputBase extends FlowNode {
  constructor(label: string, paramName: string, dgType: string, slotName = 'value', extraProps: Record<string, any> = {}) {
    super(label);
    this.dgNodeType = 'input';
    this.dgOutputType = dgType;
    this.properties = {paramName, defaultValue: '', ...extraProps};
    (this as unknown as {color: string}).color = COLOR_INPUT;
    this.addOutput(slotName, new ClassicPreset.Output(getSocket(dgType), slotName));
    // Inline value editor — built lazily on first render (reads the FINAL
    // dgOutputType, so subclass overrides like Blob's apply).
    this.addControl('value', new InputValueControl(this));
  }
}

export class TableInputNode extends InputBase {
  constructor() { super('Table Input', 'df', 'dataframe', 'table'); }
}

export class ColumnInputNode extends InputBase {
  constructor() {
    super('Column Input', 'col', 'column', 'column', {typeFilter: '', semTypeFilter: ''});
  }
}

export class ColumnListInputNode extends InputBase {
  constructor() {
    super('Column List Input', 'cols', 'column_list', 'columns', {typeFilter: '', semTypeFilter: ''});
  }
}

export class StringInputNode extends InputBase {
  constructor() {
    super('String Input', 'text', 'string', 'value',
      {nullable: false, caption: '', choices: '', semType: ''});
  }
}

export class NumberInputNode extends InputBase {
  constructor() {
    super('Number Input', 'value', 'double', 'value',
      {nullable: false, caption: '', min: '', max: '', showSlider: false});
    this.properties['defaultValue'] = 0;
  }
}

export class IntInputNode extends InputBase {
  constructor() {
    super('Int Input', 'n', 'int', 'value',
      {nullable: false, caption: '', min: '', max: '', showSlider: false});
    this.properties['defaultValue'] = 0;
  }
}

export class BooleanInputNode extends InputBase {
  constructor() {
    super('Boolean Input', 'flag', 'bool', 'value', {nullable: false, caption: ''});
    this.properties['defaultValue'] = false;
  }
}

export class DateTimeInputNode extends InputBase {
  constructor() {
    super('DateTime Input', 'date', 'datetime', 'value', {nullable: false, caption: ''});
  }
}

export class FileInputNode extends InputBase {
  constructor() {
    super('File Input', 'file', 'file', 'value', {nullable: false, caption: ''});
  }
}

export class MapInputNode extends InputBase {
  constructor() {
    super('Map Input', 'params', 'map', 'value', {nullable: false, caption: ''});
  }
}

export class DynamicInputNode extends InputBase {
  constructor() {
    super('Dynamic Input', 'value', 'dynamic', 'value', {nullable: false, caption: ''});
  }
}

export class StringListInputNode extends InputBase {
  constructor() { super('String List Input', 'items', 'string_list', 'value'); }
}

export class BlobInputNode extends InputBase {
  constructor() {
    super('Blob Input', 'data', 'byte_array', 'value', {nullable: false, caption: ''});
    // The DG annotation type is `blob`, even though the slot type maps to `byte_array`.
    this.dgOutputType = 'blob';
  }
}
