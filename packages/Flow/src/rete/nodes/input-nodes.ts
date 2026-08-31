/** Input nodes — emit `//input:` lines. Qualifiers live in `node.properties`; each node
 *  carries an inline value editor mirrored in the side panel (see `utils/input-values.ts`). */

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
    // Built lazily on first render so subclass dgOutputType overrides (Blob) apply.
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

/** A String Input pre-tagged `semType: Molecule` — the panel's value editor becomes
 *  Chem's compact molecule input, while the NODE body embeds a real inplace sketcher
 *  (`buildInlineSketcherEditor` in utils/input-values.ts); emits an ordinary string input line. */
export class MoleculeInputNode extends InputBase {
  constructor() {
    super('Sketcher Input', 'molecule', 'string', 'molecule',
      {nullable: false, caption: '', choices: '', semType: 'Molecule'});
  }
}

/** The macromolecule counterpart of {@link MoleculeInputNode} — `semType: Macromolecule`
 *  routes the value editor to Helm's. */
export class HelmInputNode extends InputBase {
  constructor() {
    super('Helm Input', 'sequence', 'string', 'sequence',
      {nullable: false, caption: '', choices: '', semType: 'Macromolecule'});
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
