/** Utility nodes (helpers, type-conversions) and constant nodes.
 *
 * Most nodes have NO inline widgets — properties are edited in the side
 * panel. The exception is `ConstStringNode`, which keeps a `text` control on
 * the node body for quick editing (mirrors the LiteGraph behavior). */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';

const COLOR_UTILITY = '#78909C';
const COLOR_INFO = '#66BB6A';
const COLOR_WARNING = '#FFA726';
const COLOR_VIEW = '#42A5F5';
const COLOR_CONST = '#66BB6A';

// ---------- helpers ----------

export class SelectColumnNode extends FlowNode {
  constructor() {
    super('Select Column');
    this.dgNodeType = 'utility';
    this.properties = {columnName: ''};
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    this.addOutput('column', new ClassicPreset.Output(getSocket('column'), 'column'));
  }
}

export class SelectColumnsNode extends FlowNode {
  constructor() {
    super('Select Columns');
    this.dgNodeType = 'utility';
    this.properties = {columnNames: ''};
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    this.addOutput('columns', new ClassicPreset.Output(getSocket('column_list'), 'columns'));
  }
}

export class AddTableViewNode extends FlowNode {
  constructor() {
    super('Add Table View');
    this.dgNodeType = 'utility';
    this.properties = {};
    (this as unknown as {color: string}).color = COLOR_VIEW;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    this.addOutput('view', new ClassicPreset.Output(getSocket('view'), 'view'));
  }
}

export class LogNode extends FlowNode {
  constructor() {
    super('Log');
    this.dgNodeType = 'utility';
    this.properties = {label: ''};
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
  }
}

export class InfoNode extends FlowNode {
  constructor() {
    super('Info');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_INFO;
    this.addInput('message', new ClassicPreset.Input(getSocket('string'), 'message'));
  }
}

export class WarningNode extends FlowNode {
  constructor() {
    super('Warning');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_WARNING;
    this.addInput('message', new ClassicPreset.Input(getSocket('string'), 'message'));
  }
}

export class ToStringNode extends FlowNode {
  constructor() {
    super('ToString');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
    this.addOutput('text', new ClassicPreset.Output(getSocket('string'), 'text'));
  }
}

export class FromJsonNode extends FlowNode {
  constructor() {
    super('FromJSON');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('json', new ClassicPreset.Input(getSocket('string'), 'json'));
    this.addOutput('value', new ClassicPreset.Output(getSocket('dynamic'), 'value'));
  }
}

export class ToJsonNode extends FlowNode {
  constructor() {
    super('ToJSON');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
    this.addOutput('json', new ClassicPreset.Output(getSocket('string'), 'json'));
  }
}

// ---------- constants ----------

/** The only node that keeps an inline widget — a text control for the value. */
export class ConstStringNode extends FlowNode {
  constructor() {
    super('String');
    this.dgNodeType = 'utility';
    this.properties = {value: ''};
    (this as unknown as {color: string}).color = COLOR_CONST;
    this.addOutput('value', new ClassicPreset.Output(getSocket('string'), 'value'));

    const control = new ClassicPreset.InputControl('text', {
      initial: '',
      change: (v) => {this.properties['value'] = v ?? '';},
    });
    this.addControl('value', control);
  }
}

export class ConstIntNode extends FlowNode {
  constructor() {
    super('Int');
    this.dgNodeType = 'utility';
    this.properties = {value: 0};
    (this as unknown as {color: string}).color = COLOR_CONST;
    this.addOutput('value', new ClassicPreset.Output(getSocket('int'), 'value'));
  }
}

export class ConstDoubleNode extends FlowNode {
  constructor() {
    super('Double');
    this.dgNodeType = 'utility';
    this.properties = {value: 0.0};
    (this as unknown as {color: string}).color = COLOR_CONST;
    this.addOutput('value', new ClassicPreset.Output(getSocket('double'), 'value'));
  }
}

export class ConstBoolNode extends FlowNode {
  constructor() {
    super('Boolean');
    this.dgNodeType = 'utility';
    this.properties = {value: false};
    (this as unknown as {color: string}).color = COLOR_CONST;
    this.addOutput('value', new ClassicPreset.Output(getSocket('bool'), 'value'));
  }
}

export class ConstListNode extends FlowNode {
  constructor() {
    super('List');
    this.dgNodeType = 'utility';
    this.properties = {value: ''};
    (this as unknown as {color: string}).color = COLOR_CONST;
    this.addOutput('value', new ClassicPreset.Output(getSocket('list'), 'value'));
  }
}
