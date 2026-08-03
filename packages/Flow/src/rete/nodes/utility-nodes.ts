/** Utility nodes (helpers, type-conversions) and constant nodes.
 *
 * Most nodes have NO inline widgets — properties are edited in the side
 * panel. The exception is `ConstStringNode`, which keeps a `text` control on
 * the node body for quick editing (mirrors the LiteGraph behavior). */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';

const COLOR_UTILITY = categoricalColor(CAT.gray);
const COLOR_INFO = categoricalColor(CAT.green);
const COLOR_WARNING = categoricalColor(CAT.orange);
const COLOR_VIEW = categoricalColor(CAT.blue);
const COLOR_CONST = categoricalColor(CAT.green);

// ---------- helpers ----------

export class SelectColumnNode extends FlowNode {
  constructor() {
    super('Select Column');
    this.dgNodeType = 'utility';
    this.properties = {columnName: ''};
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    this.addOutput('column', new ClassicPreset.Output(getSocket('column'), 'column'));
    // Needs a table to select from and a column name to select.
    this.requiredInputs = ['table'];
    this.requiredProps = ['columnName'];
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
    // Needs a table to select from and at least one column name.
    this.requiredInputs = ['table'];
    this.requiredProps = ['columnNames'];
  }
}

/** Resolves an open table by name — `grok.shell.tableByName(name)`. Stands in
 *  for the platform `ResolveTable` function in imported creation scripts. */
export class SelectTableNode extends FlowNode {
  constructor() {
    super('Select Table');
    this.dgNodeType = 'utility';
    this.properties = {tableName: ''};
    (this as unknown as {color: string}).color = COLOR_UTILITY;
    this.addOutput('table', new ClassicPreset.Output(getSocket('dataframe'), 'table'));
    // Resolves a table by name — that name must be set.
    this.requiredProps = ['tableName'];
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
    // Needs a table to open a view for.
    this.requiredInputs = ['table'];
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
      change: (v) => {
        this.properties['value'] = v ?? '';
        // Title mirrors the value (refreshes on the next node re-render —
        // forcing one here would remount the control and steal focus).
        this.label = constLabel('String', v);
      },
    });
    this.addControl('value', control);
  }
}

/** Constant nodes title themselves after their value; empty values fall back
 *  to the type name so a freshly added node isn't titled `const: `. Long
 *  values (e.g. column-name lists) are truncated — the full value lives in
 *  `properties.value`, the title is just a hint. */
export function constLabel(kind: string, value: unknown): string {
  const s = String(value ?? '').trim();
  if (s === '') return kind;
  return `const: ${s.length > 28 ? `${s.slice(0, 27)}…` : s}`;
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
