import {LGraphNode, LiteGraph} from 'litegraph.js';
import {getSlotColor} from '../types/type-map';

// --- Select Column ---

class SelectColumnNode extends LGraphNode {
  static title = 'Select Column';
  static desc = 'Gets a column from a dataframe by name';
  dgNodeType = 'utility';

  constructor() {
    super('Select Column');
    this.properties = {columnName: ''};

    const inSlot = this.addInput('table', 'dataframe');
    inSlot.color_on = getSlotColor('dataframe');
    inSlot.color_off = getSlotColor('dataframe');

    const outSlot = this.addOutput('column', 'column');
    outSlot.color_on = getSlotColor('column');
    outSlot.color_off = getSlotColor('column');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Select Columns ---

class SelectColumnsNode extends LGraphNode {
  static title = 'Select Columns';
  static desc = 'Gets multiple columns from a dataframe by names (comma-separated)';
  dgNodeType = 'utility';

  constructor() {
    super('Select Columns');
    this.properties = {columnNames: ''};

    const inSlot = this.addInput('table', 'dataframe');
    inSlot.color_on = getSlotColor('dataframe');
    inSlot.color_off = getSlotColor('dataframe');

    const outSlot = this.addOutput('columns', 'column_list');
    outSlot.color_on = getSlotColor('column_list');
    outSlot.color_off = getSlotColor('column_list');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Add Table View ---

class AddTableViewNode extends LGraphNode {
  static title = 'Add Table View';
  static desc = 'Opens a dataframe in a new table view using grok.shell.addTableView()';
  dgNodeType = 'utility';

  constructor() {
    super('Add Table View');
    this.properties = {};

    const inSlot = this.addInput('table', 'dataframe');
    inSlot.color_on = getSlotColor('dataframe');
    inSlot.color_off = getSlotColor('dataframe');

    const outSlot = this.addOutput('view', 'view');
    outSlot.color_on = getSlotColor('view');
    outSlot.color_off = getSlotColor('view');

    this.color = '#42A5F5';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Log (console.log) ---

class LogNode extends LGraphNode {
  static title = 'Log';
  static desc = 'Logs a value to the console using console.log()';
  dgNodeType = 'utility';

  constructor() {
    super('Log');
    this.properties = {label: ''};

    const inSlot = this.addInput('value', 'dynamic');
    inSlot.color_on = getSlotColor('dynamic');
    inSlot.color_off = getSlotColor('dynamic');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Info (grok.shell.info) ---

class InfoNode extends LGraphNode {
  static title = 'Info';
  static desc = 'Shows an info message using grok.shell.info()';
  dgNodeType = 'utility';

  constructor() {
    super('Info');
    this.properties = {};

    const inSlot = this.addInput('message', 'string');
    inSlot.color_on = getSlotColor('string');
    inSlot.color_off = getSlotColor('string');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Warning (grok.shell.warning) ---

class WarningNode extends LGraphNode {
  static title = 'Warning';
  static desc = 'Shows a warning message using grok.shell.warning()';
  dgNodeType = 'utility';

  constructor() {
    super('Warning');
    this.properties = {};

    const inSlot = this.addInput('message', 'string');
    inSlot.color_on = getSlotColor('string');
    inSlot.color_off = getSlotColor('string');

    this.color = '#FFA726';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- ToString ---

class ToStringNode extends LGraphNode {
  static title = 'ToString';
  static desc = 'Converts any value to a string using .toString()';
  dgNodeType = 'utility';

  constructor() {
    super('ToString');
    this.properties = {};

    const inSlot = this.addInput('value', 'dynamic');
    inSlot.color_on = getSlotColor('dynamic');
    inSlot.color_off = getSlotColor('dynamic');

    const outSlot = this.addOutput('text', 'string');
    outSlot.color_on = getSlotColor('string');
    outSlot.color_off = getSlotColor('string');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Constant String (EXCEPTION: keeps widget on node) ---

class ConstStringNode extends LGraphNode {
  static title = 'String';
  static desc = 'A constant string value';
  dgNodeType = 'utility';

  constructor() {
    super('String');
    this.properties = {value: ''};

    const outSlot = this.addOutput('value', 'string');
    outSlot.color_on = getSlotColor('string');
    outSlot.color_off = getSlotColor('string');

    this.addWidget('text', 'Value', '', (v: any) => {
      this.properties['value'] = v;
    }, {property: 'value'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Constant Int ---

class ConstIntNode extends LGraphNode {
  static title = 'Int';
  static desc = 'A constant integer value';
  dgNodeType = 'utility';

  constructor() {
    super('Int');
    this.properties = {value: 0};

    const outSlot = this.addOutput('value', 'int');
    outSlot.color_on = getSlotColor('int');
    outSlot.color_off = getSlotColor('int');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Constant Double ---

class ConstDoubleNode extends LGraphNode {
  static title = 'Double';
  static desc = 'A constant double (floating point) value';
  dgNodeType = 'utility';

  constructor() {
    super('Double');
    this.properties = {value: 0.0};

    const outSlot = this.addOutput('value', 'double');
    outSlot.color_on = getSlotColor('double');
    outSlot.color_off = getSlotColor('double');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Constant Boolean ---

class ConstBoolNode extends LGraphNode {
  static title = 'Boolean';
  static desc = 'A constant boolean value';
  dgNodeType = 'utility';

  constructor() {
    super('Boolean');
    this.properties = {value: false};

    const outSlot = this.addOutput('value', 'bool');
    outSlot.color_on = getSlotColor('bool');
    outSlot.color_off = getSlotColor('bool');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- Constant List ---

class ConstListNode extends LGraphNode {
  static title = 'List';
  static desc = 'A constant list of values (comma-separated)';
  dgNodeType = 'utility';

  constructor() {
    super('List');
    this.properties = {value: ''};

    const outSlot = this.addOutput('value', 'list');
    outSlot.color_on = getSlotColor('list');
    outSlot.color_off = getSlotColor('list');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- FromJSON ---

class FromJsonNode extends LGraphNode {
  static title = 'FromJSON';
  static desc = 'Parses a JSON string into an object using JSON.parse()';
  dgNodeType = 'utility';

  constructor() {
    super('FromJSON');
    this.properties = {};

    const inSlot = this.addInput('json', 'string');
    inSlot.color_on = getSlotColor('string');
    inSlot.color_off = getSlotColor('string');

    const outSlot = this.addOutput('value', 'dynamic');
    outSlot.color_on = getSlotColor('dynamic');
    outSlot.color_off = getSlotColor('dynamic');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

// --- ToJSON ---

class ToJsonNode extends LGraphNode {
  static title = 'ToJSON';
  static desc = 'Converts a value to a JSON string using JSON.stringify()';
  dgNodeType = 'utility';

  constructor() {
    super('ToJSON');
    this.properties = {};

    const inSlot = this.addInput('value', 'dynamic');
    inSlot.color_on = getSlotColor('dynamic');
    inSlot.color_off = getSlotColor('dynamic');

    const outSlot = this.addOutput('json', 'string');
    outSlot.color_on = getSlotColor('string');
    outSlot.color_off = getSlotColor('string');

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 100);
  }
}

/** Register all utility node types */
export function registerUtilityNodes(): void {
  LiteGraph.registerNodeType('Utilities/Select Column', SelectColumnNode);
  LiteGraph.registerNodeType('Utilities/Select Columns', SelectColumnsNode);
  LiteGraph.registerNodeType('Utilities/Add Table View', AddTableViewNode);
  LiteGraph.registerNodeType('Utilities/Log', LogNode);
  LiteGraph.registerNodeType('Utilities/Info', InfoNode);
  LiteGraph.registerNodeType('Utilities/Warning', WarningNode);
  LiteGraph.registerNodeType('Utilities/ToString', ToStringNode);
  LiteGraph.registerNodeType('Constants/String', ConstStringNode);
  LiteGraph.registerNodeType('Constants/Int', ConstIntNode);
  LiteGraph.registerNodeType('Constants/Double', ConstDoubleNode);
  LiteGraph.registerNodeType('Constants/Boolean', ConstBoolNode);
  LiteGraph.registerNodeType('Constants/List', ConstListNode);
  LiteGraph.registerNodeType('Utilities/FromJSON', FromJsonNode);
  LiteGraph.registerNodeType('Utilities/ToJSON', ToJsonNode);
}
