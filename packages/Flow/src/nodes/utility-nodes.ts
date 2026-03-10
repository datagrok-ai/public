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

    this.addWidget('text', 'Column Name', '', (v: any) => {
      this.properties['columnName'] = v;
    }, {property: 'columnName'});

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 180);
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

    this.addWidget('text', 'Column Names', '', (v: any) => {
      this.properties['columnNames'] = v;
    }, {property: 'columnNames'});

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 200);
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

    this.addWidget('text', 'Label', '', (v: any) => {
      this.properties['label'] = v;
    }, {property: 'label'});

    this.color = '#78909C';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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
    this.size[0] = Math.max(this.size[0], 140);
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
    this.size[0] = Math.max(this.size[0], 140);
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
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Constant String ---

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

    this.addWidget('number', 'Value', 0, (v: any) => {
      this.properties['value'] = Math.round(v);
    }, {precision: 0, property: 'value'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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

    this.addWidget('number', 'Value', 0, (v: any) => {
      this.properties['value'] = v;
    }, {precision: 4, property: 'value'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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

    this.addWidget('toggle', 'Value', false, (v: any) => {
      this.properties['value'] = v;
    }, {property: 'value'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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

    this.addWidget('text', 'Values', '', (v: any) => {
      this.properties['value'] = v;
    }, {property: 'value'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 180);
  }
}

/** Register all utility node types */
export function registerUtilityNodes(): void {
  LiteGraph.registerNodeType('Utilities/Select Column', SelectColumnNode);
  LiteGraph.registerNodeType('Utilities/Select Columns', SelectColumnsNode);
  LiteGraph.registerNodeType('Utilities/Log', LogNode);
  LiteGraph.registerNodeType('Utilities/Info', InfoNode);
  LiteGraph.registerNodeType('Utilities/Warning', WarningNode);
  LiteGraph.registerNodeType('Utilities/ToString', ToStringNode);
  LiteGraph.registerNodeType('Constants/String', ConstStringNode);
  LiteGraph.registerNodeType('Constants/Int', ConstIntNode);
  LiteGraph.registerNodeType('Constants/Double', ConstDoubleNode);
  LiteGraph.registerNodeType('Constants/Boolean', ConstBoolNode);
  LiteGraph.registerNodeType('Constants/List', ConstListNode);
}
