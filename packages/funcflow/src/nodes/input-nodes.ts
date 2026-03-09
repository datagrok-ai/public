import {LGraphNode, LiteGraph} from 'litegraph.js';
import {getSlotColor} from '../types/type-map';

/** Base class for all input (source) nodes - these become //input: lines in generated scripts */
class BaseInputNode extends LGraphNode {
  static category = 'Inputs';
  dgNodeType = 'input';

  constructor(title: string, outputType: string, defaultParamName: string) {
    super(title);
    this.properties = {
      paramName: defaultParamName,
      defaultValue: '',
      description: '',
    };

    const slot = this.addOutput('value', outputType);
    slot.color_on = getSlotColor(outputType);
    slot.color_off = getSlotColor(outputType);

    this.addWidget('text', 'Param Name', defaultParamName, (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
  }
}

// --- Table Input ---

class TableInputNode extends LGraphNode {
  static title = 'Table Input';
  static desc = 'External dataframe input parameter';
  dgNodeType = 'input';
  dgOutputType = 'dataframe';

  constructor() {
    super('Table Input');
    this.properties = {paramName: 'df', defaultValue: '', description: ''};

    const slot = this.addOutput('table', 'dataframe');
    slot.color_on = getSlotColor('dataframe');
    slot.color_off = getSlotColor('dataframe');

    this.addWidget('text', 'Param Name', 'df', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Column Input ---

class ColumnInputNode extends LGraphNode {
  static title = 'Column Input';
  static desc = 'External column input parameter';
  dgNodeType = 'input';
  dgOutputType = 'column';

  constructor() {
    super('Column Input');
    this.properties = {paramName: 'col', defaultValue: '', description: '', typeFilter: '', semTypeFilter: ''};

    const slot = this.addOutput('column', 'column');
    slot.color_on = getSlotColor('column');
    slot.color_off = getSlotColor('column');

    this.addWidget('text', 'Param Name', 'col', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('combo', 'Type Filter', '', (v: any) => {
      this.properties['typeFilter'] = v;
    }, {values: ['', 'numerical', 'categorical', 'string', 'int', 'double', 'bool'], property: 'typeFilter'});
    this.addWidget('text', 'SemType Filter', '', (v: any) => {
      this.properties['semTypeFilter'] = v;
    }, {property: 'semTypeFilter'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 180);
  }
}

// --- String Input ---

class StringInputNode extends LGraphNode {
  static title = 'String Input';
  static desc = 'External string input parameter';
  dgNodeType = 'input';
  dgOutputType = 'string';

  constructor() {
    super('String Input');
    this.properties = {paramName: 'text', defaultValue: '', description: ''};

    const slot = this.addOutput('value', 'string');
    slot.color_on = getSlotColor('string');
    slot.color_off = getSlotColor('string');

    this.addWidget('text', 'Param Name', 'text', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('text', 'Default', '', (v: any) => {
      this.properties['defaultValue'] = v;
    }, {property: 'defaultValue'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Number Input ---

class NumberInputNode extends LGraphNode {
  static title = 'Number Input';
  static desc = 'External number (double) input parameter';
  dgNodeType = 'input';
  dgOutputType = 'double';

  constructor() {
    super('Number Input');
    this.properties = {paramName: 'value', defaultValue: 0, description: ''};

    const slot = this.addOutput('value', 'double');
    slot.color_on = getSlotColor('double');
    slot.color_off = getSlotColor('double');

    this.addWidget('text', 'Param Name', 'value', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('number', 'Default', 0, (v: any) => {
      this.properties['defaultValue'] = v;
    }, {precision: 3, property: 'defaultValue'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Int Input ---

class IntInputNode extends LGraphNode {
  static title = 'Int Input';
  static desc = 'External integer input parameter';
  dgNodeType = 'input';
  dgOutputType = 'int';

  constructor() {
    super('Int Input');
    this.properties = {paramName: 'n', defaultValue: 0, description: ''};

    const slot = this.addOutput('value', 'int');
    slot.color_on = getSlotColor('int');
    slot.color_off = getSlotColor('int');

    this.addWidget('text', 'Param Name', 'n', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('number', 'Default', 0, (v: any) => {
      this.properties['defaultValue'] = Math.round(v);
    }, {precision: 0, property: 'defaultValue'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Boolean Input ---

class BooleanInputNode extends LGraphNode {
  static title = 'Boolean Input';
  static desc = 'External boolean input parameter';
  dgNodeType = 'input';
  dgOutputType = 'bool';

  constructor() {
    super('Boolean Input');
    this.properties = {paramName: 'flag', defaultValue: false, description: ''};

    const slot = this.addOutput('value', 'bool');
    slot.color_on = getSlotColor('bool');
    slot.color_off = getSlotColor('bool');

    this.addWidget('text', 'Param Name', 'flag', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('toggle', 'Default', false, (v: any) => {
      this.properties['defaultValue'] = v;
    }, {property: 'defaultValue'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- String List Input ---

class StringListInputNode extends LGraphNode {
  static title = 'String List Input';
  static desc = 'External string list input parameter';
  dgNodeType = 'input';
  dgOutputType = 'string_list';

  constructor() {
    super('String List Input');
    this.properties = {paramName: 'items', defaultValue: '', description: ''};

    const slot = this.addOutput('value', 'string_list');
    slot.color_on = getSlotColor('string_list');
    slot.color_off = getSlotColor('string_list');

    this.addWidget('text', 'Param Name', 'items', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('text', 'Default (comma-sep)', '', (v: any) => {
      this.properties['defaultValue'] = v;
    }, {property: 'defaultValue'});

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 200);
  }
}

/** Register all input node types */
export function registerInputNodes(): void {
  LiteGraph.registerNodeType('Inputs/Table Input', TableInputNode);
  LiteGraph.registerNodeType('Inputs/Column Input', ColumnInputNode);
  LiteGraph.registerNodeType('Inputs/String Input', StringInputNode);
  LiteGraph.registerNodeType('Inputs/Number Input', NumberInputNode);
  LiteGraph.registerNodeType('Inputs/Int Input', IntInputNode);
  LiteGraph.registerNodeType('Inputs/Boolean Input', BooleanInputNode);
  LiteGraph.registerNodeType('Inputs/String List Input', StringListInputNode);
}
