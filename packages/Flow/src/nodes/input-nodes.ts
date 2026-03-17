import {LGraphNode, LiteGraph} from 'litegraph.js';
import {getSlotColor} from '../types/type-map';

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

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Column List Input ---

class ColumnListInputNode extends LGraphNode {
  static title = 'Column List Input';
  static desc = 'External column list input parameter';
  dgNodeType = 'input';
  dgOutputType = 'column_list';

  constructor() {
    super('Column List Input');
    this.properties = {paramName: 'cols', defaultValue: '', description: '', typeFilter: '', semTypeFilter: ''};

    const slot = this.addOutput('columns', 'column_list');
    slot.color_on = getSlotColor('column_list');
    slot.color_off = getSlotColor('column_list');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
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
    this.properties = {
      paramName: 'text', defaultValue: '', description: '',
      nullable: false, caption: '', choices: '', semType: '',
    };

    const slot = this.addOutput('value', 'string');
    slot.color_on = getSlotColor('string');
    slot.color_off = getSlotColor('string');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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
    this.properties = {
      paramName: 'value', defaultValue: 0, description: '',
      nullable: false, caption: '', min: '', max: '', showSlider: false,
    };

    const slot = this.addOutput('value', 'double');
    slot.color_on = getSlotColor('double');
    slot.color_off = getSlotColor('double');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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
    this.properties = {
      paramName: 'n', defaultValue: 0, description: '',
      nullable: false, caption: '', min: '', max: '', showSlider: false,
    };

    const slot = this.addOutput('value', 'int');
    slot.color_on = getSlotColor('int');
    slot.color_off = getSlotColor('int');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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
    this.properties = {
      paramName: 'flag', defaultValue: false, description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'bool');
    slot.color_on = getSlotColor('bool');
    slot.color_off = getSlotColor('bool');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- DateTime Input ---

class DateTimeInputNode extends LGraphNode {
  static title = 'DateTime Input';
  static desc = 'External datetime input parameter';
  dgNodeType = 'input';
  dgOutputType = 'datetime';

  constructor() {
    super('DateTime Input');
    this.properties = {
      paramName: 'date', defaultValue: '', description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'datetime');
    slot.color_on = getSlotColor('datetime');
    slot.color_off = getSlotColor('datetime');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- File Input ---

class FileInputNode extends LGraphNode {
  static title = 'File Input';
  static desc = 'External file input parameter';
  dgNodeType = 'input';
  dgOutputType = 'file';

  constructor() {
    super('File Input');
    this.properties = {
      paramName: 'file', defaultValue: '', description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'file');
    slot.color_on = getSlotColor('file');
    slot.color_off = getSlotColor('file');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Map Input ---

class MapInputNode extends LGraphNode {
  static title = 'Map Input';
  static desc = 'External map (key-value) input parameter';
  dgNodeType = 'input';
  dgOutputType = 'map';

  constructor() {
    super('Map Input');
    this.properties = {
      paramName: 'params', defaultValue: '', description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'map');
    slot.color_on = getSlotColor('map');
    slot.color_off = getSlotColor('map');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Dynamic Input ---

class DynamicInputNode extends LGraphNode {
  static title = 'Dynamic Input';
  static desc = 'External dynamically-typed input parameter';
  dgNodeType = 'input';
  dgOutputType = 'dynamic';

  constructor() {
    super('Dynamic Input');
    this.properties = {
      paramName: 'value', defaultValue: '', description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'dynamic');
    slot.color_on = getSlotColor('dynamic');
    slot.color_off = getSlotColor('dynamic');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
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

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Blob Input ---

class BlobInputNode extends LGraphNode {
  static title = 'Blob Input';
  static desc = 'External binary blob input parameter';
  dgNodeType = 'input';
  dgOutputType = 'blob';

  constructor() {
    super('Blob Input');
    this.properties = {
      paramName: 'data', defaultValue: '', description: '',
      nullable: false, caption: '',
    };

    const slot = this.addOutput('value', 'byte_array');
    slot.color_on = getSlotColor('byte_array');
    slot.color_off = getSlotColor('byte_array');

    this.color = '#66BB6A';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

/** Register all input node types */
export function registerInputNodes(): void {
  LiteGraph.registerNodeType('Inputs/Table Input', TableInputNode);
  LiteGraph.registerNodeType('Inputs/Column Input', ColumnInputNode);
  LiteGraph.registerNodeType('Inputs/Column List Input', ColumnListInputNode);
  LiteGraph.registerNodeType('Inputs/String Input', StringInputNode);
  LiteGraph.registerNodeType('Inputs/Number Input', NumberInputNode);
  LiteGraph.registerNodeType('Inputs/Int Input', IntInputNode);
  LiteGraph.registerNodeType('Inputs/Boolean Input', BooleanInputNode);
  LiteGraph.registerNodeType('Inputs/DateTime Input', DateTimeInputNode);
  LiteGraph.registerNodeType('Inputs/File Input', FileInputNode);
  LiteGraph.registerNodeType('Inputs/Map Input', MapInputNode);
  LiteGraph.registerNodeType('Inputs/Dynamic Input', DynamicInputNode);
  LiteGraph.registerNodeType('Inputs/String List Input', StringListInputNode);
  LiteGraph.registerNodeType('Inputs/Blob Input', BlobInputNode);
}
