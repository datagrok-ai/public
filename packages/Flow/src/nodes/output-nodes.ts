import {LGraphNode, LiteGraph} from 'litegraph.js';
import {getSlotColor} from '../types/type-map';

// --- Table Output ---

class TableOutputNode extends LGraphNode {
  static title = 'Table Output';
  static desc = 'Marks a result dataframe for the script output';
  dgNodeType = 'output';
  dgOutputType = 'dataframe';

  constructor() {
    super('Table Output');
    this.properties = {paramName: 'result', description: ''};

    const slot = this.addInput('table', 'dataframe');
    slot.color_on = getSlotColor('dataframe');
    slot.color_off = getSlotColor('dataframe');

    this.addWidget('text', 'Param Name', 'result', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});

    this.color = '#EF5350';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Value Output ---

class ValueOutputNode extends LGraphNode {
  static title = 'Value Output';
  static desc = 'Marks a result value for the script output';
  dgNodeType = 'output';
  dgOutputType = 'dynamic';

  constructor() {
    super('Value Output');
    this.properties = {paramName: 'result', outputType: 'double', description: ''};

    const slot = this.addInput('value', 'dynamic');
    slot.color_on = getSlotColor('dynamic');
    slot.color_off = getSlotColor('dynamic');

    this.addWidget('text', 'Param Name', 'result', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.addWidget('combo', 'Type', 'double', (v: any) => {
      this.properties['outputType'] = v;
    }, {values: ['string', 'int', 'double', 'bool', 'dataframe', 'column', 'object'], property: 'outputType'});

    this.color = '#EF5350';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

/** Register all output node types */
export function registerOutputNodes(): void {
  LiteGraph.registerNodeType('Outputs/Table Output', TableOutputNode);
  LiteGraph.registerNodeType('Outputs/Value Output', ValueOutputNode);
}
