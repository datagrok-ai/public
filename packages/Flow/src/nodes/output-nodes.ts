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

  private typeWidget: any;

  constructor() {
    super('Value Output');
    this.properties = {paramName: 'result', outputType: 'double', description: ''};

    const slot = this.addInput('value', 'dynamic');
    slot.color_on = getSlotColor('dynamic');
    slot.color_off = getSlotColor('dynamic');

    this.addWidget('text', 'Param Name', 'result', (v: any) => {
      this.properties['paramName'] = v;
    }, {property: 'paramName'});
    this.typeWidget = this.addWidget('combo', 'Type', 'double', (v: any) => {
      this.properties['outputType'] = v;
    }, {values: [
      'string', 'int', 'double', 'bool',
      'dataframe', 'column', 'column_list',
      'object', 'dynamic', 'list',
      'view', 'viewer', 'widget',
      'graphics', 'grid_cell_renderer', 'filter',
      'map', 'datetime', 'blob', 'funccall',
    ], property: 'outputType'});

    this.color = '#EF5350';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }

  /** Auto-detect type from connected source node's output slot */
  onConnectionsChange(
    type: number,
    slotIndex: number,
    isConnected: boolean,
    _link_info: any,
    _input_info: any,
  ): void {
    if (type !== 1 || slotIndex !== 0) return; // only care about input slot 0
    if (!isConnected || !this.graph) return;

    const inp = this.inputs[0];
    if (!inp || inp.link == null) return;

    const links = (this.graph as any).links;
    const link = links?.[inp.link];
    if (!link) return;

    const sourceNode = this.graph.getNodeById(link.origin_id);
    if (!sourceNode?.outputs) return;

    const sourceSlot = sourceNode.outputs[link.origin_slot];
    if (!sourceSlot?.type) return;

    const detectedType = String(sourceSlot.type);
    // Only auto-set if it's a known concrete type (not dynamic/object)
    if (detectedType && detectedType !== 'dynamic' && detectedType !== 'object' && detectedType !== '*') {
      this.properties['outputType'] = detectedType;
      if (this.typeWidget)
        this.typeWidget.value = detectedType;
    }
  }
}

/** Register all output node types */
export function registerOutputNodes(): void {
  LiteGraph.registerNodeType('Outputs/Table Output', TableOutputNode);
  LiteGraph.registerNodeType('Outputs/Value Output', ValueOutputNode);
}
