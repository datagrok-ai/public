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

    this.color = '#EF5350';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

// --- Value Output ---

class ValueOutputNode extends LGraphNode {
  static title = 'Value Output';
  static desc = 'Marks a result value for the script output';
  dgNodeType = 'output';

  constructor() {
    super('Value Output');
    this.properties = {paramName: 'result', outputType: 'double', description: ''};

    const slot = this.addInput('value', 'dynamic');
    slot.color_on = getSlotColor('dynamic');
    slot.color_off = getSlotColor('dynamic');

    this.color = '#EF5350';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }

  /** Auto-detect type from connected source node's output slot */
  onConnectionsChange(
    type: number,
    slotIndex: number,
    isConnected: boolean,
    _link_info: any,
    _input_info: any,
  ): void {
    if (type !== 1 || slotIndex !== 0) return;
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
    if (detectedType && detectedType !== 'dynamic' && detectedType !== 'object' && detectedType !== '*')
      this.properties['outputType'] = detectedType;
  }
}

/** Register all output node types */
export function registerOutputNodes(): void {
  LiteGraph.registerNodeType('Outputs/Table Output', TableOutputNode);
  LiteGraph.registerNodeType('Outputs/Value Output', ValueOutputNode);
}
