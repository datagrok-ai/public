/** Breakpoint node — pauses script execution in debug mode */
import {LiteGraph, LGraphNode} from 'litegraph.js';
import {dgTypeToSlotType, getSlotColor} from '../types/type-map';

export class BreakpointNode extends LGraphNode {
  static title = 'Breakpoint';
  static desc = 'Pauses execution in debug mode until Continue is clicked';

  dgNodeType = 'utility';

  constructor() {
    super('Breakpoint');
    const slotType = dgTypeToSlotType('dynamic');
    const slotColor = getSlotColor(slotType);

    const inp = this.addInput('in', slotType);
    inp.color_on = slotColor;
    inp.color_off = slotColor;

    const out = this.addOutput('out', slotType);
    out.color_on = slotColor;
    out.color_off = slotColor;

    this.color = '#F44336';
    this.bgcolor = '#ffffff';
    this.properties = {enabled: true, label: ''};
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 120);
  }
}

export function registerBreakpointNode(): void {
  LiteGraph.registerNodeType('Debug/Breakpoint', BreakpointNode);
}
