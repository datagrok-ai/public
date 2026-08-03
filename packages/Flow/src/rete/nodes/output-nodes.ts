/** Output nodes — emit `//output:` annotation lines in the generated script.
 *
 * Each has one input socket; the user wires their final result to it.
 * `ValueOutputNode` lets the user pick the declared output type from a combo
 * in the property panel; the FlowEditor watches new connections and
 * auto-updates this when a non-dynamic source is connected. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';

const COLOR_OUTPUT = categoricalColor(CAT.red);

export class TableOutputNode extends FlowNode {
  constructor() {
    super('Table Output');
    this.dgNodeType = 'output';
    this.dgOutputType = 'dataframe';
    this.properties = {paramName: 'result'};
    (this as unknown as {color: string}).color = COLOR_OUTPUT;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
  }
}

export class ValueOutputNode extends FlowNode {
  constructor() {
    super('Value Output');
    this.dgNodeType = 'output';
    this.properties = {paramName: 'result', outputType: 'double'};
    (this as unknown as {color: string}).color = COLOR_OUTPUT;
    // Input socket is dynamic — accepts anything; declared type comes from
    // the `outputType` property and is auto-updated on connect via
    // `FlowEditor.attachOutputAutoType`.
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
  }
}
