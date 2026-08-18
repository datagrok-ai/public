/** Output nodes — emit `//output:` lines; each has one input socket for the final result. */

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
    // Dynamic input; the declared type comes from `outputType`, auto-updated on connect.
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
  }
}
