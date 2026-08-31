/** Breakpoint node — pure pass-through; in debug mode the emitted code fires `breakpoint-hit`
 *  and awaits a continue event, in normal runs it's skipped. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';

export class BreakpointNode extends FlowNode {
  constructor() {
    super('Breakpoint');
    this.dgNodeType = 'utility';
    this.properties = {enabled: true, label: ''};
    (this as unknown as {color: string}).color = categoricalColor(CAT.red);
    this.addInput('in', new ClassicPreset.Input(getSocket('dynamic'), 'in'));
    this.addOutput('out', new ClassicPreset.Output(getSocket('dynamic'), 'out'));
  }
}
