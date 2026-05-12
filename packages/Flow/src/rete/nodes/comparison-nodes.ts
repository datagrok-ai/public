/** Comparison and logical operator nodes — emit inline JS expressions. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';

const COLOR_COMPARE = '#AB47BC';

abstract class CompareBase extends FlowNode {
  constructor(label: string, lType: string, rType: string, lLabel: string, rLabel: string) {
    super(label);
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_COMPARE;
    this.addInput(lLabel, new ClassicPreset.Input(getSocket(lType), lLabel));
    this.addInput(rLabel, new ClassicPreset.Input(getSocket(rType), rLabel));
    this.addOutput('result', new ClassicPreset.Output(getSocket('bool'), 'result'));
  }
}

export class EqualsNode extends CompareBase {
  constructor() { super('Equals (==)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class NotEqualsNode extends CompareBase {
  constructor() { super('Not Equals (!=)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class GreaterThanNode extends CompareBase {
  constructor() { super('Greater Than (>)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class GreaterOrEqualNode extends CompareBase {
  constructor() { super('Greater Or Equal (>=)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class LessThanNode extends CompareBase {
  constructor() { super('Less Than (<)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class LessOrEqualNode extends CompareBase {
  constructor() { super('Less Or Equal (<=)', 'dynamic', 'dynamic', 'a', 'b'); }
}
export class ContainsNode extends CompareBase {
  constructor() { super('Contains', 'string', 'string', 'text', 'substring'); }
}
export class StartsWithNode extends CompareBase {
  constructor() { super('Starts With', 'string', 'string', 'text', 'prefix'); }
}
export class EndsWithNode extends CompareBase {
  constructor() { super('Ends With', 'string', 'string', 'text', 'suffix'); }
}

export class IsNullNode extends FlowNode {
  constructor() {
    super('Is Null');
    this.dgNodeType = 'utility';
    (this as unknown as {color: string}).color = COLOR_COMPARE;
    this.addInput('value', new ClassicPreset.Input(getSocket('dynamic'), 'value'));
    this.addOutput('result', new ClassicPreset.Output(getSocket('bool'), 'result'));
  }
}
