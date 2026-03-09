import {LGraphNode, LiteGraph} from 'litegraph.js';
import {getSlotColor} from '../types/type-map';

// --- Equals (==) ---

class EqualsNode extends LGraphNode {
  static title = 'Equals (==)';
  static desc = 'Checks if two values are equal';
  dgNodeType = 'utility';

  constructor() {
    super('Equals (==)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Not Equals (!=) ---

class NotEqualsNode extends LGraphNode {
  static title = 'Not Equals (!=)';
  static desc = 'Checks if two values are not equal';
  dgNodeType = 'utility';

  constructor() {
    super('Not Equals (!=)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Greater Than (>) ---

class GreaterThanNode extends LGraphNode {
  static title = 'Greater Than (>)';
  static desc = 'Checks if a > b';
  dgNodeType = 'utility';

  constructor() {
    super('Greater Than (>)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Greater Or Equal (>=) ---

class GreaterOrEqualNode extends LGraphNode {
  static title = 'Greater Or Equal (>=)';
  static desc = 'Checks if a >= b';
  dgNodeType = 'utility';

  constructor() {
    super('Greater Or Equal (>=)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 180);
  }
}

// --- Less Than (<) ---

class LessThanNode extends LGraphNode {
  static title = 'Less Than (<)';
  static desc = 'Checks if a < b';
  dgNodeType = 'utility';

  constructor() {
    super('Less Than (<)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- Less Or Equal (<=) ---

class LessOrEqualNode extends LGraphNode {
  static title = 'Less Or Equal (<=)';
  static desc = 'Checks if a <= b';
  dgNodeType = 'utility';

  constructor() {
    super('Less Or Equal (<=)');
    this.properties = {};

    const inA = this.addInput('a', 'dynamic');
    inA.color_on = getSlotColor('dynamic');
    inA.color_off = getSlotColor('dynamic');

    const inB = this.addInput('b', 'dynamic');
    inB.color_on = getSlotColor('dynamic');
    inB.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 180);
  }
}

// --- String Contains ---

class StringContainsNode extends LGraphNode {
  static title = 'Contains';
  static desc = 'Checks if a string contains a substring';
  dgNodeType = 'utility';

  constructor() {
    super('Contains');
    this.properties = {};

    const inStr = this.addInput('text', 'string');
    inStr.color_on = getSlotColor('string');
    inStr.color_off = getSlotColor('string');

    const inSub = this.addInput('substring', 'string');
    inSub.color_on = getSlotColor('string');
    inSub.color_off = getSlotColor('string');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- String StartsWith ---

class StartsWithNode extends LGraphNode {
  static title = 'Starts With';
  static desc = 'Checks if a string starts with a prefix';
  dgNodeType = 'utility';

  constructor() {
    super('Starts With');
    this.properties = {};

    const inStr = this.addInput('text', 'string');
    inStr.color_on = getSlotColor('string');
    inStr.color_off = getSlotColor('string');

    const inPrefix = this.addInput('prefix', 'string');
    inPrefix.color_on = getSlotColor('string');
    inPrefix.color_off = getSlotColor('string');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- String EndsWith ---

class EndsWithNode extends LGraphNode {
  static title = 'Ends With';
  static desc = 'Checks if a string ends with a suffix';
  dgNodeType = 'utility';

  constructor() {
    super('Ends With');
    this.properties = {};

    const inStr = this.addInput('text', 'string');
    inStr.color_on = getSlotColor('string');
    inStr.color_off = getSlotColor('string');

    const inSuffix = this.addInput('suffix', 'string');
    inSuffix.color_on = getSlotColor('string');
    inSuffix.color_off = getSlotColor('string');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 160);
  }
}

// --- IsNull ---

class IsNullNode extends LGraphNode {
  static title = 'Is Null';
  static desc = 'Checks if a value is null or undefined';
  dgNodeType = 'utility';

  constructor() {
    super('Is Null');
    this.properties = {};

    const inVal = this.addInput('value', 'dynamic');
    inVal.color_on = getSlotColor('dynamic');
    inVal.color_off = getSlotColor('dynamic');

    const out = this.addOutput('result', 'bool');
    out.color_on = getSlotColor('bool');
    out.color_off = getSlotColor('bool');

    this.color = '#AB47BC';
    this.bgcolor = '#ffffff';
    this.size = this.computeSize();
    this.size[0] = Math.max(this.size[0], 140);
  }
}

/** Register all comparison node types */
export function registerComparisonNodes(): void {
  LiteGraph.registerNodeType('Comparisons/Equals (==)', EqualsNode);
  LiteGraph.registerNodeType('Comparisons/Not Equals (!=)', NotEqualsNode);
  LiteGraph.registerNodeType('Comparisons/Greater Than (>)', GreaterThanNode);
  LiteGraph.registerNodeType('Comparisons/Greater Or Equal (>=)', GreaterOrEqualNode);
  LiteGraph.registerNodeType('Comparisons/Less Than (<)', LessThanNode);
  LiteGraph.registerNodeType('Comparisons/Less Or Equal (<=)', LessOrEqualNode);
  LiteGraph.registerNodeType('Comparisons/Contains', StringContainsNode);
  LiteGraph.registerNodeType('Comparisons/Starts With', StartsWithNode);
  LiteGraph.registerNodeType('Comparisons/Ends With', EndsWithNode);
  LiteGraph.registerNodeType('Comparisons/Is Null', IsNullNode);
}
