import {LGraphNode, IWidget} from 'litegraph.js';
import * as DG from 'datagrok-api/dg';
import {dgTypeToSlotType, getNodeColors, getSlotColor} from '../types/type-map';
import {getRole, getFuncQualifiedName} from '../utils/dart-proxy-utils';

/** Primitive types that get inline widgets on nodes */
const WIDGET_TYPES: Record<string, {widgetType: string; defaultVal: any}> = {
  'string': {widgetType: 'text', defaultVal: ''},
  'int': {widgetType: 'number', defaultVal: 0},
  'double': {widgetType: 'number', defaultVal: 0},
  'num': {widgetType: 'number', defaultVal: 0},
  'bool': {widgetType: 'toggle', defaultVal: false},
};

/**
 * Creates a LiteGraph node constructor for a DG.Func.
 * Each function becomes its own registered node type.
 */
export function createFuncNodeClass(func: DG.Func): {new(): LGraphNode} {
  const role = getRole(func);
  const colors = getNodeColors(role);
  const qualifiedName = getFuncQualifiedName(func);
  const funcInputs = func.inputs;
  const funcOutputs = func.outputs;

  class FuncNode extends LGraphNode {
    static title = func.name;
    static desc = func.description || '';

    dgFunc: DG.Func;
    dgFuncName: string;
    dgRole: string | null;
    inputWidgets: Record<string, IWidget>;

    constructor() {
      super(func.name);
      this.dgFunc = func;
      this.dgFuncName = qualifiedName;
      this.dgRole = role;
      this.inputWidgets = {};

      this.color = colors.color;
      this.bgcolor = colors.bgcolor;

      for (const inp of funcInputs) {
        const slotType = dgTypeToSlotType(inp.propertyType);
        const slot = this.addInput(inp.name, slotType);
        slot.color_on = getSlotColor(slotType);
        slot.color_off = getSlotColor(slotType);

        const wInfo = WIDGET_TYPES[inp.propertyType];
        if (wInfo) {
          const defaultVal = inp.defaultValue ?? wInfo.defaultVal;
          const propKey = `_input_${inp.name}`;
          const opts: any = {property: propKey};
          if (wInfo.widgetType === 'number')
            opts.precision = inp.propertyType === 'int' ? 0 : 3;
          const w = this.addWidget(
            wInfo.widgetType as IWidget['type'],
            inp.name,
            defaultVal,
            (v: any) => {this.properties[propKey] = v;},
            opts,
          );
          this.inputWidgets[inp.name] = w;
          this.properties[`_input_${inp.name}`] = defaultVal;
        }
      }

      for (const out of funcOutputs) {
        const slotType = dgTypeToSlotType(out.propertyType);
        const slot = this.addOutput(out.name, slotType);
        slot.color_on = getSlotColor(slotType);
        slot.color_off = getSlotColor(slotType);
      }

      this.size = this.computeSize();
      this.size[0] = Math.max(this.size[0], 180);
    }

    getInputValue(name: string, slotIndex: number): any {
      if (this.isInputConnected(slotIndex))
        return this.getInputData(slotIndex);
      return this.properties[`_input_${name}`];
    }

    hasHardcodedValue(name: string): boolean {
      return this.properties[`_input_${name}`] !== undefined;
    }

    onConnectionsChange(
      type: number,
      slotIndex: number,
      isConnected: boolean,
    ): void {
      if (type === 1) {
        const inp = this.inputs[slotIndex];
        if (inp && this.inputWidgets[inp.name])
          (this.inputWidgets[inp.name] as any)._hidden = isConnected;
      }
    }
  }

  return FuncNode;
}
