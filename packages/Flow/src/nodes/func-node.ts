import {LiteGraph, LGraphNode} from 'litegraph.js';
import * as DG from 'datagrok-api/dg';
import {dgTypeToSlotType, getNodeColors, getSlotColor} from '../types/type-map';
import {getRole, getFuncQualifiedName, getFuncDisplayName} from '../utils/dart-proxy-utils';

/** Primitive types that get default values stored in properties */
const PRIMITIVE_DEFAULTS: Record<string, any> = {
  'string': '',
  'int': 0,
  'double': 0,
  'num': 0,
  'bool': false,
};

/**
 * Creates a LiteGraph node constructor for a DG.Func.
 * Each function becomes its own registered node type.
 */
export function createFuncNodeClass(func: DG.Func): {new(): LGraphNode} {
  const role = getRole(func);
  const colors = getNodeColors(role);
  const qualifiedName = getFuncQualifiedName(func);
  const displayName = getFuncDisplayName(func);
  const funcInputs = func.inputs;
  const funcOutputs = func.outputs;

  class FuncNode extends LGraphNode {
    static title = displayName;
    static desc = func.description || '';

    dgFunc: DG.Func;
    dgFuncName: string;
    dgRole: string | null;

    constructor() {
      super(displayName);
      this.dgFunc = func;
      this.dgFuncName = qualifiedName;
      this.dgRole = role;

      this.color = colors.color;
      this.bgcolor = colors.bgcolor;

      for (const inp of funcInputs) {
        const slotType = dgTypeToSlotType(inp.propertyType);
        const slot = this.addInput(inp.name, slotType);
        slot.color_on = getSlotColor(slotType);
        slot.color_off = getSlotColor(slotType);

        // Store default values in properties for primitive types
        if (inp.propertyType in PRIMITIVE_DEFAULTS) {
          const defaultVal = inp.defaultValue ?? PRIMITIVE_DEFAULTS[inp.propertyType];
          this.properties[`_input_${inp.name}`] = defaultVal;
        }
      }

      // Pass-through outputs first: mirror each input as an output for ordering control.
      // Placed first so they visually align with the corresponding input slots.
      // When a function mutates its input (e.g. addNewColumn modifies a table),
      // connecting the pass-through output to the next node enforces execution order.
      this.properties['_passthroughCount'] = funcInputs.length;
      for (const inp of funcInputs) {
        const slotType = dgTypeToSlotType(inp.propertyType);
        const slot = this.addOutput(`${inp.name} \u2192`, slotType);
        slot.color_on = getSlotColor(slotType);
        slot.color_off = getSlotColor(slotType);
      }

      // Real outputs after pass-throughs, with arrow shape to distinguish them
      for (const out of funcOutputs) {
        const slotType = dgTypeToSlotType(out.propertyType);
        const slot = this.addOutput(out.name, slotType);
        slot.color_on = getSlotColor(slotType);
        slot.color_off = getSlotColor(slotType);
        slot.shape = LiteGraph.SQUARE_SHAPE;
      }

      this.size = this.computeSize();
      this.size[0] = Math.max(this.size[0], 140);
    }

    getInputValue(name: string, slotIndex: number): any {
      if (this.isInputConnected(slotIndex))
        return this.getInputData(slotIndex);
      return this.properties[`_input_${name}`];
    }

    hasHardcodedValue(name: string): boolean {
      return this.properties[`_input_${name}`] !== undefined;
    }
  }

  return FuncNode;
}
