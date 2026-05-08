/** Builds a FuncNode instance per DG.Func.
 *
 * Each func node has:
 *   - one input slot per DG parameter, typed via `dgTypeToSlotType`
 *   - **pass-through outputs** mirroring each input (label `→`), used to
 *     enforce execution ordering for mutating functions
 *   - real outputs after the pass-throughs, one per declared DG output
 *
 * Default values for primitive (string/int/double/bool) inputs are stored on
 * the node's `inputValues` map and shown in the property panel — never as
 * inline widgets, so connecting a slot replaces them at compile time. */

import {ClassicPreset} from 'rete';
import * as DG from 'datagrok-api/dg';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {dgTypeToSlotType, getNodeColors} from '../../types/type-map';
import {getRole, getFuncQualifiedName, getFuncDisplayName} from '../../utils/dart-proxy-utils';

const PRIMITIVE_DEFAULTS: Record<string, unknown> = {
  string: '',
  int: 0,
  double: 0,
  num: 0,
  bool: false,
};

export class FuncNode extends FlowNode {
  constructor(func: DG.Func) {
    const role = getRole(func);
    const colors = getNodeColors(role);
    const qualifiedName = getFuncQualifiedName(func);
    const displayName = getFuncDisplayName(func) || func.name;

    super(displayName);
    this.dgNodeType = 'func';
    this.dgFunc = func;
    this.dgFuncName = qualifiedName;
    this.dgRole = role;
    (this as unknown as {color: string; bgcolor: string}).color = colors.color;
    (this as unknown as {color: string; bgcolor: string}).bgcolor = colors.bgcolor;

    const funcInputs = func.inputs;
    const funcOutputs = func.outputs;

    // 1. Inputs.
    for (const inp of funcInputs) {
      const slotType = dgTypeToSlotType(inp.propertyType);
      this.addInput(inp.name, new ClassicPreset.Input(getSocket(slotType), inp.name));

      if (inp.propertyType in PRIMITIVE_DEFAULTS) {
        const def = (inp as unknown as {defaultValue?: unknown}).defaultValue ?? PRIMITIVE_DEFAULTS[inp.propertyType];
        this.inputValues[inp.name] = def;
      }
    }

    // 2. Pass-through outputs — one per input, in input order, label `→`.
    this.passthroughCount = funcInputs.length;
    for (const inp of funcInputs) {
      const slotType = dgTypeToSlotType(inp.propertyType);
      const ptKey = `${inp.name}__pt`;
      this.addOutput(ptKey, new ClassicPreset.Output(getSocket(slotType), '→'));
    }

    // 3. Real outputs after the pass-throughs.
    for (const out of funcOutputs) {
      const slotType = dgTypeToSlotType(out.propertyType);
      this.addOutput(out.name, new ClassicPreset.Output(getSocket(slotType), out.name));
    }
  }

  /** Look up the underlying input name corresponding to a pass-through key. */
  static passthroughInputName(ptKey: string): string | null {
    return ptKey.endsWith('__pt') ? ptKey.slice(0, -'__pt'.length) : null;
  }
}
