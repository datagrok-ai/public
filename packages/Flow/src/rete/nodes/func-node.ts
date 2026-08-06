/** FuncNode — one node per DG.Func: an input slot per parameter, pass-through outputs
 *  mirroring each input (execution ordering for mutators), then the real outputs. */

import {ClassicPreset} from 'rete';
import * as DG from 'datagrok-api/dg';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {
  dgTypeToSlotType, getNodeColors, categorizeBySignature, domainCategory, isStringListType,
} from '../../types/type-map';
import {
  getRole, getPackageName, getFuncQualifiedName, getFuncDisplayName, isInputOptional,
  getParamDescription, getParamDisplayName, getParamDefault,
} from '../../utils/dart-proxy-utils';
import {
  nodeHiddenInputsOf, hiddenOutputsOf, funcWrapperOf, wrapperProperties, funcValidatorOf, inputCaptionOf,
} from '../../utils/func-input-overrides';
import {isLiteralChoiceList} from '../../utils/choice-refs';

const PRIMITIVE_DEFAULTS: Record<string, unknown> = {
  string: '',
  int: 0,
  double: 0,
  num: 0,
  bool: false,
};

/** A non-nullable choice input has no empty option — DG renders the first choice, so seeding it
 *  keeps stored and shown in agreement. Skipped for reference choices (`Pkg:func()`), which resolve later. */
export function impliedChoiceDefault(prop: DG.Property): string | undefined {
  try {
    if (String(prop.propertyType) !== 'string' || isInputOptional(prop)) return undefined;
    const choices: unknown = prop.choices;
    if (!Array.isArray(choices) || !isLiteralChoiceList(choices)) return undefined;
    const first = choices.map((c) => String(c)).find((c) => c.length > 0);
    return first;
  } catch {
    return undefined; // Dart proxy access can throw
  }
}

/** The dataframe input a column param defaults to: same numeric suffix (`keys2` → `table2`), else the first. */
export function defaultTableParam(columnParam: string, dataframeParams: string[]): string {
  const suffix = /(\d+)$/.exec(columnParam)?.[1];
  if (suffix !== undefined) {
    const matched = dataframeParams.find((d) => d.endsWith(suffix));
    if (matched) return matched;
  }
  return dataframeParams[0];
}

export class FuncNode extends FlowNode {
  constructor(func: DG.Func) {
    const role = getRole(func);
    // A wrapped function builds the node from the wrapper's exposed inputs — real DG.Property
    // objects, so everything goes through the same code path as real params.
    const wrapper = funcWrapperOf(func);
    const effectiveInputs = wrapper ? wrapperProperties(wrapper) : func.inputs;
    const inputTypes = effectiveInputs.map((p) => String(p.propertyType));
    const category = domainCategory(getPackageName(func), inputTypes) ?? categorizeBySignature(
      inputTypes,
      func.outputs.map((p) => String(p.propertyType)),
      role);
    const colors = getNodeColors(role, func.name, category);
    const qualifiedName = getFuncQualifiedName(func);
    const displayName = getFuncDisplayName(func) || func.name;

    super(displayName);
    this.dgNodeType = 'func';
    this.dgFunc = func;
    this.dgFuncName = qualifiedName;
    this.dgRole = role;
    this.dgPackageName = getPackageName(func);
    (this as unknown as {color: string; bgcolor: string}).color = colors.color;
    (this as unknown as {color: string; bgcolor: string}).bgcolor = colors.bgcolor;

    // Hidden inputs stay fully data-carrying — only rendering skips them.
    this.hiddenInputs = nodeHiddenInputsOf(func);
    this.hiddenOutputs = hiddenOutputsOf(func);
    this.extraValidator = funcValidatorOf(func);
    this.funcWrapper = wrapper ?? undefined;
    const funcInputs = effectiveInputs;
    const funcOutputs = func.outputs;

    const dataframeParams = funcInputs
      .filter((p) => String(p.propertyType) === 'dataframe')
      .map((p) => p.name);

    // Slot key = property name (identity); the label shows the declared caption (display only).
    for (const inp of funcInputs) {
      const slotType = dgTypeToSlotType(inp.propertyType);
      this.addInput(inp.name, new ClassicPreset.Input(getSocket(slotType),
        inputCaptionOf(func, inp.name) ?? getParamDisplayName(inp)));
      const inpDesc = getParamDescription(inp);
      if (inpDesc) this.inputDescriptions[inp.name] = inpDesc;

      if (inp.propertyType in PRIMITIVE_DEFAULTS) {
        // String-encoded defaults are coerced to the declared type ('false' must not compile to `true`).
        // A REQUIRED numeric with no declared default seeds null, not 0 — a zero would read
        // as "set" and hide the missing requirement.
        const declared = getParamDefault(inp) ?? impliedChoiceDefault(inp);
        let def = declared ?? PRIMITIVE_DEFAULTS[inp.propertyType];
        if (inp.propertyType === 'bool' && typeof def === 'string')
          def = def.toLowerCase() === 'true';
        else if (typeof def === 'string' && def !== '' &&
                 ['int', 'double', 'num'].includes(String(inp.propertyType)) && !isNaN(Number(def)))
          def = Number(def);
        if (declared === undefined && !isInputOptional(inp) &&
            ['int', 'double', 'num'].includes(String(inp.propertyType)))
          def = null;
        this.inputValues[inp.name] = def;
      } else if ((inp.propertyType === 'column' || inp.propertyType === 'column_list') &&
                 dataframeParams.length > 0) {
        // Unconnected column values compile to `table.col(...)` against the associated
        // dataframe input; table-less funcs stay connection-only.
        this.inputValues[inp.name] = '';
        if (!this.properties['columnTables']) this.properties['columnTables'] = {};
        (this.properties['columnTables'] as Record<string, string>)[inp.name] =
          defaultTableParam(inp.name, dataframeParams);
      } else if (isStringListType(inp.propertyType)) {
        // Comma-separated inline value; the compiler emits a trimmed JS array.
        this.inputValues[inp.name] = '';
      } else if (String(inp.propertyType) === 'list') {
        // Plain `list` (incl. `list<string>` params) edits as a native DG List input holding a JS array.
        this.inputValues[inp.name] = [];
      }
    }

    // Pass-through outputs, one per input in input order — labeled "<Input> →" because
    // an anonymous `→` made users hunt for which row carries the table onward.
    this.passthroughCount = funcInputs.length;
    for (const inp of funcInputs) {
      const slotType = dgTypeToSlotType(inp.propertyType);
      const ptKey = `${inp.name}__pt`;
      this.addOutput(ptKey, new ClassicPreset.Output(getSocket(slotType),
        `${inputCaptionOf(func, inp.name) ?? getParamDisplayName(inp)} →`));
    }

    // Real outputs after the pass-throughs.
    for (const out of funcOutputs) {
      const slotType = dgTypeToSlotType(out.propertyType);
      this.addOutput(out.name, new ClassicPreset.Output(getSocket(slotType), out.name));
      const outDesc = getParamDescription(out);
      if (outDesc) this.outputDescriptions[out.name] = outDesc;
    }

    // Annotation override: a leading (dataframe, column) pair is forced required even when
    // declared nullable — platform mis-annotation is rampant, and such a function is
    // meaningless without its data.
    const leadingTableColumn = funcInputs.length >= 2 &&
      String(funcInputs[0].propertyType) === 'dataframe' &&
      String(funcInputs[1].propertyType) === 'column';
    this.requiredInputs = funcInputs
      .filter((p, i) => {
        if (leadingTableColumn && i < 2) return true;
        if (isInputOptional(p)) return false;
        const t = String(p.propertyType);
        if (t === 'bool' || t === 'list' || isStringListType(p.propertyType)) return false;
        return getParamDefault(p) === undefined;
      })
      .map((p) => p.name);
  }

  static passthroughInputName(ptKey: string): string | null {
    return ptKey.endsWith('__pt') ? ptKey.slice(0, -'__pt'.length) : null;
  }
}
