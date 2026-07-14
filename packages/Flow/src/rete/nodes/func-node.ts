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
import {
  dgTypeToSlotType, getNodeColors, categorizeBySignature, domainCategory, isStringListType,
} from '../../types/type-map';
import {
  getRole, getPackageName, getFuncQualifiedName, getFuncDisplayName, isInputOptional,
  getParamDescription, getParamDisplayName, getParamDefault,
} from '../../utils/dart-proxy-utils';
import {hiddenInputsOf} from '../../utils/func-input-overrides';

const PRIMITIVE_DEFAULTS: Record<string, unknown> = {
  string: '',
  int: 0,
  double: 0,
  num: 0,
  bool: false,
};

/** Pick the dataframe input a column/column-list param defaults to: the one
 *  sharing its numeric suffix (JoinTables: `keys2` → `table2`), else the first
 *  dataframe input. Used to seed the `columnTables` association on a fresh node
 *  and as the compiler's fallback for graphs that pre-date the association. */
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
    const inputTypes = func.inputs.map((p) => String(p.propertyType));
    // Domain (chem/bio) wins over the signature-based task category — but only
    // for operations on data (not pure sources/queries), matching the toolbox
    // grouping — so a cheminformatics/bioinformatics node reads its domain from
    // its color.
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

    // Hidden inputs (HIDDEN_FUNC_INPUTS) stay fully data-carrying — socket,
    // seeded value, pass-through, compile, script import/emit — but the node
    // component and the property panel don't render them.
    this.hiddenInputs = hiddenInputsOf(func);
    const funcInputs = func.inputs;
    const funcOutputs = func.outputs;

    // Dataframe input param names — column/column-list inputs resolve their
    // `table.col(...)` against one of these (see `columnTables` below).
    const dataframeParams = funcInputs
      .filter((p) => String(p.propertyType) === 'dataframe')
      .map((p) => p.name);

    // 1. Inputs. The slot key is the property name (identity — used for
    // connections, inputValues, compilation); the label shows the caption when
    // one is declared (display only).
    for (const inp of funcInputs) {
      const slotType = dgTypeToSlotType(inp.propertyType);
      this.addInput(inp.name, new ClassicPreset.Input(getSocket(slotType), getParamDisplayName(inp)));
      const inpDesc = getParamDescription(inp);
      if (inpDesc) this.inputDescriptions[inp.name] = inpDesc;

      if (inp.propertyType in PRIMITIVE_DEFAULTS) {
        // Seed the declared default (`defaultValue ?? initialValue`, unquoted),
        // else the type's zero value. String-encoded defaults are coerced to
        // the declared type so the compiler emits correct literals ('false'
        // must not compile to `true`). A REQUIRED numeric with no declared
        // default seeds null, not 0 — a zero would read as "set" and hide the
        // missing requirement (a blank string stays detectable as-is; the key
        // must exist either way or the panel renders no editor for it).
        const declared = getParamDefault(inp);
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
        // Column / column-list inputs are editable in the property panel as a
        // column name (or comma-separated list). When unconnected, the compiler
        // turns the value into `table.col(...)` against the associated dataframe
        // input — so users don't need a separate Select Column(s) node. We only
        // enable this when the func has a dataframe input to resolve against;
        // the table-less case stays connection-only (a marginal scenario).
        this.inputValues[inp.name] = '';
        if (!this.properties['columnTables']) this.properties['columnTables'] = {};
        (this.properties['columnTables'] as Record<string, string>)[inp.name] =
          defaultTableParam(inp.name, dataframeParams);
      } else if (isStringListType(inp.propertyType)) {
        // string_list / list<string> are editable inline as a comma-separated
        // string; the compiler turns the value into a JS array of trimmed,
        // non-empty strings (so users needn't wire a String List Input node).
        this.inputValues[inp.name] = '';
      } else if (String(inp.propertyType) === 'list') {
        // Plain `list` (incl. `list<string>` params, which the platform reports
        // as propertyType 'list' + propertySubType 'string') is edited with a
        // native DG List input (`ui.input.forProperty`) whose value is a JS
        // array; the compiler emits it as an array literal (empty → omitted).
        this.inputValues[inp.name] = [];
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
      const outDesc = getParamDescription(out);
      if (outDesc) this.outputDescriptions[out.name] = outDesc;
    }

    // Every input that isn't optional (Dart `isOptional` / nullable /
    // `{optional: true}`) and has no declared default must be satisfied —
    // connected, or filled in the panel — before the node can run; drives the
    // "Needs input" hint and every run gate (runnable set, live runs, rerun).
    // Exempt: bool (a checkbox always holds a value) and list-likes (an empty
    // list is a value).
    this.requiredInputs = funcInputs
      .filter((p) => {
        if (isInputOptional(p)) return false;
        const t = String(p.propertyType);
        if (t === 'bool' || t === 'list' || isStringListType(p.propertyType)) return false;
        return getParamDefault(p) === undefined;
      })
      .map((p) => p.name);
  }

  /** Look up the underlying input name corresponding to a pass-through key. */
  static passthroughInputName(ptKey: string): string | null {
    return ptKey.endsWith('__pt') ? ptKey.slice(0, -'__pt'.length) : null;
  }
}
