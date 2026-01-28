import { CONTROL_EXPR, LOOP } from "./constants";
import { ARG_INPUT_KEYS, getIVP, getViewersSpec } from "./scripting-tools";
import { getLookupsInfoData, getOptions } from "./shared-utils";

declare var code: string;
declare var result: any;

const ivp = getIVP(code);
const lookupsOptions = ivp.inputsLookup ? [getLookupsInfoData(ivp.inputsLookup) as any] : [];
const argOptions = ARG_INPUT_KEYS.map((key) => getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG));
const initsOptions = [...ivp.inits.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.INITS));
const paramsOptions = ivp.params ? [...ivp.params.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.PARAMS)) : [];
const loopOptions = ivp.loop ? [getOptions(LOOP.COUNT_NAME, ivp.loop.count, CONTROL_EXPR.LOOP)] : [];
const inputs = ([...lookupsOptions, ...argOptions, ...initsOptions, ...paramsOptions, ...loopOptions]);
const outputs = [{
  name: 'df',
  type: 'dataframe',
  options: {
    caption: ivp.name,
    viewer: getViewersSpec(ivp),
  }
}];
result = JSON.stringify({inputs, outputs});
