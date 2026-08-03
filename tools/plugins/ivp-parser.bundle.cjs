"use strict";
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};
var __copyProps = (to, from, except, desc) => {
  if (from && typeof from === "object" || typeof from === "function") {
    for (let key of __getOwnPropNames(from))
      if (!__hasOwnProp.call(to, key) && key !== except)
        __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
  }
  return to;
};
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// plugins/ivp-parser.entry.mjs
var ivp_parser_entry_exports = {};
__export(ivp_parser_entry_exports, {
  MAX_LINE_CHART: () => MAX_LINE_CHART,
  STAGE_COL_NAME: () => STAGE_COL_NAME,
  getIVP: () => getIVP
});
module.exports = __toCommonJS(ivp_parser_entry_exports);

// ../libraries/compute-utils/node_modules/diff-grok/dist/src/scripting-tools/constants.js
var CONTROL_TAG = "#";
var CONTROL_TAG_LEN = CONTROL_TAG.length;
var MAX_LINE_CHART = 4;
var CONTROL_EXPR;
(function(CONTROL_EXPR2) {
  CONTROL_EXPR2["NAME"] = "#name";
  CONTROL_EXPR2["TAGS"] = "#tags";
  CONTROL_EXPR2["DESCR"] = "#description";
  CONTROL_EXPR2["DIF_EQ"] = "#equations";
  CONTROL_EXPR2["EXPR"] = "#expressions";
  CONTROL_EXPR2["ARG"] = "#argument";
  CONTROL_EXPR2["INITS"] = "#inits";
  CONTROL_EXPR2["CONSTS"] = "#constants";
  CONTROL_EXPR2["PARAMS"] = "#parameters";
  CONTROL_EXPR2["TOL"] = "#tolerance";
  CONTROL_EXPR2["LOOP"] = "#loop";
  CONTROL_EXPR2["UPDATE"] = "#update";
  CONTROL_EXPR2["RUN_ON_OPEN"] = "#meta.runOnOpen";
  CONTROL_EXPR2["RUN_ON_INPUT"] = "#meta.runOnInput";
  CONTROL_EXPR2["OUTPUT"] = "#output";
  CONTROL_EXPR2["COMMENT"] = "#comment";
  CONTROL_EXPR2["SOLVER"] = "#meta.solver";
  CONTROL_EXPR2["INPUTS"] = "#meta.inputs";
})(CONTROL_EXPR || (CONTROL_EXPR = {}));
var LOOP;
(function(LOOP2) {
  LOOP2[LOOP2["MIN_LINES_COUNT"] = 1] = "MIN_LINES_COUNT";
  LOOP2[LOOP2["COUNT_IDX"] = 0] = "COUNT_IDX";
  LOOP2["COUNT_NAME"] = "_count";
  LOOP2[LOOP2["MIN_COUNT"] = 1] = "MIN_COUNT";
})(LOOP || (LOOP = {}));
var UPDATE;
(function(UPDATE2) {
  UPDATE2[UPDATE2["MIN_LINES_COUNT"] = 1] = "MIN_LINES_COUNT";
  UPDATE2[UPDATE2["DURATION_IDX"] = 0] = "DURATION_IDX";
  UPDATE2["DURATION"] = "_duration";
})(UPDATE || (UPDATE = {}));
var SOLVER_OPTIONS_RANGES = /* @__PURE__ */ new Map([
  ["maxTime", { min: 1, max: 1e4 }],
  ["scale", { min: 0.5, max: 1 }]
]);

// ../libraries/compute-utils/node_modules/diff-grok/dist/src/scripting-tools/model-error.js
var ModelError = class extends Error {
  helpUrl;
  toHighlight = void 0;
  constructor(message, helpUrl, toHighlight) {
    super(message);
    this.helpUrl = helpUrl;
    this.toHighlight = toHighlight;
  }
  getHelpUrl() {
    return this.helpUrl;
  }
  getToHighlight() {
    return this.toHighlight;
  }
};

// ../libraries/compute-utils/node_modules/diff-grok/dist/src/scripting-tools/scripting-tools.js
var CONTROL_SEP = ":";
var COMMA = ",";
var EQUAL_SIGN = "=";
var DIV_SIGN = "/";
var SERVICE = "_";
var BRACE_OPEN = "{";
var BRACE_CLOSE = "}";
var ANNOT_SEPAR = ";";
var DEFAULT_TOL = "0.00005";
var DEFAULT_SOLVER_SETTINGS = "{}";
var COMMENT_SEQ = "//";
var STAGE_COL_NAME = `_Stage`;
var MATH_FUNCS = ["pow", "sin", "cos", "tan", "asin", "acos", "atan", "sqrt", "exp", "log", "sinh", "cosh", "tanh"];
var POW_IDX = MATH_FUNCS.indexOf("pow");
var MATH_CONSTS = ["PI", "E"];
var SCRIPTING;
(function(SCRIPTING2) {
  SCRIPTING2["ARG_NAME"] = "name";
  SCRIPTING2["COUNT"] = "count";
  SCRIPTING2["DURATION"] = "duration";
})(SCRIPTING || (SCRIPTING = {}));
var ERROR_LINK;
(function(ERROR_LINK2) {
  ERROR_LINK2["MAIN_DOCS"] = "/help/compute/diff-studio";
  ERROR_LINK2["CORE_BLOCKS"] = "/help/compute/diff-studio#core-blocks";
  ERROR_LINK2["COMPS_SYNTAX"] = "/help/compute/diff-studio#model-components-and-syntax";
  ERROR_LINK2["LOOP"] = "/help/compute/diff-studio#cyclic-processes";
  ERROR_LINK2["UPDATE"] = "/help/compute/diff-studio#multistage-processes";
  ERROR_LINK2["UI_OPTS"] = "/help/compute/diff-studio#user-interface-options";
  ERROR_LINK2["LOOP_VS_UPDATE"] = "/help/compute/diff-studio#advanced-features";
  ERROR_LINK2["SOLVER_CONFIG"] = "/help/compute/diff-studio#solver-configuration";
  ERROR_LINK2["MODEL_PARAMS"] = "/help/compute/diff-studio#model-parameters";
})(ERROR_LINK || (ERROR_LINK = {}));
var ERROR_MSG;
(function(ERROR_MSG2) {
  ERROR_MSG2["CTRL_EXPR"] = 'Unsupported control expression with the tag **"#"**';
  ERROR_MSG2["ARG"] = "'The **#argument** block must consist of 3 lines specifying initial and final time, and solution grid step.";
  ERROR_MSG2["LOOP"] = "The **#loop** block must contain at least one line.";
  ERROR_MSG2["COUNT"] = "Incorrect loop count";
  ERROR_MSG2["LOOP_VS_UPDATE"] = "The **#loop** and **'#update'** blocks cannot be used simultaneously.";
  ERROR_MSG2["UPDATE_LINES_COUNT"] = "The **'#update'** block must contain at least one line.";
  ERROR_MSG2["DURATION"] = "Incorrect update duration";
  ERROR_MSG2["BRACES"] = " Missing one of the braces (**{**, **}**).";
  ERROR_MSG2["COLON"] = 'Incorrect position of **":"**.';
  ERROR_MSG2["CASE_INSENS"] = "Non-unique name (case-insensitive): use different caption for ";
  ERROR_MSG2["MISSING_INIT"] = "Correct the **#inits** block.";
  ERROR_MSG2["UNDEF_NAME"] = "Model name missing. Specify the model name in the **#name** block.";
  ERROR_MSG2["UNDEF_DEQS"] = "Differential equation(s) are required for this model. Add equation(s) under the **#equations** block.";
  ERROR_MSG2["UNDEF_INITS"] = "Initial conditions are required for this model. Add initial conditions under the **#inits** block.";
  ERROR_MSG2["UNDEF_ARG"] = "Argument specification is required for this model. Specify an argument, its range, and a grid step in the **#argument** block.";
  ERROR_MSG2["CORRECT_ARG_LIM"] = "Correct limits in the **#argument** block.";
  ERROR_MSG2["INTERVAL"] = "Incorrect range for";
  ERROR_MSG2["NEGATIVE_STEP"] = "Solution grid step must be positive. Correct the **#argument** block.";
  ERROR_MSG2["INCOR_STEP"] = "Grid step must less than the length of solution interval. Correct the **#argument** block.";
  ERROR_MSG2["MISS_COLON"] = 'Missing **":"**';
  ERROR_MSG2["NAN"] = "is not a valid number. Correct the line";
  ERROR_MSG2["SERVICE_START"] = 'Variable names must not begin with **"_"**.';
  ERROR_MSG2["REUSE_NAME"] = "Variable reuse (case-insensitive): rename ";
  ERROR_MSG2["SOLVER"] = "Incorrect solver options. Correct the **#meta.solver** line.";
})(ERROR_MSG || (ERROR_MSG = {}));
var ANNOT;
(function(ANNOT2) {
  ANNOT2["NAME"] = "//name:";
  ANNOT2["DESCR"] = "//description:";
  ANNOT2["TAGS"] = "//tags:";
  ANNOT2["LANG"] = "//language: javascript";
  ANNOT2["DOUBLE_INPUT"] = "//input: double";
  ANNOT2["INT_INPUT"] = "//input: int";
  ANNOT2["OUTPUT"] = "//output: dataframe df";
  ANNOT2["EDITOR"] = "//editor: Compute:RichFunctionViewEditor";
  ANNOT2["SIDEBAR"] = "//sidebar: @compute";
  ANNOT2["CAPTION"] = "caption:";
  ANNOT2["ARG_INIT"] = "{caption: Initial; category: Argument}";
  ANNOT2["ARG_FIN"] = "{caption: Final; category: Argument}";
  ANNOT2["ARG_STEP"] = "{caption: Step; category: Argument}";
  ANNOT2["INITS"] = "category: Initial values";
  ANNOT2["PARAMS"] = "category: Parameters";
})(ANNOT || (ANNOT = {}));
var SCRIPT;
(function(SCRIPT2) {
  SCRIPT2["CONSTS"] = "// constants";
  SCRIPT2["ODE_COM"] = "// the problem definition";
  SCRIPT2["ODE"] = "let odes = {";
  SCRIPT2["SOLVER_COM"] = "// solve the problem";
  SCRIPT2["SOLVER"] = "const solver = await grok.functions.eval('DiffStudio:solveEquations');";
  SCRIPT2["PREPARE"] = "let call = solver.prepare({problem: odes, options: opts});";
  SCRIPT2["CALL"] = "await call.call();";
  SCRIPT2["OUTPUT"] = "let df = call.getParamValue('df');";
  SCRIPT2["SPACE2"] = "  ";
  SCRIPT2["SPACE4"] = "    ";
  SCRIPT2["SPACE6"] = "      ";
  SCRIPT2["SPACE8"] = "        ";
  SCRIPT2["FUNC_VALS"] = "// extract function values";
  SCRIPT2["EVAL_EXPR"] = "// evaluate expressions";
  SCRIPT2["COMP_OUT"] = "// compute output";
  SCRIPT2["MATH_FUNC_COM"] = "// used Math-functions";
  SCRIPT2["MATH_CONST_COM"] = "// used Math-constants";
  SCRIPT2["ONE_STAGE_COM"] = "\n// one stage solution";
  SCRIPT2["ONE_STAGE_BEGIN"] = "let _oneStage = async (";
  SCRIPT2["ONE_STAGE_END"] = ") => {";
  SCRIPT2["ASYNC_OUTPUT"] = "let df = await _oneStage(";
  SCRIPT2["RETURN_OUTPUT"] = "return call.getParamValue('df');";
  SCRIPT2["EMPTY_OUTPUT"] = "let df = DG.DataFrame.create();";
  SCRIPT2["APPEND"] = "df.append(";
  SCRIPT2["SOLUTION_DF_COM"] = "// solution dataframe";
  SCRIPT2["LOOP_INTERVAL_COM"] = "// loop interval";
  SCRIPT2["LOOP_INTERVAL"] = "_interval";
  SCRIPT2["LAST_IDX"] = "_lastIdx";
  SCRIPT2["UPDATE_COM"] = "// update ";
  SCRIPT2["CUSTOM_OUTPUT_COM"] = "// create custom output";
  SCRIPT2["CUSTOM_COLUMNS"] = "let _columns = [";
  SCRIPT2["ONE_STAGE"] = "await _oneStage(";
  SCRIPT2["SEGMENT_COM"] = "// add segment category";
})(SCRIPT || (SCRIPT = {}));
function getStartOfProblemDef(lines) {
  const linesCount = lines.length;
  let idx = 0;
  while (!lines[idx].startsWith(CONTROL_TAG)) {
    ++idx;
    if (idx === linesCount)
      throw new ModelError(ERROR_MSG.UNDEF_NAME, ERROR_LINK.CORE_BLOCKS);
  }
  return idx;
}
function getBlocks(lines, start) {
  const linesCount = lines.length;
  let beg = start;
  let idx = start + 1;
  const blocks = [];
  while (true) {
    if (idx === linesCount) {
      blocks.push({ begin: beg, end: idx });
      break;
    }
    if (lines[idx].startsWith(CONTROL_TAG)) {
      blocks.push({ begin: beg, end: idx });
      beg = idx;
    }
    ++idx;
  }
  return blocks;
}
function concatMultilineFormulas(source) {
  const res = [source[0]];
  let idxRes = 0;
  let idxSource = 1;
  const size = source.length;
  while (idxSource < size) {
    const str = source[idxSource];
    if (!str.includes(EQUAL_SIGN))
      res[idxRes] += str;
    else {
      res.push(str);
      ++idxRes;
    }
    ++idxSource;
  }
  return res;
}
function getUsedMathIds(text, mathIds) {
  const res = [];
  const size = mathIds.length;
  for (let i = 0; i < size; ++i) {
    if (text.includes(`${mathIds[i]}`))
      res.push(i);
  }
  return res;
}
function getDifEquations(lines) {
  const deqs = /* @__PURE__ */ new Map();
  const names = [];
  let divIdx = 0;
  let eqIdx = 0;
  for (const line of lines) {
    if (line === void 0)
      throw new ModelError(ERROR_MSG.UNDEF_DEQS, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.DIF_EQ);
    divIdx = line.indexOf(DIV_SIGN);
    eqIdx = line.indexOf(EQUAL_SIGN);
    const name = line.slice(line.indexOf("d") + 1, divIdx).replace(/[() ]/g, "");
    names.push(name);
    deqs.set(name, line.slice(eqIdx + 1).trim());
  }
  return {
    equations: deqs,
    solutionNames: names
  };
}
function getExpressions(lines) {
  const exprs = /* @__PURE__ */ new Map();
  let eqIdx = 0;
  for (const line of lines) {
    if (line === void 0)
      break;
    eqIdx = line.indexOf(EQUAL_SIGN);
    exprs.set(line.slice(0, eqIdx).replaceAll(" ", ""), line.slice(eqIdx + 1).trim());
  }
  return exprs;
}
function getInput(line) {
  const str = line.slice(line.indexOf(EQUAL_SIGN) + 1).trim();
  const braceIdx = str.indexOf(BRACE_OPEN);
  let val;
  if (braceIdx === -1) {
    val = Number(str);
    if (isNaN(val))
      throw new ModelError(`**${str}** ${ERROR_MSG.NAN}
 **${line}**.`, ERROR_LINK.COMPS_SYNTAX, line);
    return {
      value: val,
      annot: null
    };
  }
  val = Number(str.slice(0, braceIdx));
  if (isNaN(val))
    throw new ModelError(`**${str.slice(0, braceIdx)}** ${ERROR_MSG.NAN}: **${line}**`, ERROR_LINK.COMPS_SYNTAX, line);
  return {
    value: val,
    annot: str.slice(braceIdx)
  };
}
function getArg(lines) {
  if (lines.length !== 4)
    throw new ModelError(ERROR_MSG.ARG, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.ARG);
  const sepIdx = lines[0].indexOf(CONTROL_SEP);
  if (sepIdx === -1)
    throw new ModelError(`${ERROR_MSG.MISS_COLON}. Correct the line **${lines[0]}**`, ERROR_LINK.CORE_BLOCKS, lines[0]);
  const commaIdx = lines[0].indexOf(COMMA);
  return {
    name: (commaIdx > -1 ? lines[0].slice(sepIdx + 1, commaIdx) : lines[0].slice(sepIdx + 1)).trim(),
    initial: getInput(lines[1]),
    final: getInput(lines[2]),
    step: getInput(lines[3]),
    updateName: commaIdx > -1 ? lines[0].slice(commaIdx + 1).trim() : void 0
  };
}
function getEqualities(lines, begin, end) {
  const source = concatMultilineFormulas(lines.slice(begin, end));
  const eqs = /* @__PURE__ */ new Map();
  let eqIdx = 0;
  for (const line of source) {
    if (line === void 0)
      break;
    eqIdx = line.indexOf(EQUAL_SIGN);
    eqs.set(line.slice(0, eqIdx).replaceAll(" ", ""), getInput(line));
  }
  return eqs;
}
function getLoop(lines, begin, end) {
  if (begin >= end)
    throw new ModelError(ERROR_MSG.LOOP, ERROR_LINK.LOOP, CONTROL_EXPR.LOOP);
  const source = concatMultilineFormulas(lines.slice(begin, end));
  const size = source.length;
  if (size < LOOP.MIN_LINES_COUNT)
    throw new ModelError(ERROR_MSG.LOOP, ERROR_LINK.LOOP, CONTROL_EXPR.LOOP);
  const count = getInput(source[LOOP.COUNT_IDX]);
  if (count.value < LOOP.MIN_COUNT)
    throw new ModelError(`${ERROR_MSG.COUNT}: **${count.value}**.`, ERROR_LINK.LOOP, source[LOOP.COUNT_IDX]);
  return { count, updates: source.slice(LOOP.COUNT_IDX + 1) };
}
function getUpdate(lines, begin, end) {
  const colonIdx = lines[begin].indexOf(CONTROL_SEP);
  if (colonIdx === -1) {
    throw new ModelError(`${ERROR_MSG.MISS_COLON}. Correct the line: **${lines[begin]}**`, ERROR_LINK.UPDATE, lines[begin]);
  }
  const source = concatMultilineFormulas(lines.slice(begin + 1, end));
  const size = source.length;
  if (size < UPDATE.MIN_LINES_COUNT)
    throw new ModelError(ERROR_MSG.UPDATE_LINES_COUNT, ERROR_LINK.UPDATE, CONTROL_EXPR.UPDATE);
  if (source[UPDATE.DURATION_IDX] == void 0)
    throw new ModelError(ERROR_MSG.UPDATE_LINES_COUNT, ERROR_LINK.UPDATE, CONTROL_EXPR.UPDATE);
  const eqIdx = source[UPDATE.DURATION_IDX].indexOf(EQUAL_SIGN);
  if (eqIdx === -1) {
    throw new ModelError(`Missing **${EQUAL_SIGN}**. Correct the line: **${source[UPDATE.DURATION_IDX]}**`, ERROR_LINK.UPDATE, source[UPDATE.DURATION_IDX]);
  }
  return {
    name: lines[begin].slice(colonIdx + 1),
    durationFormula: source[UPDATE.DURATION_IDX].slice(eqIdx + 1).trim(),
    updates: source.slice(UPDATE.DURATION_IDX + 1)
  };
}
function getOutput(lines, begin, end) {
  const res = /* @__PURE__ */ new Map();
  let line;
  let token;
  let eqIdx;
  let openIdx;
  let closeIdx;
  let colonIdx;
  for (let i = begin; i < end; ++i) {
    line = lines[i].trim();
    if (line.length < 1)
      continue;
    eqIdx = line.indexOf(EQUAL_SIGN);
    openIdx = line.indexOf(BRACE_OPEN);
    closeIdx = line.indexOf(BRACE_CLOSE);
    colonIdx = line.indexOf(CONTROL_SEP);
    if (openIdx * closeIdx <= 0) {
      throw new ModelError(`${ERROR_MSG.BRACES} Correct the line **${line}** in the **${CONTROL_EXPR.OUTPUT}** block.`, ERROR_LINK.UI_OPTS, line);
    }
    if (openIdx * colonIdx <= 0) {
      throw new ModelError(`${ERROR_MSG.COLON} Correct the line **${line}** in the **#output** block.`, ERROR_LINK.UI_OPTS, line);
    }
    if (eqIdx === -1) {
      if (openIdx === -1) {
        token = line.trim();
        res.set(token, {
          caption: token,
          formula: null
        });
      } else {
        token = line.slice(0, openIdx).trim();
        res.set(token, {
          caption: line.slice(colonIdx + 1, closeIdx).trim(),
          formula: null
        });
      }
    } else if (openIdx === -1) {
      token = line.slice(0, eqIdx).trim();
      res.set(token, {
        caption: token,
        formula: line.slice(eqIdx + 1)
      });
    } else {
      token = line.slice(0, eqIdx).trim();
      res.set(token, {
        caption: line.slice(colonIdx + 1, closeIdx),
        formula: line.slice(eqIdx + 1, openIdx)
      });
    }
  }
  return res;
}
function getIVP(text) {
  let name;
  let tags = null;
  let descr = null;
  let deqs;
  let exprs = null;
  let arg;
  let inits;
  let consts = null;
  let params = null;
  let tolerance = DEFAULT_TOL;
  let loop = null;
  const updates = [];
  const metas = [];
  let outputs = null;
  let solverSettings = DEFAULT_SOLVER_SETTINGS;
  let inputsLookup = null;
  const lines = text.replaceAll("	", " ").split("\n").map((s) => {
    const idx = s.indexOf(COMMENT_SEQ);
    return s.slice(0, idx >= 0 ? idx : void 0).trimStart();
  }).filter((s) => s !== "");
  const start = getStartOfProblemDef(lines);
  const blocks = getBlocks(lines, start);
  for (const block of blocks) {
    const firstLine = lines[block.begin];
    if (firstLine.startsWith(CONTROL_EXPR.NAME)) {
      name = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.TAGS)) {
      tags = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.DESCR)) {
      descr = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.DIF_EQ)) {
      deqs = getDifEquations(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    } else if (firstLine.startsWith(CONTROL_EXPR.EXPR)) {
      exprs = getExpressions(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    } else if (firstLine.startsWith(CONTROL_EXPR.ARG)) {
      arg = getArg(concatMultilineFormulas(lines.slice(block.begin, block.end)));
    } else if (firstLine.startsWith(CONTROL_EXPR.INITS)) {
      inits = getEqualities(lines, block.begin + 1, block.end);
    } else if (firstLine.startsWith(CONTROL_EXPR.CONSTS)) {
      consts = getEqualities(lines, block.begin + 1, block.end);
    } else if (firstLine.startsWith(CONTROL_EXPR.PARAMS)) {
      params = getEqualities(lines, block.begin + 1, block.end);
    } else if (firstLine.startsWith(CONTROL_EXPR.TOL)) {
      tolerance = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.SOLVER)) {
      solverSettings = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.LOOP)) {
      loop = getLoop(lines, block.begin + 1, block.end);
    } else if (firstLine.startsWith(CONTROL_EXPR.UPDATE)) {
      updates.push(getUpdate(lines, block.begin, block.end));
    } else if (firstLine.startsWith(CONTROL_EXPR.OUTPUT)) {
      outputs = getOutput(lines, block.begin + 1, block.end);
    } else if (firstLine.startsWith(CONTROL_EXPR.INPUTS)) {
      inputsLookup = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    } else if (firstLine.startsWith(CONTROL_EXPR.COMMENT)) {
    } else
      metas.push(firstLine.slice(CONTROL_TAG_LEN));
  }
  if (loop !== null && updates.length > 0)
    throw new ModelError(ERROR_MSG.LOOP_VS_UPDATE, ERROR_LINK.LOOP_VS_UPDATE, CONTROL_EXPR.LOOP);
  const ivp = {
    name,
    tags,
    descr,
    deqs,
    exprs,
    arg,
    inits,
    consts,
    params,
    tolerance,
    usedMathFuncs: getUsedMathIds(text, MATH_FUNCS),
    usedMathConsts: getUsedMathIds(text, MATH_CONSTS),
    loop,
    updates: updates.length === 0 ? null : updates,
    metas,
    outputs,
    solverSettings,
    inputsLookup
  };
  checkCorrectness(ivp);
  return ivp;
}
function checkSolverSettings(line) {
  const settings = /* @__PURE__ */ new Map();
  const openBraceIdx = line.indexOf(BRACE_OPEN);
  const closeBraceIdx = line.indexOf(BRACE_CLOSE);
  let sepIdx;
  if (openBraceIdx < 0 || closeBraceIdx < 0)
    throw new ModelError(`${ERROR_MSG.BRACES}. Correct the line **${line}**.`, ERROR_LINK.SOLVER_CONFIG, line);
  for (const item of line.slice(openBraceIdx + 1, closeBraceIdx).split(ANNOT_SEPAR)) {
    sepIdx = item.indexOf(CONTROL_SEP);
    if (sepIdx > 1)
      settings.set(item.slice(0, sepIdx).trim(), item.slice(sepIdx + 1).trim());
  }
  SOLVER_OPTIONS_RANGES.forEach((range, opt) => {
    if (settings.has(opt)) {
      const val = Number(settings.get(opt));
      if (val < range.min || val > range.max) {
        throw new ModelError(`${ERROR_MSG.SOLVER}: **${opt}** must be in the range **${range.min}..${range.max}**.`, ERROR_LINK.SOLVER_CONFIG, line);
      }
    }
  });
}
function checkCorrectness(ivp) {
  if (ivp.name === void 0)
    throw new ModelError(ERROR_MSG.UNDEF_NAME, ERROR_LINK.CORE_BLOCKS);
  if (ivp.deqs === void 0 || ivp.deqs.equations.size === 0)
    throw new ModelError(ERROR_MSG.UNDEF_DEQS, ERROR_LINK.CORE_BLOCKS);
  if (ivp.inits === void 0 || ivp.inits.size === 0)
    throw new ModelError(ERROR_MSG.UNDEF_INITS, ERROR_LINK.CORE_BLOCKS);
  if (ivp.arg === void 0)
    throw new ModelError(ERROR_MSG.UNDEF_ARG, ERROR_LINK.CORE_BLOCKS);
  ivp.deqs.equations.forEach((ignore, name) => {
    if (!ivp.inits.has(name)) {
      throw new ModelError(`Initial value for **${name}** is missing. ${ERROR_MSG.MISSING_INIT}`, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.INITS);
    }
  });
  const usedNames = [];
  let lowCase;
  if (ivp.outputs !== null) {
    ivp.outputs.forEach((val) => {
      lowCase = val.caption.toLowerCase();
      if (usedNames.includes(lowCase))
        throw new ModelError(`${ERROR_MSG.CASE_INSENS}**${val.caption}**.`, ERROR_LINK.UI_OPTS, CONTROL_EXPR.OUTPUT);
      else
        usedNames.push(lowCase);
    });
  } else {
    const usedNames2 = [ivp.arg.name];
    ivp.deqs.solutionNames.forEach((name) => {
      lowCase = name.toLowerCase();
      if (usedNames2.includes(lowCase))
        throw new ModelError(`${ERROR_MSG.CASE_INSENS}**${name}**.`, ERROR_LINK.UI_OPTS);
      else
        usedNames2.push(lowCase);
    });
  }
  if (ivp.arg.initial.value >= ivp.arg.final.value) {
    throw new ModelError(`${ERROR_MSG.INTERVAL} **${ivp.arg.name}**. ${ERROR_MSG.CORRECT_ARG_LIM}`, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.ARG);
  }
  if (ivp.arg.step.value <= 0) {
    throw new ModelError(`${ERROR_MSG.NEGATIVE_STEP}`, ERROR_LINK.CORE_BLOCKS);
  }
  if (ivp.arg.step.value > ivp.arg.final.value - ivp.arg.initial.value)
    throw new ModelError(ERROR_MSG.INCOR_STEP, ERROR_LINK.CORE_BLOCKS);
  const scriptInputs = [];
  let current;
  ivp.inits.forEach((_, key) => {
    if (key[0] === SERVICE) {
      throw new ModelError(`${ERROR_MSG.SERVICE_START} Correct **"${key}"** in the **#inits** block.`, ERROR_LINK.CORE_BLOCKS);
    }
    current = key.toLocaleLowerCase();
    if (scriptInputs.includes(current)) {
      throw new ModelError(`${ERROR_MSG.REUSE_NAME} **"${key}"** in the **#inits** block.`, ERROR_LINK.CORE_BLOCKS);
    }
    scriptInputs.push(current);
  });
  if (ivp.params !== null) {
    ivp.params.forEach((_, key) => {
      if (key[0] === SERVICE) {
        throw new ModelError(`${ERROR_MSG.SERVICE_START}: correct **"${key}"** in the **${CONTROL_EXPR.PARAMS}** block.`, ERROR_LINK.MODEL_PARAMS);
      }
      current = key.toLocaleLowerCase();
      if (scriptInputs.includes(current)) {
        throw new ModelError(`${ERROR_MSG.REUSE_NAME} **"${key}"** in the **${CONTROL_EXPR.PARAMS}** block.`, ERROR_LINK.MODEL_PARAMS);
      }
      scriptInputs.push(current);
    });
  }
  checkSolverSettings(ivp.solverSettings);
}
// Annotate the CommonJS export names for ESM import in node:
0 && (module.exports = {
  MAX_LINE_CHART,
  STAGE_COL_NAME,
  getIVP
});
