// Scripting tools for the Initial Value Problem (IVP) solver

import {CONTROL_TAG, DF_NAME, CONTROL_EXPR, DOSING} from './constants';

// Scripting specific constants
const CONTROL_SEP = ':';
const EQUAL_SIGN = '=';
const DIV_SIGN = '/';
const SERVICE = '_';
const BRACE = '{';
const DEFAULT_TOL = '0.00005';

/** Solver package name */
const ODES_PACKAGE = 'ODEs';

/** Numerical solver function */
const SOLVER_FUNC = 'solve';

/** Elementary math tools */
const MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh',];
const POW_IDX = MATH_FUNCS.indexOf('pow');
const MATH_CONSTS = ['PI', 'E', ];

/** Numerical input specification */
type Input = {
  value: number,
  annot: string | null,
};

/** Argument of IVP specification */
type Arg = {
  name: string,
  start: Input,
  finish: Input,
  step: Input
};

/** Differential equations specification */
type DifEqs = {
  equations: Map<string, string>,
  solutionNames: string[]
};

/** Dosing specification */
type Dosing = {
  dose: Input,
  count: Input,
  updates: string[],
};

/** Initial Value Problem (IVP) specification type */
type IVP = {
  name: string,
  tags: string | null,
  descr: string | null,
  deqs: DifEqs,
  exprs: Map<string, string> | null,
  arg: Arg,
  inits: Map<string, Input>,
  consts: Map<string, Input> | null,
  params: Map<string, Input> | null,
  tolerance: string,
  usedMathFuncs: number[],
  usedMathConsts: number[],
  dosing: Dosing | null,
};

/** Specific error messeges */
enum ERROR_MSG {
  ODES = 'incorrect definition of the system of ordinary differential equations',
  CONTROL_EXPR = `unsupported control expression with the tag '${CONTROL_TAG}'`,
  ARG = 'incorrect argument specification',
  INITS = 'incorrect initial values specification',
  DOSING = 'incorrect dosing specification',
  DOSE = 'incorrect dose',
  DOSES_COUNT = 'incorrect doses count',
}

/** Datagrok annatations */
enum ANNOT {
  NAME = '//name:',
  DESCR = '//description:',
  TAGS = '//tags:',
  LANG = '//language: javascript',
  DOUBLE_INPUT = '//input: double',
  INT_INPUT = '//input: int',
  OUTPUT = `//output: dataframe ${DF_NAME}`,
  EDITOR = '//editor: Compute:RichFunctionViewEditor',
  CAPTION = 'caption:',
  ARG_INIT = '{caption: Initial; category: Argument}',
  ARG_FIN = '{caption: Final; category: Argument}',
  ARG_STEP = '{caption: Step; category: Argument}',
  INITS = 'category: Initial values',
  PARAMS = 'category: Parameters',
  VIEWERS = 'viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)'
}

/** JS-scripting components */
enum SCRIPT {
  CONSTS = '// constants',
  ODE_COM = '// the problem definition',
  ODE = 'let odes = {',
  SOLVER_COM = '// solve the problem',
  SOLVER = `const solver = await grok.functions.eval('${ODES_PACKAGE}:${SOLVER_FUNC}');`,
  PREPARE = 'let call = solver.prepare({problem: odes});',
  CALL = 'await call.call();',
  OUTPUT = `let ${DF_NAME} = call.getParamValue('${DF_NAME}');`,
  SPACE2 = '  ',
  SPACE4 = '    ',
  SPACE6 = '      ',
  SPACE8 = '        ',
  FUNC_VALS = '// extract function values',
  EVAL_EXPR = '// evaluate expressions',
  COMP_OUT = '// compute output',
  MATH_FUNC_COM = '// used Math-functions',
  MATH_CONST_COM = '// used Math-constants',
  ONE_STAGE_COM = '\n// one stage solution',
  ONE_STAGE_BEGIN = 'let _oneStage = async (',
  ONE_STAGE_END = ') => {',
  ASYNC_OUTPUT = `let ${DF_NAME} = await _oneStage(`,
  RETURN_OUTPUT = `return call.getParamValue('${DF_NAME}');`,
  EMPTY_OUTPUT = `let ${DF_NAME} = DG.DataFrame.create();`,
  APPEND_ASYNC = `${DF_NAME}.append(await _oneStage(`,
  SOLUTION_DF_COM = '// solution dataframe',
  DOSE_INTERVAL_COM = '// dosing interval',
  DOSE_INTERVAL = `${SERVICE}interval`,
  LAST_IDX = `${SERVICE}lastIdx`,
}

/** Limits of the problem specification */
type Block = {
  begin: number,
  end: number
}

/** Get strat of the problem skipping note-lines */
function getStartOfProblemDef(lines: string[]): number {
  const linesCount = lines.length;
  let idx = 0;
 
  while (!lines[idx].startsWith(CONTROL_TAG)) {
    ++idx;
    
    if (idx === linesCount)
      throw new Error(ERROR_MSG.ODES);
  }

  return idx;
}

/** Get limits of blocks specifying IVP */
function getBlocks(lines: string[], start: number): Block[] {
  const linesCount = lines.length;
  let beg = start;
  let idx = start + 1;
  const blocks = [] as Block[];

  while (true) {
    if (idx === linesCount) {
      blocks.push({begin: beg, end: idx});
      break;
    }

    if (lines[idx].startsWith(CONTROL_TAG)) {
      blocks.push({begin: beg, end: idx});
      beg = idx;
    }

    ++idx;
  }

  return blocks;
}

/** Concatenate multi-line formulas & expressions */
function concatMultilineFormulas(source: string[]): string[] {
  const res: string[] = [ source[0] ];

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

/** Get indeces of math functions used in the current text */
function getUsedMathIds(text: string, mathIds: string[]) {
  const res = [] as number[];
  const size = mathIds.length;  

  for (let i = 0; i < size; ++i)
    if (text.includes(`${mathIds[i]}`))
      res.push(i);

  return res;  
}

/** Get differential equations */
function getDifEquations(lines: string[]): DifEqs {
  const deqs = new Map<string, string>();
  const names = [] as string[];

  let divIdx = 0;
  let eqIdx = 0;

  for (const line of lines) {
    divIdx = line.indexOf(DIV_SIGN);
    eqIdx = line.indexOf(EQUAL_SIGN);
    const name = line.slice(line.indexOf('d') + 1, divIdx).replace(/[() ]/g, '');
    names.push(name);
    deqs.set(name, line.slice(eqIdx + 1).trim());
  }

  return {
    equations: deqs,
    solutionNames: names
  };
}

/** Get expressions of IVP */
function getExpressions(lines: string[]): Map<string, string> {
  const exprs = new Map<string, string>();  
  let eqIdx = 0;

  for (const line of lines) {    
    eqIdx = line.indexOf(EQUAL_SIGN);
    exprs.set(line.slice(0, eqIdx).replace(' ', ''), line.slice(eqIdx + 1).trim());
  }

  return exprs;
}

/** Get input specification */
function getInput(line: string) : Input {
  const str = line.slice(line.indexOf(EQUAL_SIGN) + 1).trim();
  const braceIdx = str.indexOf(BRACE);

  if (braceIdx === -1)
    return {
      value: Number(str),
      annot: null,
    };
  
  return {
    value: Number(str.slice(0, braceIdx)),
    annot: str.slice(braceIdx),
  }
}

/** Get argument (independent variable) of IVP */
function getArg(lines: string[]): Arg {
  if (lines.length !== 4)
    throw new Error(ERROR_MSG.ARG);

  return {
    name: lines[0].slice(lines[0].indexOf(CONTROL_SEP) + 1).trim(),
    start: getInput(lines[1]),
    finish: getInput(lines[2]),
    step: getInput(lines[3]),
  }
}

/** Get equlities specification */
function getEqualities(lines: string[], begin: number, end: number): Map<string, Input> {
  const source = concatMultilineFormulas(lines.slice(begin, end));
  const eqs = new Map<string, Input>();
  let eqIdx = 0;

  for (const line of source) {    
    eqIdx = line.indexOf(EQUAL_SIGN);
    eqs.set(line.slice(0, eqIdx).replace(' ', ''), getInput(line.slice(eqIdx + 1).trim()));
  }

  return eqs;
}

/** Get dosing specification */
function getDosing(lines: string[], begin: number, end: number): Dosing {
  const source = concatMultilineFormulas(lines.slice(begin, end));
  const size = source.length;

  if (size < DOSING.MIN_LINES_COUNT)
    throw new Error(ERROR_MSG.DOSING);

  const dose = getInput(source[DOSING.DOSE_IDX]);
  if (dose.value < DOSING.MIN_DOSE)
    throw new Error(ERROR_MSG.DOSE);
  
  const count = getInput(source[DOSING.COUNT_IDX]);
  if (count.value <  DOSING.MIN_DOSES_COUNT)
    throw new Error(ERROR_MSG.DOSES_COUNT);

  return {dose: dose, count: count, updates: source.slice(DOSING.COUNT_IDX + 1)};
}

/** Get initial value problem specification given in the text */
export function getIVP(text: string): IVP {
  // The current Initial Value Problem (IVP) specification
  let name: string;
  let tags: string | null = null;
  let descr: string | null = null;
  let deqs: DifEqs;
  let exprs: Map<string, string> | null = null;
  let arg: Arg;
  let inits: Map<string, Input>;
  let consts: Map<string, Input> | null = null;
  let params: Map<string, Input> | null = null;
  let tolerance = DEFAULT_TOL;
  let dosing: Dosing | null = null;

  // 0. Split text into lines
  const lines = text.split('\n').filter((s) => s !== '').map((s) => s.trimStart());

  // 1. Skip first lines without the control tag  
  const start = getStartOfProblemDef(lines);

  // 2. Get blocks limits
  const blocks = getBlocks(lines, start);

  // 3. Process blocks
  for (const block of blocks) {
    const firstLine = lines[block.begin];

    if (firstLine.startsWith(CONTROL_EXPR.NAME)) { // the 'name' block      
      name = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    }
    else if (firstLine.startsWith(CONTROL_EXPR.TAGS)) { // the 'tags' block      
      tags = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    }
    else if (firstLine.startsWith(CONTROL_EXPR.DESCR)) { // the 'description' block      
      descr = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
    }
    else if (firstLine.startsWith(CONTROL_EXPR.DIF_EQ)) { // the 'differential equations' block      
      deqs = getDifEquations(concatMultilineFormulas(lines.slice(block.begin + 1, block.end))); 
    }
    else if (firstLine.startsWith(CONTROL_EXPR.EXPR)) { // the 'expressions' block
      exprs = getExpressions(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    }
    else if (firstLine.startsWith(CONTROL_EXPR.ARG)) { // the 'argument' block
      arg = getArg(concatMultilineFormulas(lines.slice(block.begin, block.end)));
    }
    else if (firstLine.startsWith(CONTROL_EXPR.INITS)) { // the 'initial values' block
      inits = getEqualities(lines, block.begin + 1, block.end);
    }
    else if (firstLine.startsWith(CONTROL_EXPR.CONSTS)) { // the 'constants' block
      consts = getEqualities(lines, block.begin + 1, block.end);
    }
    else if (firstLine.startsWith(CONTROL_EXPR.PARAMS)) { // the 'parameters' block
      params = getEqualities(lines, block.begin + 1, block.end);
    }
    else if (firstLine.startsWith(CONTROL_EXPR.TOL)) { // the 'tolerance' block
      tolerance = firstLine.slice( firstLine.indexOf(CONTROL_SEP) + 1).trim();
    }
    else if (firstLine.startsWith(CONTROL_EXPR.DOSING)) { // the 'dosing' block
      dosing = getDosing(lines, block.begin + 1, block.end);
    }
    else // error: unsupported control expression 
    {
      //console.log(firstLine);
      throw new Error(ERROR_MSG.CONTROL_EXPR);
    }
  }

  if (inits!.size !== deqs!.solutionNames.length)
    throw new Error(ERROR_MSG.INITS);

  return {
    name: name!,
    tags: tags,
    descr: descr,
    deqs: deqs!,
    exprs: exprs,
    arg: arg!,
    inits: inits!,
    consts: consts,
    params: params,
    tolerance: tolerance,
    usedMathFuncs: getUsedMathIds(text, MATH_FUNCS),
    usedMathConsts: getUsedMathIds(text, MATH_CONSTS),
    dosing: dosing
  };
} // getIVP

/** Return input specification, required for annotation generating */
function getInputSpec(inp: Input): string {
  if (inp.annot)
    return `${inp.value} ${inp.annot}`;

  return `${inp.value}`;
}

/** Generate annotation lines */
function getAnnot(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  const res = [] as string[];

  // the 'name' line
  res.push(`${ANNOT.NAME} ${ivp.name}`);

  // the 'tags' line
  if (ivp.tags)
    res.push(`${ANNOT.TAGS} ${ivp.tags}`);

  // the 'description' line
  if (ivp.descr)
    res.push(`${ANNOT.DESCR} ${ivp.descr}`); 

  // the 'language' line
  res.push(ANNOT.LANG);

  // the 'dosing' lines
  if (ivp.dosing) {
    res.push(`${ANNOT.DOUBLE_INPUT} ${DOSING.DOSE} = ${getInputSpec(ivp.dosing.dose)}`);
    res.push(`${ANNOT.INT_INPUT} ${DOSING.COUNT} = ${getInputSpec(ivp.dosing.count)}`);
  }

  // argument lines
  const arg = ivp.arg;
  const t0 = `${SERVICE}${arg.name}0`;
  const t1 = `${SERVICE}${arg.name}1`;
  const h = `${SERVICE}h`
  res.push(`${ANNOT.DOUBLE_INPUT} ${t0} = ${getInputSpec(arg.start)}`);
  res.push(`${ANNOT.DOUBLE_INPUT} ${t1} = ${getInputSpec(arg.finish)}`);
  res.push(`${ANNOT.DOUBLE_INPUT} ${h} = ${getInputSpec(arg.step)}`);

  // initial values lines
  ivp.inits.forEach((val, key) => res.push(`${ANNOT.DOUBLE_INPUT} ${key} = ${getInputSpec(val)}`));

  // parameters lines
  if (ivp.params !== null)
    ivp.params.forEach((val, key) => res.push(`${ANNOT.DOUBLE_INPUT} ${key} = ${getInputSpec(val)}`));

  // the 'output' line
  if (toAddViewers)
    res.push(`${ANNOT.OUTPUT} {${ANNOT.CAPTION} ${ivp.name}; ${ANNOT.VIEWERS}}`);
  else
    res.push(`${ANNOT.OUTPUT} {${ANNOT.CAPTION} ${ivp.name}`);

  // the 'editor' line
  if (toAddEditor)
    res.push(ANNOT.EDITOR);

  return res;
} // getAnnot

/** Returns math functions arguments string */
function getMathArg(funcIdx: number): string {
  return (funcIdx > POW_IDX) ? '(x)' : '(x, y)';
}

/** Return main body of JS-script: basic variant */
function getScriptMainBodyBasic(ivp: IVP): string[] {  
  const res = [] as string[];

  // 1. Constants lines
  if (ivp.consts !== null) {
    res.push('');
    res.push(SCRIPT.CONSTS);
    ivp.consts.forEach((val, key) => res.push(`const ${key} = ${val.value};`));
  }

  // 2. The problem definition lines
  res.push('');
  res.push(SCRIPT.ODE_COM);
  res.push(SCRIPT.ODE);
  res.push(`${SCRIPT.SPACE4}name: '${ivp.name}',`);

  // 2.1) argument
  const t = ivp.arg.name;
  const t0 = `${SERVICE}${t}0`;
  const t1 = `${SERVICE}${t}1`;
  const h = `${SERVICE}h`;
  res.push(`${SCRIPT.SPACE4}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);

  const names = ivp.deqs.solutionNames;

  // 2.2) initial values
  res.push(`${SCRIPT.SPACE4}initial: [${names.join(', ')}],`);

  // 2.3) the right-hand side of the problem
  res.push(`${SCRIPT.SPACE4}func: (${t}, ${SERVICE}y, ${SERVICE}output) => {`);

  res.push(`${SCRIPT.SPACE6}${SCRIPT.FUNC_VALS}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SPACE6}const ${name} = ${SERVICE}y[${idx}];`));

  if (ivp.exprs !== null) {
    res.push(`\n${SCRIPT.SPACE6}${SCRIPT.EVAL_EXPR}`);
    ivp.exprs.forEach((val, key, map) => res.push(`${SCRIPT.SPACE6}const ${key} = ${val};`));
  }

  res.push(`\n${SCRIPT.SPACE6}${SCRIPT.COMP_OUT}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SPACE6}${SERVICE}output[${idx}] = ${ivp.deqs.equations.get(name)};`));

  res.push(`${SCRIPT.SPACE4}},`);

  // 2.4) final lines of the problem specification
  res.push(`${SCRIPT.SPACE4}tolerance: ${ivp.tolerance},`);
  res.push(`${SCRIPT.SPACE4}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
  res.push('};');

  // 3. Math functions
  if (ivp.usedMathFuncs.length > 0) {
    res.push('');
    res.push(SCRIPT.MATH_FUNC_COM);
    ivp.usedMathFuncs.forEach((i) => res.push(`const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
  }

  // 4. Math constants
  if (ivp.usedMathConsts.length > 0) {
    res.push('');
    res.push(SCRIPT.MATH_CONST_COM);
    ivp.usedMathConsts.forEach((i) => res.push(`const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
  }

  // 5. The 'call solver' lines
  res.push('');
  res.push(SCRIPT.SOLVER_COM);
  res.push(SCRIPT.SOLVER);
  res.push(SCRIPT.PREPARE);
  res.push(SCRIPT.CALL);
  res.push(SCRIPT.OUTPUT);

  return res;
} // getScriptMainBodyBasic

/** Return function for JS-script */
function getScriptFunc(ivp: IVP, funcParamsNames: string): string[] {  
  const res = [] as string[];

  // 0. Function declaration
  res.push(SCRIPT.ONE_STAGE_COM);
  res.push(`${SCRIPT.ONE_STAGE_BEGIN}${funcParamsNames}${SCRIPT.ONE_STAGE_END}`);

  // 1. Constants lines
  if (ivp.consts !== null) {
    res.push('');
    res.push(SCRIPT.SPACE2 + SCRIPT.CONSTS);
    ivp.consts.forEach((val, key) => res.push(`${SCRIPT.SPACE2}const ${key} = ${val.value};`));
  }

  // 2. The problem definition lines
  res.push(SCRIPT.SPACE2 + SCRIPT.ODE_COM);
  res.push(SCRIPT.SPACE2 + SCRIPT.ODE);
  res.push(`${SCRIPT.SPACE6}name: '${ivp.name}',`);

  // 2.1) argument
  const t = ivp.arg.name;
  const t0 = `${SERVICE}${t}0`;
  const t1 = `${SERVICE}${t}1`;
  const h = `${SERVICE}h`;
  res.push(`${SCRIPT.SPACE6}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);

  const names = ivp.deqs.solutionNames;

  // 2.2) initial values
  res.push(`${SCRIPT.SPACE6}initial: [${names.join(', ')}],`);

  // 2.3) the right-hand side of the problem
  res.push(`${SCRIPT.SPACE6}func: (${t}, ${SERVICE}y, ${SERVICE}output) => {`);

  res.push(`${SCRIPT.SPACE8}${SCRIPT.FUNC_VALS}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SPACE8}const ${name} = ${SERVICE}y[${idx}];`));

  if (ivp.exprs !== null) {
    res.push(`\n${SCRIPT.SPACE8}${SCRIPT.EVAL_EXPR}`);
    ivp.exprs.forEach((val, key, map) => res.push(`${SCRIPT.SPACE8}const ${key} = ${val};`));
  }

  res.push(`\n${SCRIPT.SPACE8}${SCRIPT.COMP_OUT}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SPACE8}${SERVICE}output[${idx}] = ${ivp.deqs.equations.get(name)};`));

  res.push(`${SCRIPT.SPACE6}},`);

  // 2.4) final lines of the problem specification
  res.push(`${SCRIPT.SPACE6}tolerance: ${ivp.tolerance},`);
  res.push(`${SCRIPT.SPACE6}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
  res.push(`${SCRIPT.SPACE2}};`);

  // 3. Math functions
  if (ivp.usedMathFuncs.length > 0) {
    res.push('');
    res.push(SCRIPT.SPACE2 + SCRIPT.MATH_FUNC_COM);
    ivp.usedMathFuncs.forEach((i) => res.push(`${SCRIPT.SPACE2}const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
  }
  
  // 4. Math constants
  if (ivp.usedMathConsts.length > 0) {
    res.push('');
    res.push(SCRIPT.SPACE2 + SCRIPT.MATH_CONST_COM);
    ivp.usedMathConsts.forEach((i) => res.push(`${SCRIPT.SPACE2}const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
  }  

  // 5. The 'call solver' lines
  res.push('');
  res.push(SCRIPT.SPACE2 + SCRIPT.SOLVER_COM);
  res.push(SCRIPT.SPACE2 + SCRIPT.SOLVER);
  res.push(SCRIPT.SPACE2 + SCRIPT.PREPARE);
  res.push(SCRIPT.SPACE2 + SCRIPT.CALL);
  res.push(SCRIPT.SPACE2 + SCRIPT.RETURN_OUTPUT);

  // 6. Close the function
  res.push('};');

  return res;
} // getScriptFunc

/** Return main body of JS-script: dosing case */
function getScriptMainBodyDosingCase(ivp: IVP): string[] {
  const funcParamsNames = getFuncParamsNames(ivp);
  const res = getScriptFunc(ivp, funcParamsNames);
  res.push('');
  //res.push(`${SCRIPT.ASYNC_OUTPUT}${funcParamsNames});`);

  res.push(SCRIPT.SOLUTION_DF_COM);
  const dfNames = getSolutionDfColsNames(ivp);

  res.push(`let ${DF_NAME} = DG.DataFrame.fromColumns([`);
  dfNames.forEach((name) => res.push(`${SCRIPT.SPACE2}DG.Column.fromFloat32Array('${name}', []),`));
  res.push(`]);`);
  res.push(`${DF_NAME}.name = '${ivp.name}';`);
  res.push('');

  res.push(SCRIPT.DOSE_INTERVAL_COM);
  res.push(`const ${SCRIPT.DOSE_INTERVAL} = ${SERVICE}${ivp.arg.name}1 - ${SERVICE}${ivp.arg.name}0;`);
  res.push('');
  res.push(`let ${SCRIPT.LAST_IDX} = 0`);

  res.push(SCRIPT.SOLVER_COM);
  res.push(`for (let ${SERVICE}idx = 0; ${SERVICE}idx < ${DOSING.COUNT}; ++${SERVICE}idx) {`);
  ivp.dosing!.updates.forEach((upd) => res.push(`${SCRIPT.SPACE2}${upd};`));
  res.push(`${SCRIPT.SPACE2}${SCRIPT.APPEND_ASYNC}${funcParamsNames}), true);`);
  res.push(`${SCRIPT.SPACE2}${SERVICE}${ivp.arg.name}0 = ${SERVICE}${ivp.arg.name}1;`);
  res.push(`${SCRIPT.SPACE2}${SERVICE}${ivp.arg.name}1 += ${SCRIPT.DOSE_INTERVAL};`);
  res.push(`${SCRIPT.SPACE2}${SCRIPT.LAST_IDX} = ${DF_NAME}.rowCount - 1;`);

  dfNames.forEach((name, idx) => { 
    if (idx !== 0)
      res.push(`${SCRIPT.SPACE2}${name} = ${DF_NAME}.get('${name}', ${SCRIPT.LAST_IDX});`);
  });

  res.push('};');

  return res;
} // getScriptMainBodyDosingCase

/** Return main body of JS-script */
function getScriptMainBody(ivp: IVP): string[] {
  if (ivp.dosing)
    return  getScriptMainBodyDosingCase(ivp);
  
  return getScriptMainBodyBasic(ivp);
}

/** Return JS-script lines */
export function getScriptLines(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  return getAnnot(ivp, toAddViewers, toAddEditor).concat(getScriptMainBody(ivp));
}

/** Return parameters of JS-script */
export function getScriptParams(ivp: IVP): Record<string, number> {
  const res = {} as Record<string, number>;

  if (ivp.dosing) {
    res[DOSING.DOSE] = ivp.dosing.dose.value;
    res[DOSING.COUNT] = ivp.dosing.count.value;
  }

  const arg = ivp.arg;

  res[`${SERVICE}${arg.name}0`] = arg.start.value;
  res[`${SERVICE}${arg.name}1`] = arg.finish.value;
  res[`${SERVICE}h`] = arg.step.value;

  ivp.inits.forEach((val, key) => res[key] = val.value);

  if (ivp.params)
    ivp.params.forEach((val, key) => res[key] = val.value);

  return res;
}

/** Return func parameters names string */
function getFuncParamsNames(ivp: IVP): string {
  const names = [] as string [];

  const arg = ivp.arg.name;

  names.push(`${SERVICE}${arg}0`);
  names.push(`${SERVICE}${arg}1`);
  names.push(`${SERVICE}h`);  

  ivp.inits.forEach((val, key) => names.push(key));

  if (ivp.params)
    ivp.params.forEach((val, key) => names.push(key));

  return names.join(', ');
}

/** Return solution dataframe columns names*/
function getSolutionDfColsNames(ivp: IVP): string[] {
  const res = [] as string[];

  res.push(ivp.arg.name);

  ivp.inits.forEach((val, key) => res.push(key));

  return res;
}


// TODO: to remove the following debugging lines:
/*const TEMPLATE_BASIC = `${CONTROL_EXPR.NAME}: Template 
${CONTROL_EXPR.DIF_EQ}:
  dy/dt = -y + sin(t) / t

${CONTROL_EXPR.ARG}: t
  initial = 0.01 {units: sec; caption: Initial; category: Time} [Initial value of t]
  final = 15.0
  step = 0.001{units: sec; caption: Step; category: Time} [Step of numerical solution]

${CONTROL_EXPR.INITS}:  
  y = 0`;


const TEMPLATE_ADVANCED = `NOTES. This is an advanced template. Modify it. 
Use multi-line formulas if needed.
Add new equations, expressions, constants & parameters.
Edit these header lines if required.

${CONTROL_EXPR.NAME}: Advanced
${CONTROL_EXPR.DESCR}: 2D ordinary differential equations system sample
${CONTROL_EXPR.DIF_EQ}:
  dx/dt = E1 * y + sin(t)

  dy/dt = E2 * x - pow(t, 5)

${CONTROL_EXPR.EXPR}:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2

${CONTROL_EXPR.ARG}: t
  start = 0.0
  finish = 2.0
  step = 0.01

${CONTROL_EXPR.INITS}:
  x = 2.0
  y = 0.0

${CONTROL_EXPR.CONSTS}:
  C1 = 1.0
  C2 = 3.0

${CONTROL_EXPR.PARAMS}:
  P1 = 1.0
  P2 = -1.0

${CONTROL_EXPR.TOL}: 0.00005`;*/

/*import {TEMPLATES} from './constants';

const lines = TEMPLATES.PK_PD;

console.log(lines);
console.log('==============================================================================================');

const ivp = getIVP(lines);
console.log(ivp);

console.log('==============================================================================================');
console.log(getScriptLines(ivp));*/
