// Scripting tools for the Initial Value Problem (IVP) solver

const CONTROL_TAG = '#';
const CONTROL_SEP = ':';
const EQUAL_SIGN = '=';
const DIV_SIGN = '/';
const SERVICE = '_';
const DEFAULT_TOL = '0.00005';
export const DF_NAME = 'df';
const ODES_PACKAGE = 'ODEs';
const SOLVER_FUNC = 'solve';

type Arg = {
  name: string,
  start: number,
  finish: number,
  step: number
};

type DifEqs = {
  equations: Map<string, string>,
  solutionNames: string[]
};

/** Initial Value Problem (IVP) specification type */
type IVP = {
  name: string,
  deqs: DifEqs,
  exprs: Map<string, string> | null,
  arg: Arg,
  inits: Map<string, number>,
  consts: Map<string, number> | null,
  params: Map<string, number> | null,
  tolerance: string,
  usedMathFuncs: number[],
  usedMathConsts: number[],
};

export enum CONTROL_EXPR {
  NAME = `${CONTROL_TAG}name`,
  DIF_EQ = `${CONTROL_TAG}equations`,
  EXPR = `${CONTROL_TAG}expressions`,
  ARG = `${CONTROL_TAG}argument`,
  INITS = `${CONTROL_TAG}inits`,
  CONSTS = `${CONTROL_TAG}constants`,
  PARAMS = `${CONTROL_TAG}parameters`,
  TOL = `${CONTROL_TAG}tolerance`,
};

enum ERROR_MSG {
  ODES = 'incorrect definition of the system of ordinary differential equations',
  CONTROL_EXPR = `unsupported control expression with the tag '${CONTROL_TAG}'`,
  ARG = 'incorrect argument specification',
  INITS = 'incorrect initial values specification',
}

enum ANNOT {
  NAME = '//name:',
  TAGS = '//tags: model',
  LANG = '//language: javascript',
  INPUT = '//input: double',
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

enum SCRIPT {
  CONSTS = '\n// constants',
  ODE_COM = '\n// the problem definition',
  ODE = 'let odes = {',
  SOLVER_COM = '\n// solve the problem',
  SOLVER = `const solver = await grok.functions.eval('${ODES_PACKAGE}:${SOLVER_FUNC}');`,
  PREPARE = 'let call = solver.prepare({problem: odes});',
  CALL = 'await call.call();',
  OUTPUT = `let ${DF_NAME} = call.getParamValue('${DF_NAME}');`,
  SPACE = '    ',
  SUBSPACE = '      ',
  FUNC_VALS = '// extract function values',
  EVAL_EXPR = '// evaluate expressions',
  COMP_OUT = '// compute output',
  MATH_FUNC_COM = '\n// used Math-functions',
  MATH_CONST_COM = '\n// used Math-constants',
}

type Block = {
  begin: number,
  end: number
}

/** */
const MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh',];

/** */
const POW_IDX = 0;

/** */
const MATH_CONSTS = ['PI', 'E', ];

/**  */
function getStartOfProblemDef(lines: string[]): number {
  const linesCount = lines.length;
  let idx = 0;
 
  while (!lines[idx].startsWith(CONTROL_TAG)) {
    ++idx;
    
    if (idx === linesCount)
      throw new Error(ERROR_MSG.ODES);
  }

  let i = idx;
  
  while (i < linesCount) {
    //console.log(`line ${i}:${lines[i]}`);
    ++i;
  }

  return idx;
}

/** */
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

/** */
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

/** */
function getUsedMathIds(text: string, mathIds: string[]) {
  const res = [] as number[];
  const size = mathIds.length;  

  for (let i = 0; i < size; ++i)
    if (text.includes(`${mathIds[i]}`))
      res.push(i);

  return res;  
}

/** */
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

/**  */
function getExpressions(lines: string[]): Map<string, string> {
  const exprs = new Map<string, string>();  
  let eqIdx = 0;

  for (const line of lines) {    
    eqIdx = line.indexOf(EQUAL_SIGN);
    exprs.set(line.slice(0, eqIdx).replace(' ', ''), line.slice(eqIdx + 1).trim());
  }

  return exprs;
}

/** */
function getArg(lines: string[]): Arg {
  if (lines.length !== 4)
    throw new Error(ERROR_MSG.ARG);

  return {
    name: lines[0].slice(lines[0].indexOf(CONTROL_SEP) + 1).trim(),
    start: Number(lines[1].slice(lines[1].indexOf(EQUAL_SIGN) + 1).trim()),
    finish: Number(lines[2].slice(lines[2].indexOf(EQUAL_SIGN) + 1).trim()),
    step: Number(lines[3].slice(lines[3].indexOf(EQUAL_SIGN) + 1).trim())
  }
}

/**  */
function getEqualities(lines: string[]): Map<string, number> {
  const eqs = new Map<string, number>();  
  let eqIdx = 0;

  for (const line of lines) {    
    eqIdx = line.indexOf(EQUAL_SIGN);
    eqs.set(line.slice(0, eqIdx).replace(' ', ''), Number(line.slice(eqIdx + 1).trim()));
  }

  return eqs;
}

/** */
export function getIVP(text: string): IVP {
  // The current Initial Value Problem (IVP) specification
  let name: string;
  let deqs: DifEqs;
  let exprs: Map<string, string> | null = null;
  let arg: Arg;
  let inits: Map<string, number>;
  let consts: Map<string, number> | null = null;
  let params: Map<string, number> | null = null;
  let tolerance = DEFAULT_TOL;

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
      inits = getEqualities(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    }
    else if (firstLine.startsWith(CONTROL_EXPR.CONSTS)) { // the 'constants' block
      consts = getEqualities(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    }
    else if (firstLine.startsWith(CONTROL_EXPR.PARAMS)) { // the 'parameters' block
      params = getEqualities(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    }
    else if (firstLine.startsWith(CONTROL_EXPR.TOL)) { // the 'tolerance' block
      tolerance = firstLine.slice( firstLine.indexOf(CONTROL_SEP) + 1).trim();
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
    deqs: deqs!,
    exprs: exprs,
    arg: arg!,
    inits: inits!,
    consts: consts,
    params: params,
    tolerance: tolerance,
    usedMathFuncs: getUsedMathIds(text, MATH_FUNCS),
    usedMathConsts: getUsedMathIds(text, MATH_CONSTS),
  };
} // getIVP

/** */
function getAnnot(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  const res = [] as string[];

  // the 'name' line
  res.push(`${ANNOT.NAME} ${ivp.name}`);

  // the 'tags' line
  //res.push(ANNOT.TAGS);

  // the 'language' line
  res.push(ANNOT.LANG);

  // argument lines
  const arg = ivp.arg;
  const t0 = `${SERVICE}${arg.name}0`;
  const t1 = `${SERVICE}${arg.name}1`;
  const h = `${SERVICE}h`
  res.push(`${ANNOT.INPUT} ${t0} = ${arg.start} ${ANNOT.ARG_INIT}`);
  res.push(`${ANNOT.INPUT} ${t1} = ${arg.finish} ${ANNOT.ARG_FIN}`);
  res.push(`${ANNOT.INPUT} ${h} = ${arg.step} ${ANNOT.ARG_STEP}`);

  // initial values lines
  ivp.inits.forEach((val, key, map) => res.push(`${ANNOT.INPUT} ${key} = ${val} {${ANNOT.CAPTION} ${key}; ${ANNOT.INITS}}`));

  // parameters lines
  if (ivp.params !== null)
    ivp.params.forEach((val, key, map) => res.push(`${ANNOT.INPUT} ${key} = ${val} {${ANNOT.CAPTION} ${key}; ${ANNOT.PARAMS}}`));

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

/** */
function getMathArg(funcIdx: number): string {
  return (funcIdx > POW_IDX) ? '(x)' : '(x, y)';
}

function getScriptMainBody(ivp: IVP): string[] {  
  const res = [] as string[];

  // 1. Constants lines
  if (ivp.consts !== null) {
    res.push(SCRIPT.CONSTS);
    ivp.consts.forEach((val, key, map) => res.push(`const ${key} = ${val};`));
  }

  // 2. The problem definition lines
  res.push(SCRIPT.ODE_COM);
  res.push(SCRIPT.ODE);
  res.push(`${SCRIPT.SPACE}name: '${ivp.name}',`);

  // 2.1) argument
  const t = ivp.arg.name;
  const t0 = `${SERVICE}${t}0`;
  const t1 = `${SERVICE}${t}1`;
  const h = `${SERVICE}h`;
  res.push(`${SCRIPT.SPACE}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);

  const names = ivp.deqs.solutionNames;

  // 2.2) initial values
  res.push(`${SCRIPT.SPACE}initial: [${names.join(', ')}],`);

  // 2.3) the right-hand side of the problem
  res.push(`${SCRIPT.SPACE}func: (${t}, ${SERVICE}y, ${SERVICE}output) => {`);

  res.push(`${SCRIPT.SUBSPACE}${SCRIPT.FUNC_VALS}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SUBSPACE}const ${name} = ${SERVICE}y[${idx}];`));

  if (ivp.exprs !== null) {
    res.push(`\n${SCRIPT.SUBSPACE}${SCRIPT.EVAL_EXPR}`);
    ivp.exprs.forEach((val, key, map) => res.push(`${SCRIPT.SUBSPACE}const ${key} = ${val};`));
  }

  res.push(`\n${SCRIPT.SUBSPACE}${SCRIPT.FUNC_VALS}`);
  names.forEach((name, idx) => res.push(`${SCRIPT.SUBSPACE}${SERVICE}output[${idx}] = ${ivp.deqs.equations.get(name)};`));

  res.push(`${SCRIPT.SPACE}},`);

  // 2.4) final lines of the problem specification
  res.push(`${SCRIPT.SPACE}tolerance: ${ivp.tolerance},`);
  res.push(`${SCRIPT.SPACE}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
  res.push('};');

  // 3. Math functions
  if (ivp.usedMathFuncs.length > 0) {
    res.push(SCRIPT.MATH_FUNC_COM);
    ivp.usedMathFuncs.forEach((i) => res.push(`const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
  }

  // 4. Math constants
  if (ivp.usedMathConsts.length > 0) {
    res.push(SCRIPT.MATH_CONST_COM);
    ivp.usedMathConsts.forEach((i) => res.push(`const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
  }

  // 5. The 'call solver' lines
  res.push(SCRIPT.SOLVER_COM);
  res.push(SCRIPT.SOLVER);
  res.push(SCRIPT.PREPARE);
  res.push(SCRIPT.CALL);
  res.push(SCRIPT.OUTPUT);

  return res;
} // getScriptMainBody

/** */
export function getScriptLines(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  return getAnnot(ivp, toAddViewers, toAddEditor).concat(getScriptMainBody(ivp));
}

/** */
export function getScriptParams(ivp: IVP): Record<string, number> {
  const res = {} as Record<string, number>;

  const arg = ivp.arg;

  res[`${SERVICE}${arg.name}0`] = arg.start;
  res[`${SERVICE}${arg.name}1`] = arg.finish;
  res[`${SERVICE}h`] = arg.step;

  ivp.inits.forEach((val, key) => res[key] = val);

  if (ivp.params)
    ivp.params.forEach((val, key) => res[key] = val);

  return res;
}


/*
const ivp = getIVP(simple);
console.log(ivp);

console.log();

const scriptLines = getScriptLines(ivp);

for (const line of scriptLines)
  console.log(line);*/
