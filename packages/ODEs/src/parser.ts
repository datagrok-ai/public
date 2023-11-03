const eqs = `
Some very-y-y-y smart text
Another line
()))((()))

  #name: Complex Test
#differential equations:
  d(x)/dt = coef1 * y

  dy/d(t) = coef2 
   * x

               #expressions:
  coef1 = const1 + param1



  coef2 = const2 
   + param2
   + 0.0

#argument: t
  start = 
    0.0


  finish = 5.0



  step = 0.1

#initial values:
  x = 2.0
  y = 0.0

#constants:
  const1 = 1.0
  const2 = 3.0

#parameters:
  param1 = 1.0
  param2 = -1.0

#tolerance: 0.00005`;

const CONTROL_TAG = '#';
const CONTROL_SEP = ':';
const EQUAL_SIGN = '=';
const DIV_SIGN = '/';
const SERVICE = '_';
const DEFAULT_TOL = 0.00005;

type Arg = {
  name: string,
  start: number,
  finish: number,
  step: number
};

/** Initial Value Problem (IVP) specification type */
type IVP = {
  name: string,
  deqs: Map<string, string>,
  exprs: Map<string, string> | null,
  arg: Arg,
  inits: Map<string, number>,
  consts: Map<string, number> | null,
  params: Map<string, number> | null,
  tolerance: number,
};

enum CONTROL_EXPR {
  NAME = `${CONTROL_TAG}name`,
  DIF_EQ = `${CONTROL_TAG}differential equations`,
  EXPR = `${CONTROL_TAG}expressions`,
  ARG = `${CONTROL_TAG}argument`,
  INITS = `${CONTROL_TAG}initial values`,
  CONSTS = `${CONTROL_TAG}constants`,
  PARAMS = `${CONTROL_TAG}parameters`,
  TOL = `${CONTROL_TAG}tolerance`,
}

enum ERROR_MSG {
  ODES = 'incorrect definition of the system of ordinary differential equations',
  CONTROL_EXPR = `unsupported control expression with the tag '${CONTROL_TAG}'`,
  ARG = 'incorrect argument specification',
  INITS = 'incorrect initial values specification',
}

enum ANNOT {
  NAME = '//name:',
  LANG = '//language: javascript',
  INPUT = '//input: double',
  OUTPUT = '//input: dataframe df',
  EDITOR = '//editor: Compute:RichFunctionViewEditor',
  CAPTION = 'caption:',
  ARG_INIT = '{caption: Initial; category: Argument}',
  ARG_FIN = '{caption: Final; category: Argument}',
  ARG_STEP = '{caption: Step; category: Argument}',
  INITS = 'category: Initial values',
  PARAMS = 'category: Parameters',
  VIEWERS = 'viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)'
}

type Block = {
  begin: number,
  end: number
}

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
    console.log(`line ${i}:${lines[i]}`);
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
function getDifEquations(lines: string[]): Map<string, string> {
  const deqs = new Map<string, string>();

  let divIdx = 0;
  let eqIdx = 0;

  for (const line of lines) {
    divIdx = line.indexOf(DIV_SIGN);
    eqIdx = line.indexOf(EQUAL_SIGN);
    deqs.set(line.slice(0, divIdx).replace(/[d() ]/g, ''), line.slice(eqIdx + 1).trim());
  }

  return deqs;
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
function getIVP(text: string): IVP {
  // The current Initial Value Problem (IVP) specification
  let name: string;
  let deqs: Map<string, string>;
  let exprs: Map<string, string> | null = null;
  let arg: Arg;
  let inits: Map<string, number>;
  let consts: Map<string, number> | null = null;
  let params: Map<string, number> | null = null;
  let tolerance = DEFAULT_TOL;

  // 0. Split text into lines
  const lines = eqs.split('\n').filter((s) => s !== '').map((s) => s.trimStart());

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
      tolerance = Number(firstLine.slice( firstLine.indexOf(CONTROL_SEP) + 1).trim());
    }
    else // error: unsupported control expression 
    {
      console.log(firstLine);
      throw new Error(ERROR_MSG.CONTROL_EXPR);
    }
  }

  if (inits!.size !== deqs!.size)
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
  };
}

/** */
function getAnnot(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  const res = [] as string[];

  // the 'name' line
  res.push(`${ANNOT.NAME} ${ivp.name}`);

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
}

const ivp = getIVP(eqs);
console.log(ivp);

console.log();

const annot = getAnnot(ivp);

for (const line of annot)
  console.log(line);
