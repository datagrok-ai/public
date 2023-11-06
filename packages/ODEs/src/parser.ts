const complex = `
Some very-y-y-y smart text
Another line
()))((()))

  #name: Complex Test
#differential equations:
  d(x)/dt = coef1 * y + tanh(t)

  dy/d(t) = coef2 
   * x - pow(t, 3)

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

const simple = `
#name: Simple
#differential equations:
  dy/d(t) = mu * y + sin(PI * t) + E + smth

#expressions:
  smth = exp(t) + sqrt(t) + log(t) + pow(t, 0.4)
  
#argument: t
  start = 0.1
  finish = 5.0
  step = 0.1

#initial values:  
  y = 1

#parameters:
  mu = -1`;

const MM = `
#name: MarkovModel
#differential equations:
  dp0/dt = lambda00 * p0 + lambda10 * p1
  dp1/dt = lambda01 * p0 + lambda11 * p1 + lambda21 * p2

#expressions:
  p2 = 1 - p0 - p1
  lambda00 = -lambda
  lambda10 = mu
  lambda01 = lambda
  lambda11 = -lambda - mu
  lambda21 = mu

#argument: t 
  initial = 0.0
  final = 2.0
  step = 0.01

#initial values:
  p0 = 1.0
  p1 = 0.0

#parameters:
  lambda = 1.0
  mu = 1.0`;

const stiff = `
#name: stiff
#differential equations:
  dx/dt = -5.0 * x + 3.0 * y
  dy/dt = 100.0 * x - 301.0 * y

#argument: t
  initial = 0.0
  final = 2.0
  step = 0.1

#initial values:
  x = 52.29
  y = 83.82`;

const bio = `
#name: Bioreactor
#differential equations:

  d(FFox)/dt = -FFox_to_FFred + FFred_to_FFox

  d(KKox)/dt = -KKox_to_KKred + KKred_to_KKox

  d(FFred)/dt = FFox_to_FFred - FFred_to_FFox - FFred_to_Ffree + Ffree_to_FFred

  d(KKred)/dt = KKox_to_KKred - KKred_to_KKox - KKred_to_Kfree + Kfree_to_KKred

  d(Ffree)/dt = 2.0 * FFred_to_Ffree - 2.0 * Ffree_to_FFred - free_to_FKred + FKred_to_free

  d(Kfree)/dt = 2.0 * KKred_to_Kfree - 2.0 * Kfree_to_KKred - free_to_FKred + FKred_to_free

  d(FKred)/dt = free_to_FKred - FKred_to_free - FKred_to_FKox + FKox_to_FKred

  d(FKox)/dt = FKred_to_FKox - FKox_to_FKred

  d(MEAthiol)/dt = 2.0 * (-FFox_to_FFred + FFred_to_FFox - KKox_to_KKred
			+ KKred_to_KKox + FFred_to_Ffree + KKred_to_Kfree - Ffree_to_FFred
			- Kfree_to_KKred - FKox_to_FKred
			- kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2)
			- (MEAthiol + MEAthiolate) * (Fin + Fpermeate) / VL

  d(CO2)/dt = (Fin * CO2in - 2.0 * Fpermeate * CO2) / VL + OTR
	- 0.5 * kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2
			
  d(yO2P)/dt = -OTR * (VL / Vg) * R * T * P + yO2in * qin - yO2P * qout

  d(CYST)/dt = kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2
		- krcyst * CYST * Vres - (Fin + Fpermeate) * CYST / VL

  d(VL)/dt = Fin - Fpermeate

#expressions:

  constForklasurface = (3.932 * pow((pow(AgitatorSpeed, 3.0) * pow(AgitatorDiameter, 5.0) 
   * AgitatorPowerNumber / 2160000000000), 0.361)) / 60.0

  klasurface = pow(VL, -0.65) * constForklasurface

  MEAthiolate = MEAthiol * pow(10.0,(pH - pKa2MEA))

  qout = qin - klasurface*(yO2P*H - CO2) * VL * R * T / (P*1000.0)
  
  OTR = klasurface*(yO2P*H - CO2)

  Vg = Vtotalvessel - VL

  Fin = t < TimeToSwitch ? (0.0) : (0.025)

  Fpermeate = t < TimeToSwitch ? (0.025) : (Fin)
  
  CO2in = percentO2saturation * 7.17 / (32.0 * 100.0)

  Vres = VLinitial / VL

  MEAthiolate_t_by_Vres_t_squared = pow(MEAthiolate * Vres, 2.0)  
  
  FFox_to_FFred = k1red * FFox * Vres * MEAthiolate_t_by_Vres_t_squared

  FFred_to_FFox = k1ox * FFred * Vres

  FFred_to_Ffree = k2Fd * FFred * Vres

  Ffree_to_FFred = k2Fa * pow(Ffree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared

  KKox_to_KKred = k1red * KKox * Vres * MEAthiolate_t_by_Vres_t_squared

  KKred_to_KKox = k1ox * KKred * Vres

  KKred_to_Kfree = k2Kd * KKred * Vres

  Kfree_to_KKred = k2Ka * pow(Kfree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared

  free_to_FKred = k3FKa * Ffree * Vres * Kfree * Vres

  FKred_to_free = k3FKd * FKred * Vres

  FKred_to_FKox = k4ox * FKred * Vres * pow(CYST * Vres, 2.0)

  FKox_to_FKred = k4red * FKox * Vres * MEAthiolate_t_by_Vres_t_squared

  Vres_t_by_CO2 = Vres * CO2

  sqrt_of_Vres_t_by_CO2 = (Vres_t_by_CO2 >= 0.0) ? sqrt(Vres_t_by_CO2) : 0.0		

  MEAthiol_t_by_Vres_t_squared = pow(MEAthiol * Vres, 2.0)

#argument: t
  initial = 0.0
  final = 1000.0
  step = 0.1

#initial values:  
  FFox = 0.2
  KKox = 0.2
  FFred = 0.1
  KKred = 0.1
  Ffree = 0.0
  Kfree = 0.0
  FKred = 0.0
  FKox = 0.0
  MEAthiol = 15.0 
  CO2 = 0.12
  yO2P = 0.209
  CYST = 0.0
  VL = 7.2

#constants:
  VLinitial = 7.2
  Vtotalvessel = 10.0
  AgitatorSpeed = 400.0
  AgitatorDiameter = 6.0
  AgitatorPowerNumber = 2.1
  pH = 7.4
  k1red =  0.05604
  k1ox = 0.0108
  k2Fd =  1.35
  k2Fa =  110400000.0
  k2Kd =  0.04038
  k2Ka =  120000000.0
  k3FKa =  181200000.0
  k3FKd =  0.01188
  k4ox =  0.0108
  k4red =  0.05604
  kthiolox = 0.005
  krcyst = 0.0
  percentO2saturation = 100.0
  pKa2MEA = 8.18

#parameters:
  qin = 1.0
  yO2in = 0.21
  H = 1.072069378
  T = 300.0
  R = 0.082
  P = 1.0
  TimeToSwitch = 135.0
`;

const CONTROL_TAG = '#';
const CONTROL_SEP = ':';
const EQUAL_SIGN = '=';
const DIV_SIGN = '/';
const SERVICE = '_';
const DEFAULT_TOL = '0.00005';

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
  OUTPUT = '//output: dataframe df',
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
  OUTPUT = `let df = call.getParamValue('df');`,
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
function processed(expression: string): string {
  let proc = expression;

  //MATH_FUNCS.forEach((func) => proc = proc.replace(func, `Math.${func}`));

  return proc;
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
    deqs.set(name, processed(line.slice(eqIdx + 1).trim()));
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
    exprs.set(line.slice(0, eqIdx).replace(' ', ''), processed(line.slice(eqIdx + 1).trim()));
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
      console.log(firstLine);
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

/** */
function getScript(ivp: IVP, toAddViewers = true, toAddEditor = true): string[] {
  // 1. Create annotation
  const res = getAnnot(ivp, toAddViewers, toAddEditor);

  // 2. Constants lines
  if (ivp.consts !== null) {
    res.push(SCRIPT.CONSTS);
    ivp.consts.forEach((val, key, map) => res.push(`const ${key} = ${val};`));
  }

  // 3. The problem definition lines
  res.push(SCRIPT.ODE_COM);
  res.push(SCRIPT.ODE);
  res.push(`${SCRIPT.SPACE}name: '${ivp.name}',`);

  // 3.1) argument
  const t = ivp.arg.name;
  const t0 = `${SERVICE}${t}0`;
  const t1 = `${SERVICE}${t}1`;
  const h = `${SERVICE}h`;
  res.push(`${SCRIPT.SPACE}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);

  const names = ivp.deqs.solutionNames;

  // 3.2) initial values
  res.push(`${SCRIPT.SPACE}initial: [${names.join(', ')}],`);

  // 3.3) the right-hand side of the problem
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

  // 3.4) final lines of the problem specification
  res.push(`${SCRIPT.SPACE}tolerance: ${ivp.tolerance},`);
  res.push(`${SCRIPT.SPACE}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
  res.push('};');

  // 4. Math functions
  if (ivp.usedMathFuncs.length > 0) {
    res.push(SCRIPT.MATH_FUNC_COM);
    ivp.usedMathFuncs.forEach((i) => res.push(`const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
  }

  // 5. Math constants
  if (ivp.usedMathConsts.length > 0) {
    res.push(SCRIPT.MATH_CONST_COM);
    ivp.usedMathConsts.forEach((i) => res.push(`const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
  }

  // 6. The 'call solver' lines
  res.push(SCRIPT.SOLVER_COM);
  res.push(SCRIPT.SOLVER);
  res.push(SCRIPT.PREPARE);
  res.push(SCRIPT.CALL);
  res.push(SCRIPT.OUTPUT);

  return res;
} // getScript


/*const ivp = getIVP(bio);
console.log(ivp);

console.log();

const scriptLines = getScript(ivp);

for (const line of scriptLines)
  console.log(line);*/
