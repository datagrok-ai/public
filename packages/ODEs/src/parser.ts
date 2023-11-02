const eqs = `
Some stupid text
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

#argument: t [min]
  start = 0.0


  finish = 5.0



  step = 0.1

#initial values:
  x = 2.0 [mu, x(t)]
  y = 0.0 [mu, y(t)]

#constants:
  const1 = 1.0 [mu, constanta 1]
  const2 = 3.0 [mu, cOnSTanTA 2]

#parameters:
  param1 = 1.0 [mu, paRaM 1]
  param2 = -1.0 [mu, PArAm 2]

#tolerance: 0.00005`;

const CONTROL_TAG = '#';
const CONTROL_SEP = ':';
const EQUAL_SIGN = '='

enum CONTROL_EXPR {
  NAME = `${CONTROL_TAG}name`,
  DIF_EQ = `${CONTROL_TAG}differential equations`,
}

enum ERROR_MSG {
  INCORRECT_ODES = 'incorrect definition of the system of ordinary differential equations'
}

type Block = {
  begin: number,
  end: number
}

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
function getODEs(text: string) {
  // 0. Split text into lines
  const lines = eqs.split('\n').filter((s) => s !== '').map((s) => s.trimStart());
  const linesCount = lines.length;  

  // 1. Skip first lines without the control tag  
  let idx = 0;
 
  while (!lines[idx].startsWith(CONTROL_TAG)) {
    ++idx;
    
    if (idx === linesCount)
      throw new Error(ERROR_MSG.INCORRECT_ODES);
  }

  let i = idx;
  
  while (i < linesCount) {
    console.log(`line ${i}:${lines[i]}`);
    ++i;
  }

  // 2. Get blocks limits
  let beg = idx;
  ++idx;
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

  let sep: number;

  // 3. Process blocks
  for (const block of blocks) {
    const firstLine = lines[block.begin];

    if (firstLine.startsWith(CONTROL_EXPR.NAME)) {
      sep = firstLine.indexOf(CONTROL_SEP);
      const name = firstLine.slice(sep + 1).trim();
      console.log(name);
    }
    else if (firstLine.startsWith(CONTROL_EXPR.DIF_EQ)) {
      console.log(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
    }
  }

  console.log('end');
}

getODEs(eqs);
