let s = '{units: min; caption: Initial; category: Time; min: 10; max: 20; showSlider: true} [Initial time of simulation]';
let descr: string | null = null;

console.log(s);

let posOpen = s.indexOf('[');
let posClose;
console.log(posOpen);

if (posOpen !== -1) {
  posClose = s.indexOf(']');

  if (posClose === -1)
    throw new Error('Mistake!');

  descr = s.slice(posOpen + 1, posClose);

  s = s.slice(0, posOpen);
}

posOpen = s.indexOf('{');
posClose = s.indexOf('}');

s = s.slice(posOpen + 1, posClose);

console.log(s.split(';'));
console.log(descr);

let obj = {};

let pos: number;
let key: string;
let val;

const strToVal = (s: string) => {
  let num = Number(s);

  if (!isNaN(num))
    return num;
  
  if (s === 'true')
    return true;

  if (s === 'false')
    return false;

  return s;
};

s.split(';').forEach((str) => {
  pos = str.indexOf(':');

  if (pos === -1)
    throw new Error('Mistake');

  key = str.slice(0, pos).trim();
  val = str.slice(pos + 1).trim();

  console.log(key);
  console.log();

  // @ts-ignore
  obj[key] = strToVal(val);
});

console.log(obj);
