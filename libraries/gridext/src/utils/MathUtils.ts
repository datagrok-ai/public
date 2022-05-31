export function isStrictInt(val : number) : boolean {
  const b = Math.floor(val) === val;
  return b;
}
