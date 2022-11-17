export function isOverhang(modification: string): boolean {
  return modification.slice(-3) == '(o)';
}
