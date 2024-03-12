const OVERHANG_POSTFIX = '(o)';

export function isOverhangNucleotide(modification: string): boolean {
  return modification.endsWith(OVERHANG_POSTFIX);
}
