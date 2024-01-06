const OVERHANG_SUFFIX = '(o)';

export function isOverhangNucleotide(modification: string): boolean {
  return modification.endsWith(OVERHANG_SUFFIX);
}

// export function isOverhang(modification: string): boolean {
//   const overhangSuffix = '(o)';
//   return modification.endsWith(overhangSuffix);
// }
