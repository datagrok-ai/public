// Global flag is for replaceAll
const helmGapStartRe = /\{(\*\.)+/g;
const helmGapIntRe = /\.(\*\.)+/g;
const helmGapEndRe = /(\.\*)+\}/g;

export function removeGapsFromHelm(srcHelm: string): string {
  return srcHelm.replaceAll(helmGapStartRe, '{').replaceAll(helmGapIntRe, '.')
    .replaceAll(helmGapEndRe, '}').replace('{*}', '{}');
}
