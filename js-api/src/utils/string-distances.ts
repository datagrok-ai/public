// Dependency-free on purpose: imported by both StringUtils (helpers.ts) and headless test stubs.

/** Normalized Levenshtein distance in [0, 1]. */
export function levenshteinDistance(a: string, b: string): number {
  const longest = Math.max(a.length, b.length);
  if (longest === 0)
    return 0;
  let prev: number[] = [];
  for (let k = 0; k <= b.length; k++)
    prev.push(k);
  for (let i = 1; i <= a.length; i++) {
    const row = [i];
    for (let k = 1; k <= b.length; k++)
      row.push(Math.min(prev[k] + 1, row[k - 1] + 1, prev[k - 1] + (a[i - 1] === b[k - 1] ? 0 : 1)));
    prev = row;
  }
  return prev[b.length] / longest;
}

/** 1 − Jaro-Winkler similarity (standard definition, p = 0.1, prefix ≤ 4). */
export function jaroWinklerDistance(a: string, b: string): number {
  const j = jaro(a, b);
  const limit = Math.min(4, a.length, b.length);
  let prefix = 0;
  while (prefix < limit && a[prefix] === b[prefix])
    prefix++;
  return 1 - (j + prefix * 0.1 * (1 - j));
}

function jaro(a: string, b: string): number {
  if (a.length === 0 && b.length === 0)
    return 1;
  if (a.length === 0 || b.length === 0)
    return 0;
  const window = Math.max(0, Math.floor(Math.max(a.length, b.length) / 2) - 1);
  const aMatched: boolean[] = new Array(a.length).fill(false);
  const bMatched: boolean[] = new Array(b.length).fill(false);
  let matches = 0;
  for (let i = 0; i < a.length; i++) {
    const to = Math.min(b.length - 1, i + window);
    for (let k = Math.max(0, i - window); k <= to; k++) {
      if (!bMatched[k] && a[i] === b[k]) {
        aMatched[i] = true;
        bMatched[k] = true;
        matches++;
        break;
      }
    }
  }
  if (matches === 0)
    return 0;
  let transpositions = 0;
  let k = 0;
  for (let i = 0; i < a.length; i++) {
    if (!aMatched[i])
      continue;
    while (!bMatched[k])
      k++;
    if (a[i] !== b[k])
      transpositions++;
    k++;
  }
  return (matches / a.length + matches / b.length +
    (matches - transpositions / 2) / matches) / 3;
}
