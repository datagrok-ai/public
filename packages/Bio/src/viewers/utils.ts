
export function numbersWithinMaxDiff(a: number, b: number, maxDiff: number): boolean {
  return Math.abs(a - b) <= maxDiff;
}
