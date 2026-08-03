/**
 * Shared printing utilities for stats examples.
 *
 * Mirrors the `printProblem` / `printResult` style used in the
 * optimizer examples (`src/optimization/single-objective/examples/`).
 */

/** Top-level section header. */
export function printSection(title: string): void {
  console.log(`\n========== ${title} ==========`);
}

/** Subsection separator with a label. */
export function printCase(label: string): void {
  console.log(`\n--- ${label} ---`);
}

/** Indented multi-line block describing the problem and its expected result. */
export function printProblem(lines: readonly string[]): void {
  for (const line of lines)
    console.log(`  ${line}`);
}

/** One key-value line, indented and aligned. */
export function printKv(key: string, value: string | number, width = 20): void {
  console.log(`  ${key.padEnd(width)} ${value}`);
}

/** Format a numeric value with `digits` decimal places, or `null` as 'null'. */
export function fmt(x: number | null | undefined, digits = 4): string {
  if (x === null || x === undefined) return 'null';
  if (!Number.isFinite(x)) return String(x);
  return x.toFixed(digits);
}
