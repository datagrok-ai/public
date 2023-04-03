export function mapToFixed(ar: Float32Array | number[]) {
  return Array.from(ar).map((d) => Number(d.toFixed(3)));
}
