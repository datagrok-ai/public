import { erf } from 'mathjs';

/* ---------- Basic stats ---------- */

export function mean(x: number[]): number {
  return x.reduce((a, b) => a + b, 0) / x.length;
}

export function variance(x: number[], m?: number): number {
  const mu = m ?? mean(x);
  return x.reduce((s, v) => s + (v - mu) ** 2, 0) / (x.length - 1);
}

export function std(x: number[], m?: number): number {
  return Math.sqrt(variance(x, m));
}

/* ---------- Gaussian ---------- */

export function gaussian(x: number, mu: number, sigma: number): number {
  if (sigma === 0) return 0;
  const z = (x - mu) / sigma;
  return Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2 * Math.PI));
}

/* ---------- Welch t-test ---------- */

export function welchTTest(
  x: number[],
  y: number[]
): { t: number; p: number } {
  const mx = mean(x);
  const my = mean(y);
  const vx = variance(x, mx);
  const vy = variance(y, my);

  const nx = x.length;
  const ny = y.length;

  const t =
    (mx - my) / Math.sqrt(vx / nx + vy / ny);

  // Degrees of freedom (Welch–Satterthwaite)
  const df =
    ((vx / nx + vy / ny) ** 2) /
    ((vx * vx) / (nx * nx * (nx - 1)) +
     (vy * vy) / (ny * ny * (ny - 1)));

  const p = 2 * (1 - normalCDF(Math.abs(t)));

  return { t, p };
}

/* ---------- Normal CDF ---------- */

export function normalCDF(x: number): number {
  return 0.5 * (1 + erf(x / Math.sqrt(2)));
}

/* ---------- Pearson correlation ---------- */

export function pearson(x: number[], y: number[]): number {
  const mx = mean(x);
  const my = mean(y);

  let num = 0;
  let dx = 0;
  let dy = 0;

  for (let i = 0; i < x.length; i++) {
    const a = x[i] - mx;
    const b = y[i] - my;
    num += a * b;
    dx += a * a;
    dy += b * b;
  }

  return num / Math.sqrt(dx * dy);
}
