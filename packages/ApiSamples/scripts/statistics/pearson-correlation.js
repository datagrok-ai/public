//name: Pearson
//description: Computes Pearson correlation for two numerical columns
//language: javascript
//input: dataframe df
//input: column c1
//input: column c2
//output: double corr

let a1 = c1.getRawData();
let a2 = c1.getRawData();
let n = 0;
let sumX = 0, sumY = 0;
let sumXX = 0, sumYY = 0;
let sumXY = 0;

for (let i = 0; i < c1.length; i++) {
  const x = a1[i];
  const y = a2[i];

  if (x == null || y == null)
    continue;
  const xf = Number(x);
  const yf = Number(y);
  if (!Number.isFinite(xf) || !Number.isFinite(yf))
    continue;

  n++;
  sumX += xf;
  sumY += yf;
  sumXX += xf * xf;
  sumYY += yf * yf;
  sumXY += xf * yf;
}

if (n >= 2) {
  const num = n * sumXY - sumX * sumY;
  const den = Math.sqrt((n * sumXX - sumX * sumX) * (n * sumYY - sumY * sumY));
  corr = den === 0 ? NaN : (num / den);
}
