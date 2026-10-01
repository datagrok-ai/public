//name: FedBatchTiter
//friendlyName: Fed-batch titer (example)
//description: Synthetic fed-batch cell culture model for the example biologics CMC portfolio (EX-BIOC). Logistic growth with temperature- and feed-dependent productivity. Invented for the example; not fitted to any process.
//language: javascript
//meta.role: model
//meta.domain: Bioprocess
//input: double feedRate = 3 {caption: Feed rate; units: % v/v per day; min: 0; max: 10; category: Process}
//input: double temperature = 36.5 {caption: Temperature; units: °C; min: 30; max: 39; category: Process}
//input: int days = 14 {caption: Culture duration; units: days; min: 1; max: 30; category: Process}
//input: double seedingDensity = 0.5 {caption: Seeding density; units: 1e6 cells/mL; min: 0.1; max: 5; category: Inoculum}
//output: double titer {caption: Final titer; units: g/L}
//output: double peakVcd {caption: Peak viable cell density; units: 1e6 cells/mL}
//output: double finalViability {caption: Final viability; units: %}
//output: dataframe timeCourse {caption: Time course; viewer: Line chart(x: "Day")}

// Everything below is made up for the example: round numbers chosen to give plausible shapes.
const dt = 0.05;
const steps = Math.round(days / dt);
const f = Math.max(0, feedRate);

const muMax = 0.9 * Math.exp(-Math.pow((temperature - 36.8) / 4, 2)); // 1/day
const capacity = Math.max(2, 8 + 6 * f - 0.4 * f * f);                   // 1e6 cells/mL
const qp = 18 * (1 + 0.08 * (37 - temperature)) * (0.6 + 0.4 * Math.min(f, 5) / 5); // pg/cell/day

let live = seedingDensity;
let dead = 0;
let product = 0; // g/L
let peak = live;

const n = Math.floor(days) + 1;
const dayCol = new Float32Array(n);
const vcdCol = new Float32Array(n);
const viaCol = new Float32Array(n);
const titerCol = new Float32Array(n);
const record = (i) => {
  dayCol[i] = i;
  vcdCol[i] = live;
  viaCol[i] = 100 * live / (live + dead);
  titerCol[i] = product;
};
record(0);

for (let s = 1; s <= steps; s++) {
  const t = s * dt;
  const death = (0.02 + 0.02 * Math.max(0, t - 8)) * (temperature / 37);
  const growth = muMax * live * (1 - live / capacity);
  const dying = death * live;
  live = Math.max(0, live + (growth - dying) * dt);
  dead = Math.max(0, dead + (dying - 0.05 * dead) * dt);
  product += qp * live * 1e-3 * dt; // pg/cell/day x 1e6 cells/mL -> g/L/day
  peak = Math.max(peak, live);
  if (Math.abs(t - Math.round(t)) < dt / 2 && Math.round(t) < n)
    record(Math.round(t));
}

titer = Math.round(product * 1000) / 1000;
peakVcd = Math.round(peak * 100) / 100;
finalViability = Math.round(1000 * live / (live + dead)) / 10;
timeCourse = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('Day', dayCol),
  DG.Column.fromFloat32Array('VCD', vcdCol),
  DG.Column.fromFloat32Array('Viability', viaCol),
  DG.Column.fromFloat32Array('Titer', titerCol),
]);
timeCourse.name = 'Fed-batch time course';
