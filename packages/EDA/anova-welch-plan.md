# План: добавление Welch ANOVA в `packages/EDA/src/anova/`

**Цель:** закрыть дыру в продукте — пользователи с гетероскедастичными данными сейчас упираются в `NON_EQUAL_VARIANCES` error и не могут провести анализ. Welch ANOVA штатно работает с неравными дисперсиями.

**Принцип:** не переписать ANOVA, а закрыть дыру. Welch — добавление, а не замена. Существующий Fisher остаётся.

**Сопутствующее улучшение:** попутно чиним две численные проблемы в существующем Fisher (Welford для variance, ibeta для p-value в хвосте) — без breaking change в выдаваемых значениях для большинства данных, но с корректной работой на стиффных входах и больших F.

**Срок:** 2–3 дня на код + 1 день на тесты с NIST.

---

## 1. Изменения в UI (`anova-ui.ts`)

### 1.1. Новый input — Method

Между `Feature` и `Alpha` в существующем диалоге добавляется Radio-input для выбора метода:

```ts
const methodProp = DG.Property.fromOptions({
  name: 'method',
  inputType: 'Radio',
  choices: ['Welch', 'Fisher'],
});
const methodInput = ui.input.forProperty(methodProp, { method: 'Welch' });
// Точный API создания инпута из Property уточняется при имплементации;
// если ui.input.forProperty недоступен — собрать через ui.input.form
// с массивом из одной property.
```

Дефолт — **Welch**.

### 1.2. Tooltip метода

Биндится на сам input (не на caption):

```ts
const methodTooltip = ui.divV([
  ui.h2('Method'),
  ui.markdown(
    'Fisher — classical ANOVA. Assumes **equal variances** across groups; ' +
    'more powerful when this assumption holds.\n\n' +
    'Welch — robust to **unequal variances** across groups. ' +
    'Recommended as the default when group variances are unknown or differ.'
  ),
]);
ui.tooltip.bind(methodInput.root, methodTooltip);
```

### 1.3. Гейт Run-кнопки

Вариансный чек переезжает из тела `oneWayAnova` в UI — вызывается на изменение **любого** input'а (factor, feature, method, alpha) и при инициализации диалога.

```ts
let currentMethod: 'Welch' | 'Fisher' = 'Welch';

function updateRunButtonState(): void {
  // 1. Проверка alpha (как сейчас)
  if (significance <= SIGNIFICANCE.INFIMUM || significance >= SIGNIFICANCE.SUPREMUM) {
    runBtn.disabled = true;
    ui.tooltip.bind(runBtn, 'Alpha must be strictly between 0 and 1.');
    return;
  }

  // 2. Проверка возможности построить FactorizedData
  //    (single category, no variation — те же ошибки, что бросает конструктор)
  let varEqual: boolean;
  try {
    const uniqueCount = factor!.stats.uniqueCount;
    if (uniqueCount < 2) throw new Error(ERROR_MSG.SINGLE_FACTOR);
    const factorized = new FactorizedData(factor!, feature!, uniqueCount);
    varEqual = factorized.areVarsEqual(significance);
  } catch (err) {
    runBtn.disabled = true;
    ui.tooltip.bind(runBtn, (err as Error).message);
    return;
  }

  // 3. Гейт по комбинации (метод, равенство дисперсий)
  if (currentMethod === 'Fisher' && !varEqual) {
    runBtn.disabled = true;
    ui.tooltip.bind(
      runBtn,
      'Variances differ significantly between groups. ' +
      "Fisher's ANOVA requires equal variances — switch Method to Welch."
    );
    return;
  }

  // 4. Всё ок
  runBtn.disabled = false;
  ui.tooltip.bind(runBtn, 'Perform analysis of variances');
}
```

Привязка к input'ам — через `onValueChanged` каждого. Все четыре input'а вызывают `updateRunButtonState()`. После сборки диалога и **до** `dlg.show()` — один вызов `updateRunButtonState()` для инициализации.

### 1.4. Run-handler

Получает выбранный метод и передаёт в `oneWayAnova`:

```ts
dlg.addButton('Run', () => {
  dlg.close();
  try {
    const res = oneWayAnova(factor!, feature!, significance, {
      method: currentMethod,
      toValidate: false,  // гейт уже отработал на UI
    });
    addVizualization(df, factor!.name, feature!.name, res);
  } catch (error) {
    // существующий error-handling без изменений
  }
});
```

### 1.5. Рендеринг таблицы (`getAnovaGrid`)

Функция становится разветвлённой по `report.method`:

**Для Fisher** (`report.method === 'Fisher'`): текущая таблица 3×7 — Between/Within/Total × SS/DF/MS/F/F-critical/p-value. **Без изменений.**

**Для Welch** (`report.method === 'Welch'`): таблица 1×4 с релевантными для W-теста полями — `F`, `df₁` (целое), `df₂` (дробное), `p-value`. Никаких SS/MS/Total — концептуально их у Welch нет. `F-critical` отображается рядом с `F` как отдельная колонка или в conclusion-блоке (детали рендеринга — на имплементации).

```ts
function getAnovaGrid(report: OneWayAnovaReport): DG.Grid {
  if (report.method === 'Fisher')
    return getFisherGrid(report.anovaTable, report.fCritical, report.significance);
  else
    return getWelchGrid(report.anovaTable, report.fCritical, report.significance);
}
```

`getFisherGrid` — это текущее тело `getAnovaGrid` без изменений. `getWelchGrid` — новая функция, 1 строка × 5 колонок (`F`, `df₁`, `df₂`, `F-critical`, `p-value`), с tooltip'ами per column в стиле существующего.

### 1.6. Conclusion-блок (`addVizualization`)

Заголовок таба остаётся `'ANOVA'` (по согласованию). Conclusion остаётся таким же: `F > F-critical` → reject null. Это работает одинаково для обоих методов.

---

## 2. Изменения в математике (`anova-tools.ts`)

### 2.1. Новые типы

```ts
/** One-way ANOVA Welch result. Different shape from OneWayAnova because
 *  Welch's W-test does not have a SS/MS decomposition. */
type WelchAnova = {
  /** F-statistic, Welch's W */
  fStat: number;
  /** Degrees of freedom: numerator (k−1). Integer. */
  dfBn: number;
  /** Welch–Satterthwaite degrees of freedom: denominator. Fractional. */
  dfWn: number;
  /** p-value, P(F_{dfBn, dfWn} > fStat) */
  pValue: number;
  /** Per-group means (for display). */
  groupMeans: Float64Array;
  /** Per-group unbiased variances (for display). */
  groupVariances: Float64Array;
  /** Per-group sizes (for display). */
  groupSizes: Int32Array;
};

export type OneWayAnovaReport =
  | { method: 'Fisher'; anovaTable: OneWayAnova; fCritical: number; significance: number; }
  | { method: 'Welch';  anovaTable: WelchAnova;  fCritical: number; significance: number; };

export interface OneWayAnovaOptions {
  /** Default: 'Welch'. */
  method?: 'Fisher' | 'Welch';
  /** Default: true. Only meaningful for Fisher. */
  toValidate?: boolean;
}
```

Тип `OneWayAnova` (Fisher) **не меняется** — обратная совместимость по полям.

### 2.2. Изменение сигнатуры `oneWayAnova`

```ts
export function oneWayAnova(
  categores: CatCol,
  values: NumCol,
  alpha: number,
  opts: OneWayAnovaOptions = {},
): OneWayAnovaReport {
  checkSignificanceLevel(alpha);
  const method = opts.method ?? 'Welch';
  const toValidate = opts.toValidate ?? true;

  const uniqueCount = categores.stats.uniqueCount;
  if (uniqueCount < 2) throw new Error(ERROR_MSG.SINGLE_FACTOR);

  const factorized = new FactorizedData(categores, values, uniqueCount);

  if (method === 'Fisher') {
    if (toValidate && !factorized.areVarsEqual(alpha))
      throw new Error(ERROR_MSG.NON_EQUAL_VARIANCES);
    const anova = factorized.getOneWayAnova();
    return {
      method: 'Fisher',
      anovaTable: anova,
      fCritical: jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn),
      significance: alpha,
    };
  } else {
    const anova = factorized.getWelchAnova();
    return {
      method: 'Welch',
      anovaTable: anova,
      fCritical: jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn),
      significance: alpha,
    };
  }
}
```

**Существующий старый callsite в `anova-ui.ts`** — `oneWayAnova(factor, feature, significance)` — будет работать (дефолты применятся), но получит **Welch** вместо Fisher. Это согласованное продуктовое решение.

**Старый позиционный четвёртый аргумент** `toValidate: boolean` в существующем тесте `anova-tests.ts:48` — нужно обновить вызов на новый API: `oneWayAnova(categories!, features!, ALPHA, { toValidate: TO_VALIDATE })`. Это единственная правка в существующих тестах, кроме обновления эталонов.

### 2.3. Welford в `FactorizedData.setStats`

Заменяем `Float64Array sums` и `sumsOfSquares` на `Float64Array means` и `m2` (sum of squared deviations).

**Внутренний цикл становится:**

```ts
// Numeric branch (categories != BOOL)
for (let i = 0; i < size; ++i) {
  const cat = cats[i];
  if ((cat === categoriesNull) || (vals[i] === featuresNull)) {
    ++this.nullsCount;
    continue;
  }
  const x = vals[i];
  const n = ++subSampleSizes[cat];          // n_cat after increment
  const delta = x - means[cat];
  means[cat] += delta / n;
  const delta2 = x - means[cat];
  m2[cat] += delta * delta2;
}
```

Аналогично для bool-ветки (с bit-packing — без изменений по структуре, только тело внутреннего цикла переписывается на Welford).

Private-поля класса:
- `sums: Float64Array` → `means: Float64Array`;
- `sumsOfSquares: Float64Array` → `m2: Float64Array`;
- `subSampleSizes: Int32Array` — без изменений.

### 2.4. Обновление `getVariance` и `areVarsEqual`

Старый интерфейс через `SampleData { sum, sumOfSquares, size }` не подходит — теперь дисперсия группы вычисляется напрямую через `m2[k] / (subSampleSizes[k] - 1)`. Переписываем:

```ts
type GroupStats = { mean: number; m2: number; size: number; };

/** Unbiased sample variance from Welford accumulator. */
export function getVariance(g: GroupStats): number {
  if (g.size <= 1) return 0;
  return g.m2 / (g.size - 1);
}

function areVarsEqual(x: GroupStats, y: GroupStats, alpha: number): boolean {
  checkSignificanceLevel(alpha);
  const xVar = getVariance(x);
  const yVar = getVariance(y);
  if (xVar === 0 || yVar === 0) return xVar === yVar;
  const fStat = xVar / yVar;
  const fCrit = jStat.centralF.inv(1 - alpha, x.size - 1, y.size - 1);
  return fStat < fCrit;
}

// FactorizedData.areVarsEqual:
public areVarsEqual(alpha: number): boolean {
  const K = this.catCount;
  if (K === 1) return true;
  const first: GroupStats = { mean: this.means[0], m2: this.m2[0], size: this.subSampleSizes[0] };
  for (let i = 1; i < K; ++i) {
    const next: GroupStats = { mean: this.means[i], m2: this.m2[i], size: this.subSampleSizes[i] };
    if (!areVarsEqual(first, next, alpha)) return false;
  }
  return true;
}
```

`SampleData` тип становится не нужен — удаляется.

### 2.5. Существующий `getOneWayAnova` (Fisher) — переписать через means/m2

Численно устойчивая версия:

```ts
public getOneWayAnova(): OneWayAnova {
  let K = this.catCount;
  let N = 0;
  let nonEmpty = 0;

  // Подсчёт N и непустых категорий
  for (let i = 0; i < K; ++i) {
    if (this.subSampleSizes[i] !== 0) {
      N += this.subSampleSizes[i];
      ++nonEmpty;
    }
  }
  K = nonEmpty;

  if (K === 1) throw new Error(ERROR_MSG.SINGLE_FACTOR);
  if (N === K) throw new Error(ERROR_MSG.CATS_EQUAL_SIZE);

  // Grand mean (взвешенное по размерам групп)
  let grandMeanNumer = 0;
  for (let i = 0; i < this.catCount; ++i)
    if (this.subSampleSizes[i] !== 0)
      grandMeanNumer += this.subSampleSizes[i] * this.means[i];
  const grandMean = grandMeanNumer / N;

  // ssWn = Σ m2[k] — тождество Welford
  // ssBn = Σ n[k] · (means[k] − grandMean)²
  // ssTot = ssBn + ssWn (точное тождество, не пересчитываем независимо)
  let ssWn = 0;
  let ssBn = 0;
  for (let i = 0; i < this.catCount; ++i) {
    const n = this.subSampleSizes[i];
    if (n === 0) continue;
    ssWn += this.m2[i];
    const diff = this.means[i] - grandMean;
    ssBn += n * diff * diff;
  }
  const ssTot = ssBn + ssWn;

  if (ssWn === 0) throw new Error(ERROR_MSG.NO_FEATURE_VARIATION_WITHIN_GROUPS);

  const dfBn = K - 1;
  const dfWn = N - K;
  const dfTot = N - 1;
  const msBn = ssBn / dfBn;
  const msWn = ssWn / dfWn;
  const fStat = msBn / msWn;

  return {
    ssBn, ssWn, ssTot,
    dfBn, dfWn, dfTot,
    msBn, msWn,
    fStat,
    pValue: fSurvival(fStat, dfBn, dfWn),
  };
}
```

Формула SS-разложения **математически эквивалентна** текущей, но численно устойчива (нет разности больших близких чисел `(Σx²) − (Σx)²/N`).

### 2.6. Новый метод `getWelchAnova`

```ts
public getWelchAnova(): WelchAnova {
  // Собираем непустые группы и считаем means/variances/sizes per group.
  const groupMeans: number[] = [];
  const groupVars: number[] = [];
  const groupSizes: number[] = [];

  for (let i = 0; i < this.catCount; ++i) {
    const n = this.subSampleSizes[i];
    if (n === 0) continue;
    if (n < 2) continue;                       // Welch requires n_i ≥ 2
    const v = this.m2[i] / (n - 1);
    if (v <= 0)
      throw new Error(ERROR_MSG.NO_FEATURE_VARIATION_WITHIN_GROUPS);
    groupMeans.push(this.means[i]);
    groupVars.push(v);
    groupSizes.push(n);
  }

  const K = groupMeans.length;
  if (K < 2) throw new Error(ERROR_MSG.SINGLE_FACTOR);

  // Welch (1951) formulae.
  let W = 0;
  let wxSum = 0;
  for (let i = 0; i < K; i++) {
    const w = groupSizes[i] / groupVars[i];
    W += w;
    wxSum += w * groupMeans[i];
  }
  const wGrandMean = wxSum / W;

  let numer = 0;
  let lambda = 0;
  for (let i = 0; i < K; i++) {
    const w = groupSizes[i] / groupVars[i];
    const diff = groupMeans[i] - wGrandMean;
    numer += w * diff * diff;
    const r = 1 - w / W;
    lambda += (r * r) / (groupSizes[i] - 1);
  }
  numer /= (K - 1);
  const denom = 1 + 2 * (K - 2) / (K * K - 1) * lambda;

  const fStat = numer / denom;
  const dfBn = K - 1;
  const dfWn = (K * K - 1) / (3 * lambda);

  return {
    fStat,
    dfBn,
    dfWn,
    pValue: fSurvival(fStat, dfBn, dfWn),
    groupMeans: Float64Array.from(groupMeans),
    groupVariances: Float64Array.from(groupVars),
    groupSizes: Int32Array.from(groupSizes),
  };
}
```

### 2.7. Numerical-устойчивый p-value helper

Добавляется в `anova-tools.ts` рядом с другими утилитами:

```ts
/**
 * P(F > f | df1, df2) для центрального F-распределения.
 * Численно устойчив в хвосте — не использует `1 − cdf`.
 *
 * Identity: P(F > f) = I_x(df2/2, df1/2),  где x = df2 / (df2 + df1·f).
 * `jStat.ibeta` реализована через continued fraction Lentz'а — точность
 * сохраняется при малых x (= больших f).
 *
 * Применяется и в Fisher (getOneWayAnova), и в Welch (getWelchAnova).
 */
function fSurvival(f: number, df1: number, df2: number): number {
  if (!Number.isFinite(f) || f <= 0) return f > 0 ? 0 : 1;
  const x = df2 / (df2 + df1 * f);
  return jStat.ibeta(x, df2 / 2, df1 / 2);
}
```

**Инвариант:** в файле `anova-tools.ts` после правок не должно остаться ни одной строки `1 - jStat.centralF.cdf(...)`. Проверяется грепом при ревью.

`jStat.centralF.inv(1 - alpha, ...)` в `oneWayAnova` для вычисления `fCritical` — **не трогаем**, это inverse CDF в "нормальной" области (alpha = 0.05, не в хвосте), точности jStat хватает.

---

## 3. Изменения в тестах (`anova-tests.ts`)

### 3.1. Существующий `Correctness` тест

Один точечный апдейт сигнатуры в вызове:

```ts
// БЫЛО:
const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, TO_VALIDATE);

// СТАЛО:
const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, {
  method: 'Fisher',
  toValidate: TO_VALIDATE,
});
// Явно method: 'Fisher', потому что зашитые ожидания EXPECTED.* — Fisher-значения.
```

Все эталоны (`EXPECTED.SS_BN = 63.333` и т.д.) остаются. После правки Welford значения для этого датасета совпадут с текущими в пределах ERR=0.01 — данные не стиффные, наивная и Welford-формулы дают одно и то же до 4–5 знаков.

`ERR = 0.01` оставляем как есть.

### 3.2. Существующий `Performance` тест

Один апдейт сигнатуры (как выше). Возможно, придётся поднять `TIMEOUT` с 4000 до 5000–6000 ms — Welford примерно на 20–30% медленнее наивного цикла. Решается по факту прогона.

### 3.3. Новый тест: `Correctness — Welch on validation data`

Тот же 15-точечный датасет, что в существующем `Correctness`, но с `method: 'Welch'` и эталонами от scipy `f_oneway(equal_var=False)`:

```ts
import { FIXTURES } from './anova-fixtures';

test('Correctness: Welch on validation data', async () => {
  const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, {
    method: 'Welch',
    toValidate: false,
  });
  if (analysis.method !== 'Welch') throw new Error('Expected Welch report');
  const w = analysis.anovaTable;
  const f = FIXTURES.validation_simple.welch;

  expectClose(w.fStat, f.fStat, f.tol);
  expectClose(w.dfBn, f.dfBn, { rtol: 0, atol: 0 });          // integer, exact
  expectClose(w.dfWn, f.dfWn, { rtol: 1e-6, atol: 1e-9 });    // fractional
  expectClose(w.pValue, f.pValue, f.tol);
  expectClose(analysis.fCritical, f.fCritical, f.tol);
}, { timeout: TIMEOUT });
```

Где `expectClose(actual, expected, tol)` — helper для зональных допусков (см. §3.5).

### 3.4. NIST-валидация

Категория `'ANOVA: NIST StRD'` с 5 датасетами (SiRstv, AtmWtAg, SmLs04, SmLs07, SmLs09):

```ts
category('ANOVA: NIST StRD', () => {
  for (const datasetName of ['SiRstv', 'AtmWtAg', 'SmLs04', 'SmLs07', 'SmLs09'] as const) {
    test(`Fisher: ${datasetName}`, async () => {
      const f = FIXTURES.nist[datasetName];
      const cats = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'group', f.factor);
      const vals = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', f.value);
      const r = oneWayAnova(cats, vals, 0.05, { method: 'Fisher', toValidate: false });
      if (r.method !== 'Fisher') throw new Error();
      expectClose(r.anovaTable.fStat, f.fisher.fStat, f.fisher.tol);
      expectClose(r.anovaTable.ssBn, f.fisher.ssBn, f.fisher.tol);
      expectClose(r.anovaTable.ssWn, f.fisher.ssWn, f.fisher.tol);
    });

    test(`Welch: ${datasetName}`, async () => {
      const f = FIXTURES.nist[datasetName];
      const cats = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'group', f.factor);
      const vals = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', f.value);
      const r = oneWayAnova(cats, vals, 0.05, { method: 'Welch', toValidate: false });
      if (r.method !== 'Welch') throw new Error();
      expectClose(r.anovaTable.fStat, f.welch.fStat, f.welch.tol);
      expectClose(r.anovaTable.dfWn, f.welch.dfWn, { rtol: 1e-6, atol: 1e-9 });
    });
  }
});
```

**Зональные допуски (стратегия "c"):**

| Источник эталона | rtol | Применяется к |
|---|---|---|
| NIST certified (Fisher: SiRstv, AtmWtAg) | 1e-10 | Лёгкие датасеты |
| NIST certified (Fisher: SmLs04) | 1e-10 | Average difficulty |
| NIST certified (Fisher: SmLs07, SmLs09) | 1e-6 | Higher difficulty (13 const digits) |
| scipy (Welch: все NIST датасеты) | 1e-9 | Welch reference из scipy 1.16.0 |
| scipy (validation_simple, обе) | 1e-6 | Простой 15-точечный датасет |
| p-value в зоне < 1e-10 | rtol 1e-2 | Tail-zone — jStat.ibeta может расходиться со scipy |

Эти допуски встроены в `FIXTURES.*.tol` — генератор сам решает, какой допуск приложить к каждой паре (датасет, метод).

### 3.5. `expectClose` helper

Минимальный helper в стиле существующего `expect`:

```ts
// в начале anova-tests.ts или в отдельном test-utils.ts
type Tol = { rtol: number; atol: number };

function expectClose(actual: number, expected: number, tol: Tol, label?: string): void {
  if (Number.isNaN(actual) && Number.isNaN(expected)) return;
  if (!Number.isFinite(actual) || !Number.isFinite(expected)) {
    expect(actual, expected, label ?? `non-finite mismatch`);
    return;
  }
  const diff = Math.abs(actual - expected);
  const threshold = tol.atol + tol.rtol * Math.abs(expected);
  const ok = diff <= threshold;
  expect(ok, true,
    `${label ?? 'value'}: |${actual} - ${expected}| = ${diff} > ${threshold} ` +
    `(rtol=${tol.rtol}, atol=${tol.atol})`);
}
```

### 3.6. Канареечные тесты на численные правки

**Welford guard:**

```ts
test('Variance stable on shifted data (Welford regression guard)', async () => {
  const offset = 1e9;
  const cats = DG.Column.fromStrings('cat', [
    'A','A','A','A','A', 'B','B','B','B','B'
  ]);
  const features = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'feature', [
    offset+0.0, offset+1.0, offset+2.0, offset+1.0, offset+0.5,
    offset+10.0, offset+11.0, offset+12.0, offset+11.0, offset+10.5,
  ]);
  // На наивной формуле variance разлетится в мусор; F будет некорректным.
  // На Welford: variance ≈ 0.55 в группе A и ≈ 0.55 в группе B, F огромное.
  const r = oneWayAnova(cats, features, 0.05, { method: 'Welch', toValidate: false });
  if (r.method !== 'Welch') throw new Error();
  expect(r.anovaTable.fStat > 100, true, 'F-statistic too small — variance lost precision');
  expect(Number.isFinite(r.anovaTable.fStat), true, 'F-statistic is NaN/Infinity');
});
```

**Tail SF guard:**

```ts
test('p-value preserves precision in extreme tail (ibeta regression guard)', async () => {
  const cats = DG.Column.fromStrings('cat', [
    'A','A','A','A','A','A','A','A','A','A',
    'B','B','B','B','B','B','B','B','B','B',
  ]);
  const features = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'feature', [
    0, 0.01, 0.02, 0.01, 0, 0.03, 0.02, 0.01, 0.02, 0.01,
    1000, 1000.1, 999.9, 1000.2, 999.8, 1000.1, 999.9, 1000.0, 1000.1, 999.9,
  ]);
  const r = oneWayAnova(cats, features, 0.05, { method: 'Fisher', toValidate: false });
  if (r.method !== 'Fisher') throw new Error();
  // Naive `1 - cdf` дал бы p = 0 ровно.
  expect(r.anovaTable.pValue > 0, true, 'p-value collapsed to 0 in tail');
  expect(r.anovaTable.pValue < 1e-30, true, 'p-value too large for this separation');
});
```

Эти два теста — **регрессионные стражи**. Если кто-то откатит Welford обратно к наивной формуле или вернёт `1 - cdf` — упадут.

---

## 4. Файл эталонов (`anova-fixtures.ts`)

### 4.1. Структура

Сгенерированный TypeScript-модуль, коммитится в репозиторий. NIST .dat файлы — **не коммитятся**, скачиваются генератором при регенерации.

```ts
// anova-fixtures.ts (СГЕНЕРИРОВАНО, не редактировать вручную)
// Generated by tools/generate-anova-fixtures.py from:
//   - scipy 1.16.0
//   - NIST StRD: https://www.itl.nist.gov/div898/strd/anova/

export const FIXTURES = {
  validation_simple: {
    factor: [...],     // category labels (numeric encoding)
    value:  [...],     // feature values
    fisher: { ssBn, ssWn, ssTot, dfBn, dfWn, dfTot, msBn, msWn, fStat, pValue, fCritical, tol },
    welch:  { fStat, dfBn, dfWn, pValue, fCritical, tol },
  },
  nist: {
    SiRstv:  { factor: [...], value: [...], fisher: {...}, welch: {...} },
    AtmWtAg: { factor: [...], value: [...], fisher: {...}, welch: {...} },
    SmLs04:  { factor: [...], value: [...], fisher: {...}, welch: {...} },
    SmLs07:  { factor: [...], value: [...], fisher: {...}, welch: {...} },
    SmLs09:  { factor: [...], value: [...], fisher: {...}, welch: {...} },
  },
} as const;
```

Поле `tol: { rtol: number; atol: number }` встроено в каждый набор эталонов — генератор выставляет его по таблице из §3.4.

**Размер:** SmLs09 — самый большой (18 009 точек), целиком файл ≈ 250 KB. Это нормальный размер для embedded fixtures в TypeScript-тестах.

### 4.2. Генератор

Python-скрипт `tools/generate-anova-fixtures.py` (не коммитится в bundle, живёт рядом с тестами):

1. Скачивает NIST `.dat` файлы из `itl.nist.gov/div898/strd/anova/` (кэш в `.gitignored` папке).
2. Парсит формат NIST (factor + value колонки, sertified F, SS).
3. Прогоняет каждый датасет через:
   - `scipy.stats.f_oneway(*groups)` → Fisher F, p;
   - `scipy.stats.f_oneway(*groups, equal_var=False)` → Welch F, p, df.
4. Для Fisher: использует NIST certified `F`, `SS_between`, `SS_within` как ground truth (когда они есть в .dat файле); scipy — для p-value, который NIST не предоставляет.
5. Для Welch: использует scipy как единственный reference.
6. Применяет таблицу зональных допусков из §3.4 → каждый эталон получает свой `tol`.
7. Генерирует `anova-fixtures.ts`.

В репозитории также — короткий README с инструкцией: «для регенерации запустить `python tools/generate-anova-fixtures.py`; требует `scipy==1.16.0`; коммитить только `anova-fixtures.ts`».

---

## 5. Чек-лист имплементации

### Phase 1 — Математика (1 день)
- [ ] Добавить тип `WelchAnova` в `anova-tools.ts`
- [ ] Превратить `OneWayAnovaReport` в discriminated union
- [ ] Добавить интерфейс `OneWayAnovaOptions`
- [ ] Переписать `FactorizedData.setStats` на Welford (numeric и bool ветки)
- [ ] Обновить private-поля: `sums`→`means`, `sumsOfSquares`→`m2`
- [ ] Удалить тип `SampleData`, ввести `GroupStats`
- [ ] Переписать `getVariance`, `areVarsEqual`, `FactorizedData.areVarsEqual` под новое представление
- [ ] Переписать `FactorizedData.getOneWayAnova` через means/m2 (численно устойчиво, формула SS-разложения через grand mean)
- [ ] Добавить метод `FactorizedData.getWelchAnova`
- [ ] Добавить helper `fSurvival` рядом с другими утилитами
- [ ] Заменить `1 - jStat.centralF.cdf(...)` в `getOneWayAnova` на `fSurvival(...)`
- [ ] Расширить `oneWayAnova` опциями, дефолт method='Welch'
- [ ] Грепом проверить отсутствие `1 - jStat.*.cdf(` в финальном файле

### Phase 2 — UI (0.5 дня)
- [ ] Создать `methodInput` через `DG.Property` с `inputType: 'Radio'`
- [ ] Прикрепить tooltip к `methodInput.root`
- [ ] Реализовать `updateRunButtonState()` helper
- [ ] Подписать все 4 input'а на `onValueChanged → updateRunButtonState()`
- [ ] Вызвать `updateRunButtonState()` перед `dlg.show()`
- [ ] Передать `currentMethod` в `oneWayAnova` через новый API
- [ ] Разветвить `getAnovaGrid` на `getFisherGrid` / `getWelchGrid`
- [ ] Реализовать `getWelchGrid` (1×5 таблица: F, df₁, df₂, F-critical, p-value)
- [ ] Обновить tooltip-карту для Welch-колонок

### Phase 3 — Тесты и фикстуры (1 день)
- [ ] Написать `tools/generate-anova-fixtures.py`
- [ ] Сгенерировать `anova-fixtures.ts` с 5 NIST датасетами + validation_simple
- [ ] Обновить существующий `Correctness` тест на новый API (`{ method: 'Fisher', toValidate: ... }`)
- [ ] Обновить существующий `Performance` тест на новый API
- [ ] Добавить `Correctness: Welch on validation data` тест
- [ ] Добавить категорию `'ANOVA: NIST StRD'` с 10 тестами (5 datasets × 2 methods)
- [ ] Добавить канареечный `Welford guard` тест
- [ ] Добавить канареечный `Tail SF guard` тест
- [ ] Реализовать `expectClose` helper в `anova-tests.ts`

### Phase 4 — Финальная проверка
- [ ] Все существующие тесты ANOVA проходят (с обновлёнными вызовами)
- [ ] Все новые тесты ANOVA проходят
- [ ] Performance тест укладывается в timeout (поднять до 5000-6000ms при необходимости)
- [ ] Запуск ANOVA из UI: Welch по умолчанию, Run-кнопка корректно disables/enables при изменении данных и метода
- [ ] Грепом подтверждено: нет `1 - jStat.*.cdf(` в `anova-tools.ts`
- [ ] Грепом подтверждено: нет `sumOfSquares` в `anova-tools.ts` (Welford rename)
- [ ] Smoke test: открыть демо-датафрейм, прогнать ANOVA на heteroscedastic данных, убедиться что Welch выдаёт результат без error'а

---

## 6. Что НЕ входит

- Post-hoc тесты после Welch (Games-Howell);
- Confidence intervals для group means;
- Robust ANOVA варианты (Brown-Forsythe, James);
- Two-way ANOVA, repeated measures, MANOVA;
- Перенос Welch в sci-comp library как универсальной функции;
- Property-based / mutation тесты;
- Замена существующего Fisher как метода (он остаётся).

---

## 7. Риски и митигация

| Риск | Митигация |
|---|---|
| **Дефолт метода изменился с Fisher на Welch — потенциальный breaking change для внешних вызовов `oneWayAnova(...)`** | Единственный внутренний callsite (`anova-ui.ts` → `runOneWayAnova`) обновляется явно. Если есть внешние package callers — выявить грепом по репозиторию, обновить им вызов с явным `method: 'Fisher'`. |
| **Welford чуть медленнее наивной формулы (~20-30% на пер-элементный цикл)** | Performance тест на 1M строк проверит. TIMEOUT поднять при необходимости. |
| **jStat.ibeta может расходиться со scipy в глубоком хвосте (p < 1e-50)** | NIST-тесты не уходят в такие хвосты (F ~ 1.5–21). Tail-guard тест использует rtol-style проверку, не точное значение. |
| **`ui.input.forProperty` API может работать иначе** | Fallback на `ui.input.form` с одним property; уточняется при имплементации. |
| **NIST .dat parser может не сработать на edge case'ах форматирования** | Файлы NIST stable с 1997 г.; парсер тестируется при первой регенерации; в репозитории — только результат, не парсер. |
| **Welch's F с дробным df₂ в jStat.centralF.inv** | Поддерживается через ibeta-inverse, проверяется первым же NIST Welch-тестом. |
| **Изменение типа `OneWayAnovaReport` на union** | Все consumers (`addVizualization`, `getAnovaGrid`) обновляются в той же PR; TypeScript ловит несовместимости при сборке. |

---

## Приложение A. Канареечные тесты — что они защищают

Если в будущем кто-то «оптимизирует» код, эти тесты должны упасть:

1. **Welford guard** (`Variance stable on shifted data`) — защищает от отката `setStats` к наивной формуле `Σx² − (Σx)²/n`.
2. **Tail SF guard** (`p-value preserves precision in extreme tail`) — защищает от возврата `1 - jStat.centralF.cdf(...)` в коде.
3. **NIST SmLs07 Fisher** — защищает от потери точности на стиффных данных (mean ≫ std в группах).

Эти три теста — minimum bar для тех, кто рефакторит `anova-tools.ts` в будущем.
