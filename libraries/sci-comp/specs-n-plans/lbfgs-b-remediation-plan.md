# План устранения дефектов реализации L-BFGS-B

Документ описывает упорядоченный набор правок к TypeScript-реализации L-BFGS-B, расположенной в `src/optimization/single-objective/optimizers/lbfgs-b/`, и сопутствующие изменения в тестах. Восемь правок упорядочены от наиболее изолированных к наиболее инвазивным; каждый коммит независимо проходит CI зелёным и независимо ревьюабелен.

Документ опирается на ревью из предыдущего обсуждения. Канонические источники: Byrd, Lu, Nocedal, Zhu (1995, далее BLNZ95); Zhu, Byrd, Lu, Nocedal (1997, далее Zhu1997); Morales, Nocedal (2011); Moré, Thuente (1994); и Fortran-реализация L-BFGS-B v3.0 (`mainlb.f`, `lnsrlb.f`, `cauchy.f`, `bmv.f`, `subsm.f`, MINPACK-2 `dcsrch.f`/`dcstep.f`).

---

## Совместимость публичного API

L-BFGS-B экспортирует наружу через `index.ts` ровно четыре сущности:

- класс `LBFGSB` (extends `Optimizer<LBFGSBSettings>`),
- интерфейсы `LBFGSBSettings`, `LBFGSBLineSearchSettings`,
- тип `LBFGSBBounds`.

Плюс общий контракт `OptimizationResult = {point, value, iterations, converged, costHistory}` и `IterationCallback`.

**Ни одна из этих TS-сигнатур не меняется ни в одном коммите 1–8.** Внутренние имена `BFGSMat`, `runLineSearch`, `subspaceMin`, `cauchyPoint`, `BFGSMat.update`, `SubspaceResult.improved` — не часть публичного API: они не реэкспортированы из `lbfgs-b/index.ts` и недоступны через `@datagrok-libraries/sci-comp`.

### Семантические сдвиги, видимые через `OptimizationResult`

Несмотря на сохранение shape, **семантика поля `converged` ужесточается** в коммитах 6 и 7. Это главное изменение, которое стоит явно держать в голове:

| До правок | После правок |
|---|---|
| `converged === true` мог означать одно из четырёх: T1 (KKT), T2 (Δf), SD-fallback не дал убывания, или `stpMax = 0` | `converged === true` означает строго T1 (`‖proj_grad‖_∞ ≤ gradTolerance`) или T2 (`|Δf|/max(...) ≤ tolerance`) |

Это **усиление** контракта в пользу пользователя: после правок `result.converged === true` действительно означает близость к KKT/стационарности. Но семантически это breaking change: код вида `if (result.converged) trustPoint(result.point)` после правок может видеть `converged === false` в случаях, где раньше получал `true`. Возвращаемые `result.point` и `result.value` остаются best-so-far и осмысленны независимо от `converged` — никакого «было работает / стало не работает» не появляется, просто поле теперь честнее.

### Числовые метрики (без контрактных гарантий)

Это **не** часть формального контракта, но движется и может влиять на тесты пользователя, опирающиеся на конкретные числа:

| Поле | Эффект |
|---|---|
| `iterations`, `costHistory.length` | На bounded-задачах после commit 8 ожидается заметное снижение. На unconstrained — parity с baseline. |
| `value` / `point` | На bounded-задачах после commit 8 — лучше или равно baseline. |
| Число вызовов `onIteration` | Идёт за `iterations`. |

Имена и типы полей `extra` в `IterationCallback` (`projGradInfNorm`, `functionEvaluations`, `stepSize`, `lineSearchSteps`, `historyCount`, `activeBounds`) — сохранены без изменений.

### Изоляция от других оптимизаторов

Все правки серии локализованы внутри `optimizers/lbfgs-b/` и соответствующих тестов в `__tests__/lbfgs-b-*`. Общие модули sci-comp **не трогаются**:

| Модуль | Статус |
|---|---|
| `single-objective/optimizer.ts` (Optimizer base) | не меняется |
| `single-objective/types.ts` (`ObjectiveFunction`, `OptimizationResult`, `IterationCallback`, ...) | не меняется |
| `single-objective/penalty.ts`, `registry.ts` | не меняется |
| `single-objective/optimizers/{nelder-mead,pso,gradient-descent,adam,lbfgs}.ts` | не меняется |
| `multi-objectives/*`, `time-series/*` | не меняется |

Соседние оптимизаторы (Nelder-Mead, PSO, GradientDescent, Adam, L-BFGS) не импортируют ничего из `lbfgs-b/` и поведенчески не затронуты — их собственные `OptimizationResult` на тех же задачах **поразрядно идентичны** baseline после каждого коммита серии.

**Единственный наблюдаемый эффект на «соседей»** — общие сравнительные таблицы в `unconstrained-benchmarks.md`, `multistart-benchmarks.md`, `bounded-benchmarks.md`: столбец L-BFGS-B в них меняется (это и есть цель), а **относительный ранкинг** оптимизаторов на bounded-задачах сдвинется в пользу L-BFGS-B после commit 8. Это не регрессия для других оптимизаторов; это ожидаемая консолидация преимущества L-BFGS-B на задачах с активными bounds.

**Контрольная точка:** при ревью каждого коммита проверить, что в diff не появились правки в любом файле вне `optimizers/lbfgs-b/` и `__tests__/lbfgs-b-*`. Исключения — только: (а) `CLAUDE.md` (архитектурное дерево), (б) `*.md`-отчёты бенчмарков. Любая правка вне этих границ — повод остановиться и пересмотреть план.

---

## Сводная таблица правок

| # | Коммит | Bug | Файлы | ~LoC | TS-сигнатуры (вне `lbfgs-b/`) | `OptimizationResult` поведение |
|---|---|---|---|---|---|---|
| 0 | Baseline-снапшот бенчмарков | — | `*.md` | +reports | — | — |
| 1 | NaN/Inf cap в Moré–Thuente | #7 | `line-search.ts` | +20 | — | при патологическом NaN: было silent success, станет `converged=false` |
| 2 | `α₀` через `‖d‖_2` + условие `boxed` | #5 | `driver.ts` | +15 | — | `iterations`/`functionEvaluations` могут двигаться |
| 3 | Convergence-tests reorder + сброс `retriedMemReset` | #8 | `driver.ts` | ~30 рефактор | — | `iterations` по сошедшимся задачам без изменений; nskip/theta для воспроизводимости с SciPy |
| 4 | Curvature gate через `−gᵀs` | #6 | `bfgs-mat.ts`, `driver.ts` | +10 | — *(`BFGSMat.update` внутренний)* | numerical traces двигаются ради SciPy parity |
| 5 | Snapshot `(x,f,g)` + restore перед reset | #3 | `driver.ts` | +25 | — | без изменений (страховка на recovery-путях) |
| 6 | Снять projected-SD fallback | #4 | `driver.ts` | −20 | — | **`converged` ужесточение**: false-success path #1 закрыт |
| 7 | `stpMax ≤ 0` → recovery вместо converged | #2 | `driver.ts` | +10 | — | **`converged` ужесточение**: false-success path #2 закрыт |
| 8 | Morales–Nocedal 2011: angle-test + truncation fallback | #1 | `subspace.ts` | ~150 переписать | — *(семантика `improved` внутренняя)* | bounded: `iterations` ↓, `value` лучше или равно |

Идея ordering: **0** — baseline-замер до правок, чтобы ловить регрессии. **1–3** — точечные правки внутри одного файла каждая, никакой контракт не меняется. **4** меняет внутреннюю сигнатуру `BFGSMat.update`, но локально и невидимо снаружи. **5–7** — связка в `driver.ts`: 5 ставит snapshot-инфраструктуру; 6 снимает fallback и опирается на restore из 5; 7 переиспользует тот же путь. Коммиты 6 и 7 — единственные, ужесточающие семантику `converged`. **8** — самая большая правка, изолирована в `subspace.ts`, контракт `ws.xHat` сохраняется, поэтому driver не трогается.

Порядок 5→6 (snapshot до снятия fallback) — намеренный: первый раз делалась обратная последовательность, что приводило к forward-refs на ещё не существующие буферы. Объединять 5 и 6 в один коммит можно, но текущий порядок даёт более чистое ревью.

### Ритуал каждого коммита

В **каждом** коммите серии — следующий чеклист, без исключений:

1. Применить правку кода + связанные правки тестов **в одном коммите** (атомарность).
2. Прогнать `npm run lint-fix && npm run build && npm test` — должно быть зелёным.
3. Если коммит затрагивает архитектуру (4, 8 — точно; 3, 5 — возможно) — обновить дерево в `libraries/sci-comp/CLAUDE.md`.
4. Если коммит меняет публичный API (4: `BFGSMat.update`; 8: семантика `improved`) — обновить TSDoc и README в той же папке.

---

## Анализ существующих тестов

| Файл | Сценарий | Действие | Причина |
|---|---|---|---|
| `lbfgs-b-bounds.test.ts` | все | не трогать | Вне scope правок |
| `lbfgs-b-cauchy.test.ts` | все | не трогать | Cauchy-логика корректна |
| `registry.test.ts` | все | не трогать | Не задет |
| `lbfgs-b.test.ts` (scaffold) | все | не трогать | Только validation |
| `lbfgs-b-bfgs-mat.test.ts` | curvature-gate тесты | переписать | Тест проверяет `‖y‖²`-формулу — после commit 4 формула другая |
| `lbfgs-b-bfgs-mat.test.ts` | остальные `mat.update(...)` | обновить вызовы | Добавить 4-й параметр `negGTs` или сохранить `curvatureEps = 0` |
| `lbfgs-b-subspace.test.ts` | "Morales–Nocedal backtrack path" | переписать на 2 теста | Старый тест проверяет неправильную схему |
| `lbfgs-b-subspace.test.ts` | "model decrease", "Newton-like", "feasible" | оставить, переименовать | В unconstrained случае поведение не меняется |
| `lbfgs-b-line-search.test.ts` | NaN-bisection тесты | оставить | Существующих ≤6 бисекций — под NAN_CAP |
| `lbfgs-b-integration.test.ts` | все | оставить | Толерантности достаточно loose |

### Новые тесты к добавлению

| # | Тест | Файл | Доказывает |
|---|---|---|---|
| N1 | NaN-cap exhaustion → `status='error'` | `lbfgs-b-line-search.test.ts` | Bug #7 |
| N2 | Curvature gate с конкретным `−gᵀs` | `lbfgs-b-bfgs-mat.test.ts` | Bug #6 (новая формула) |
| N3 | Subspace: feasible projected → angle-test path | `lbfgs-b-subspace.test.ts` | Bug #1, ветка успеха |
| N4 | Subspace: bad-angle projected → truncation fallback | `lbfgs-b-subspace.test.ts` | Bug #1, ветка fallback |
| N5 | Driver: после reset `(x, f, g)` восстановлены | `lbfgs-b-integration.test.ts` | Bug #3 |
| N6 | Driver: `stpMax = 0` не ведёт к converged=true | `lbfgs-b-integration.test.ts` | Bug #2 |

N5 и N6 сложно сделать чистым unit-тестом без mock-ов; реалистичнее через интеграционную задачу с искусственно подобранной целевой функцией, которая на ранней итерации выдаёт ситуацию `slope ≥ 0` или `stpMax = 0`, а после reset находит верный путь. Конкретный setup для N6 приведён в §Commit 7.

---

## Commit 0 — Baseline-снапшот бенчмарков

### Цель

Зафиксировать численные характеристики L-BFGS-B **до** правок, чтобы на каждом последующем коммите ловить регрессии (число итераций, nfev, success rate, feasibility violation). Без baseline regression detection невозможен.

### Действия

Из ветки до правок прогнать:

```bash
npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts > /tmp/baseline-unconstrained.md
npx tsx src/optimization/single-objective/benchmarks/multistart-benchmarks.ts  > /tmp/baseline-multistart.md
npx tsx src/optimization/single-objective/benchmarks/bounded-benchmarks.ts     > /tmp/baseline-bounded.md
```

Сохранить три `.md` в `specs-n-plans/lbfgs-b-baseline/` (новая директория) или приложить к PR описанию. После каждого коммита 1–8 сравнивать только релевантные метрики из соответствующего runner-а:

| Коммит | Чем проверяется |
|---|---|
| 1 | unconstrained (NaN-cap не должен влиять на «здоровые» задачи) |
| 2 | unconstrained, multistart (другой `α₀` → другой trace) |
| 3 | все три (последний BFGS update пропущен → −1 итерация на сошедшихся задачах) |
| 4 | все три (curvature gate перерисовывает `nskip`) |
| 5 | все три (snapshot: numerical parity ожидается; функциональных изменений нет) |
| 6 | bounded (снят SD fallback → ранее «успешные» проходы могут стать неуспешными — это сигнал о замаскированном баге, **не** регрессия) |
| 7 | bounded (false success → теперь corrected) |
| 8 | bounded — главная цель; ожидается значительное улучшение |

### Тесты

Не применимо. Это документационный коммит.

---

## Commit 1 — NaN/Inf cap в `line-search.ts`

### Bug

NaN/Inf bisection в `runLineSearch`/`runLineSearchAsync` не имеет отдельного счётчика. При расходящейся целевой функции цикл бисектит к `state.stx` без верхней границы (кроме общего `nfev ≥ maxSteps`), и когда `stx = 0` это даёт «тихую сходимость» в стартовой точке — оптимизатор сообщает success в плохой точке. Дополнительно, текущая реализация может обновлять `state.stx/sty/fx/fy/gx/gy` в ходе обработки NaN-точки, что разрушает инвариант брекета Moré–Thuente.

### Изменения

В `line-search.ts`:

```ts
const NAN_CAP = 20;

// в runLineSearch:
let stp = alpha0;
let phi = phi0;
let phiPrime = phiPrime0;
let isFirst = true;
let nfev = 0;
let nanCount = 0;     // ← новый счётчик

for (let k = 0; k <= params.maxSteps; k++) {
  const tick = dcsrch(state, stp, phi, phiPrime, isFirst, /* … */);
  isFirst = false;
  stp = tick.stp;

  if (tick.task !== 'fg')
    return finalise(tick.task, stp, phi, phiPrime, nfev, evalFn);

  if (nfev >= params.maxSteps)
    return finalise('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);

  let ev = evalFn(stp);
  nfev++;
  while (!Number.isFinite(ev.phi) || !Number.isFinite(ev.phiPrime)) {
    nanCount++;
    if (nanCount >= NAN_CAP)
      return finalise('error', state.stx, state.fx, state.gx, nfev, evalFn);
    if (nfev >= params.maxSteps)
      return finalise('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);
    const shrunk = state.stx + 0.5 * (stp - state.stx);
    if (shrunk === stp || shrunk <= stpmin)
      return finalise('warning_rounding', state.stx, state.fx, state.gx, nfev, evalFn);
    stp = shrunk;
    ev = evalFn(stp);
    nfev++;
  }
  // в обычном пути обнуляем nanCount: подряд считаются только последовательные NaN
  nanCount = 0;
  phi = ev.phi;
  phiPrime = ev.phiPrime;
}
```

Аналогично в `runLineSearchAsync`.

Важно: счётчик накапливает NaN-точки **в пределах одного NaN-bisection эпизода**. На следующем `dcsrch`-tick (когда evaluator вернул конечное значение) `nanCount` обнуляется. Это семантически совпадает с Fortran-овским `iback`: `iback` накапливается до 20 *в пределах одного line search*, потом следует `RESTART_FROM_LNSRCH`.

`stx/sty/fx/fy/gx/gy` уже не обновляются на NaN-точке в текущем коде (потому что `phi/phiPrime` обновляются после while-цикла), так что инвариант брекета сохраняется автоматически.

### Тесты

**N1 — добавить:**

```ts
it('returns error after NAN_CAP consecutive non-finite evaluations', () => {
  // Целевая функция всегда возвращает NaN — бисекция не помогает.
  const phi: PhiEval = () => ({phi: NaN, phiPrime: NaN});
  const r = runLineSearch(phi, 0, -1, 1, DEFAULT_PARAMS);
  expect(r.status).toBe('error');
  expect(r.ok).toBe(false);
  expect(r.nfev).toBeLessThan(50);  // должен выйти за O(NAN_CAP), не за maxSteps
});
```

**Существующие NaN-тесты** (`'handles NaN from evaluator by bisecting toward stx'` и `'handles persistently bad region via repeated bisection'`) **не трогать**: они выполняют ≤6 бисекций каждый и под NAN_CAP проходят без изменений.

---

## Commit 2 — `α₀` через 2-норму и условие `boxed`

### Bug

В `driver.ts` начальный шаг вычисляется как `min(1, 1/max(‖d‖_∞, 1e-16))` и применяется при `iteration === 0 || retriedMemReset`. Fortran `lnsrlb.f` (строки 152–158 v3.0) использует **евклидову норму**: `dnorm = sqrt(ddot(n,d,1,d,1)); stp = min(one/dnorm, stpmx)`, и только при `iter == 0 .and. .not. boxed`. Различие даёт начальный шаг до фактора `√n` больше; на больших задачах это приводит к лишним backtrack-ам в Moré–Thuente.

Дополнительно: после `retriedMemReset` Fortran не сбрасывает `iter` и поэтому идёт в ветку `else stp = one`. TS-код повторно применяет эвристику, что излишне осторожно.

### Изменения

В начале `runSync` / `runAsync`, после `normalizeBounds`:

```ts
const {lower, upper, nbd, anyFinite} = normalizeBounds(s.bounds, n);

let boxed = anyFinite;
if (boxed) {
  for (let i = 0; i < n; i++) {
    if (nbd[i] !== BOUND_BOTH) { boxed = false; break; }
  }
}
```

`boxed` соответствует Fortran-флагу: `true` тогда и только тогда, когда **каждая** переменная имеет конечные оба ограничения.

Заменить блок выбора `alpha0`:

```ts
// было:
let alpha0 = iteration === 0 || retriedMemReset ?
  Math.min(1, 1 / Math.max(infNorm(d), 1e-16)) :
  1;
if (alpha0 > stpMax) alpha0 = stpMax;

// стало:
let alpha0: number;
if (iteration === 0 && !boxed) {
  const dNorm = Math.sqrt(dot(d, d));
  alpha0 = Math.min(1 / Math.max(dNorm, 1e-16), stpMax);
} else {
  alpha0 = Math.min(1, stpMax);
}
```

Удалить вспомогательную функцию `infNorm`, если больше не используется.

### Тесты

Новых не нужно. Интеграционные тесты `lbfgs-b-integration.test.ts` будут регрессировать через `iterations`/`functionEvaluations`, если правка введена с ошибкой; явно прописывать численные ожидания на эти счётчики не стоит — они FP-неассоциативны.

---

## Commit 3 — порядок тестов сходимости и сброс `retriedMemReset`

### Bug

В текущем драйвере порядок такой:

```
mat.update(stepS, stepY, ε)   ← BFGS update
x.set(xNew); g.set(gNew); fCur = fNew
if (fChange ≤ ftol) converged = true; break
pgNorm = projectedGradient(...)
if (pgNorm ≤ gradTol) converged = true; break
```

Fortran выполняет оба теста **до** `matupd`. Это лишний BFGS-update на финальной итерации; влияет на воспроизводимость относительно SciPy через расхождение `nskip` и `theta`. Дополнительно: `retriedMemReset` ставится в `true` после неудачи line search и **не сбрасывается** после успешной итерации, что ограничивает recovery до одной попытки за всё время выполнения.

Также отсутствует pre-loop pgnorm-check для случая, когда стартовая точка уже оптимальна (Fortran выходит до первой итерации).

### Изменения

Pre-loop check уже частично есть:

```ts
let pgNorm = projectedGradient(x, g, lower, upper, nbd, pg);
if (pgNorm <= gradTol)
  return finaliseResult(bestX, bestValue, 0, true, costHistory, costLen);
```

Это корректно. Оставить.

Перепорядочить хвост итерации:

```ts
// (h) accept step
const stp = lsResult.stp;
for (let i = 0; i < n; i++) xNew[i] = x[i] + stp * d[i];
project(xNew, lower, upper, nbd, xNew);
const fNew = lsResult.phi;

// (i) Δf relative test — ДО BFGS update
const denom = Math.max(Math.abs(fCur), Math.abs(fNew), 1);
const fChange = Math.abs(fCur - fNew) / denom;

// (m) projected-gradient test — тоже ДО BFGS update,
//     требует gNew, который у нас уже есть из last LS evaluation
const pgNormNew = projectedGradient(xNew, gNew, lower, upper, nbd, pg);

if (fChange <= ftol || pgNormNew <= gradTol) {
  // принять шаг и выйти, но без BFGS update — он уже не понадобится
  x.set(xNew); fCur = fNew; g.set(gNew);
  if (fCur < bestValue) { bestValue = fCur; bestX.set(x); }
  iteration++;
  costHistory[costLen++] = bestValue;
  converged = true;
  break;
}

// (j, k) BFGS update
for (let i = 0; i < n; i++) {
  stepS[i] = xNew[i] - x[i];
  stepY[i] = gNew[i] - g[i];
}
mat.update(stepS, stepY, CURVATURE_EPS);  // в commit 4 здесь добавится 4-й аргумент

// (l) advance state
x.set(xNew); fCur = fNew; g.set(gNew);
if (fCur < bestValue) { bestValue = fCur; bestX.set(x); }
lastStepSize = stp;
lineSearchSteps = lsResult.nfev;
iteration++;
costHistory[costLen++] = bestValue;

retriedMemReset = false;   // ← сброс после успешной итерации

// (n) eval limit
if (fEvalCount >= maxFEval) break;

// (o) callback
if (fireCallback(onIter, iteration, bestValue, bestX, {
  projGradInfNorm: pgNormNew,
  // …
})) break;
```

Аналогично для `runAsync`.

Замечание: `pgNormNew` теперь вычисляется один раз и переиспользуется в callback. Переменная `pgNorm` (на старом `g`) больше не нужна в теле цикла — её область видимости сжимается до pre-loop check.

Замечание про буфер `pg`: scratch-массив `pg` сейчас передаётся в `projectedGradient(x, g, lower, upper, nbd, pg)` как out-параметр и используется только внутри той же функции для оценки `‖proj_grad‖_∞`. После reorder последовательность такая: `projectedGradient(...) → pgNormNew → (return / break)` или `projectedGradient(...) → pgNormNew → BFGS update → callback (читает pgNormNew, не pg)`. Между записью в `pg` и его последним consumer-ом никто не пишет в этот буфер — инвариант сохраняется. Если в будущем добавится consumer, использующий сам массив `pg` (а не скаляр `pgNormNew`), эту страховку нужно пересмотреть.

### Тесты

Существующие интеграционные тесты могут показать **уменьшение** числа итераций на 1 (последний BFGS update не делается). Если тесты содержат жёсткие верхние границы вроде `expect(r.iterations).toBeLessThanOrEqual(20)` — они продолжат проходить. Если есть `toBe(20)` — нужно ослабить, но в текущих файлах таких ожиданий нет.

Тесты `'callback receives populated extra'` ожидают `extra.projGradInfNorm` присутствие; контракт сохраняется. `'stops iteration when callback returns true'` — поведение не меняется, callback по-прежнему вызывается после iteration++.

---

## Commit 4 — curvature gate через `−gᵀs`

### Bug

В `BFGSMat.update`:

```ts
if (!(yy > 0) || !(ys > curvatureEps * Math.max(1, yy)))
  return false;
```

Это вариант BLNZ95 (3.9): `sᵀy > eps · ‖y‖²`. Fortran v3.0 использует Zhu1997 (3): `sᵀy > eps · |gᵀs|`, конкретно в `mainlb.f`:

```fortran
dr = (gd - gdold) * stp        ! = sᵀy
ddum = -gdold * stp            ! = -gᵀs > 0
if (dr .le. epsmch * ddum) then
   nskip = nskip + 1; updatd = .false.; goto 888
endif
```

На плохо отмасштабированных задачах эти два теста расходятся, теряется numerical parity со SciPy.

### Изменения

`BFGSMat.update` получает 4-й параметр:

```ts
update(
  sNew: Float64Array,
  yNew: Float64Array,
  curvatureEps: number,
  negGTs: number,        // = -gᵀs · stp = -gᵀd · stp = -slope · stp
): boolean {
  const n = this.n;
  const m = this.m;

  if (sNew.length !== n || yNew.length !== n)
    throw new Error('BFGSMat.update: sNew/yNew length mismatch');

  let ys = 0;
  let yy = 0;
  let ssNew = 0;
  for (let i = 0; i < n; i++) {
    const si = sNew[i];
    const yi = yNew[i];
    ys += si * yi;
    yy += yi * yi;
    ssNew += si * si;
  }
  if (!(yy > 0) || !(ys > curvatureEps * Math.max(0, negGTs)))
    return false;
  // ... остальное без изменений
}
```

Семантика: при `curvatureEps = 0` тест становится `ys > 0` независимо от `negGTs`. Это сохраняет совместимость с тестами, которые передают `0`.

В `driver.ts` на месте `mat.update(...)`:

```ts
mat.update(stepS, stepY, CURVATURE_EPS, -slope * stp);
```

`slope` — переменная, уже вычисленная как `dot(g, d)` ранее в той же итерации. `-slope * stp = -gᵀd · stp = -gᵀ(x_{k+1} - x_k)`, что в точности соответствует Fortran-вычислению `ddum`.

### Тесты

`lbfgs-b-bfgs-mat.test.ts`, секция «Curvature gate» — переписать.

**Старый тест** `'rejects pair with sᵀy too small relative to yᵀy'` удалить.

**Новые тесты:**

```ts
describe('BFGSMat curvature gate', () => {
  it('rejects pair with sᵀy ≤ 0', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1, 0, 0]);
    const y = new Float64Array([-1, 0, 0]);
    expect(mat.update(s, y, 1e-12, 1)).toBe(false);
    expect(mat.col).toBe(0);
  });

  // N2:
  it('rejects pair when sᵀy ≤ eps · (−gᵀs)', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1e-10, 0, 0]);
    const y = new Float64Array([1, 0, 0]);
    // sᵀy = 1e-10, negGTs = 1, eps = 1e-8 → threshold = 1e-8 → reject.
    expect(mat.update(s, y, 1e-8, 1)).toBe(false);
  });

  it('accepts pair well above threshold', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1, 0, 0]);
    const y = new Float64Array([1, 0, 0]);
    // sᵀy = 1, negGTs = 0.5, eps = 1e-8 → threshold = 5e-9 → accept.
    expect(mat.update(s, y, 1e-8, 0.5)).toBe(true);
  });

  it('with curvatureEps = 0 accepts any positive sᵀy', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1e-20, 0, 0]);
    const y = new Float64Array([1, 0, 0]);
    expect(mat.update(s, y, 0, 0)).toBe(true);
  });
});
```

**Остальные тесты в `lbfgs-b-bfgs-mat.test.ts`**, передающие `mat.update(s, y, 0)`, — обновить вызов до `mat.update(s, y, 0, 0)` (или `mat.update(s, y, 0, 1)` — при `curvatureEps = 0` четвёртый параметр игнорируется). Это все тесты в секциях «Construction», «single-pair compact form», «Secant equation», «solveM inverse check», «applyW / applyWt», «Ring-buffer semantics».

---

## Commit 5 — snapshot и restore перед reset

### Bug

В `driver.ts` есть точки, где после неудачи делается `mat.reset(); continue`, но `(x, f, g)` к этому моменту уже могут быть мутированы внутри line search или приёма шага. Fortran перед `RESTART_FROM_LNSRCH` явно делает:

```fortran
call dcopy(n, t, 1, x, 1)    ! восстановить старый x
call dcopy(n, r, 1, g, 1)    ! восстановить старый g
f = fold                     ! восстановить старое f
```

Без этого retry стартует с уже частично-сдвинутой точки, что несовместимо с инвариантом «restart along steepest descent **from the same x_k**» из Zhu1997 §4.

### Изменения

Добавить буферы в начало `runSync`/`runAsync` рядом с другими preallocated arrays:

```ts
const xSave = new Float64Array(n);
const gSave = new Float64Array(n);
let fSave = 0;
```

В начало каждой итерации главного цикла, перед Cauchy-фазой:

```ts
while (iteration < maxIter) {
  if (fEvalCount >= maxFEval) break;

  // Snapshot for memory-reset retry.
  xSave.set(x);
  gSave.set(g);
  fSave = fCur;

  // ---- (a) Cauchy ------------------------------------------------
  const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchy);
  // ...
}
```

Перед каждой точкой `mat.reset(); retriedMemReset = true; continue` — добавить restore. Точек три:

1. Cauchy fail (`!cr.ok`):
   ```ts
   if (!cr.ok) {
     if (mat.col === 0) break;
     x.set(xSave); g.set(gSave); fCur = fSave;
     mat.reset();
     retriedMemReset = true;
     continue;
   }
   ```

2. Line search fail (`!lsResult.ok`):
   ```ts
   if (!lsResult.ok) {
     if (retriedMemReset || mat.col === 0) break;
     x.set(xSave); g.set(gSave); fCur = fSave;
     mat.reset();
     retriedMemReset = true;
     continue;
   }
   ```

3. Slope check (после commit 6): аналогично.

4. (После commit 7) `stpMax ≤ 0`: аналогично.

В рамках **этого** коммита (5) добавляются только snapshot-буферы и блок их инициализации в начале итерации, плюс restore в точках (1) и (2). Точки (3) и (4) приходят с коммитами 6 и 7 соответственно — там лишь дописывается одна строка с restore рядом с уже существующим `mat.reset()`.

В точке (1) — Cauchy fail — `(x, f, g)` ещё **не** мутированы (Cauchy не пишет в `x`, только в `cauchy.xc`), поэтому restore технически избыточен. Но добавить его всё равно для единообразия и чтобы не зависеть от деталей реализации Cauchy.

Замечание: `xNew`, `gNew` буферы line search — это отдельные scratch-массивы, не `x`, `g`. Текущий код мутирует `x` только в строке `x.set(xNew)` после успешного шага. Поэтому в точке (2) — Line search fail — `x`/`g` тоже формально не изменены... **но** в TS line search через замыкание получает `x` и пишет в `xNew = x + α·d`, а `gNew` через `evalGrad(xNew, gNew)`. Никакого побочного изменения `x`/`g` быть не должно. Однако: внутри `evalGrad`, когда используется finite-differences fallback, временно мутируется `p`. И `p` в FD-вызове — это `xNew`, не `x`. Значит `x` действительно неизменно.

Несмотря на это, snapshot — дешёвая страховка: 2 × n doubles на итерацию. На 1000-мерной задаче и 1000 итерациях это 16 МБ выделений (если бы выделялись каждый раз), но мы предварительно аллоцировали — это два `set(x)` на итерацию, O(n) per iteration, что ничтожно по сравнению с O(mn) самой итерации.

### Тесты

**N5 — добавить:** интеграционный тест, в котором первая итерация заведомо запускает recovery, и проверяется, что после recovery точка корректна.

```ts
it('recovers correctly after line-search failure with snapshot restore', () => {
  // Целевая функция, которая на первой итерации провоцирует line search failure
  // (например, очень крутая в одном направлении), но затем ведёт к минимуму.
  // Ожидаем converged=true и x близко к минимуму.
  const opt = new LBFGSB();
  // … подбираем сложную функцию
});
```

Конкретный choice функции потребует подбора; альтернатива — unit-тест с мокнутым `runLineSearch`, возвращающим `{ok: false}` на первой итерации. Это вариант проще, но требует exposed-точки для замены; делать только если интеграционный тест не получится скомпоновать.

---

## Commit 6 — снять projected-SD fallback

### Bug

В `driver.ts` после вычисления `slope = dot(g, d)`:

```ts
if (!(slope < 0)) {
  // Fall back to projected steepest descent.
  for (let i = 0; i < n; i++) {
    let di = -g[i];
    if ((nbd[i] === BOUND_LOWER || nbd[i] === BOUND_BOTH) && x[i] <= lower[i] && di < 0) di = 0;
    if ((nbd[i] === BOUND_UPPER || nbd[i] === BOUND_BOTH) && x[i] >= upper[i] && di > 0) di = 0;
    d[i] = di;
  }
  slope = dot(g, d);
  if (!(slope < 0)) {
    converged = true;
    break;
  }
}
```

В Fortran `lnsrlb.f` ситуация `gd ≥ 0` обрабатывается по-другому:

```fortran
if (gd .ge. zero) then
   write(6,*)' ascent direction in projection gd = ', gd
   info = -4
   return
endif
```

То есть line search возвращает ошибку `info = -4`, и в `mainlb` через `info != 0` идёт `RESTART_FROM_LNSRCH` (с reset памяти), либо `ABNORMAL_TERMINATION_IN_LNSRCH`. Никакого fallback на projected SD в каноне нет.

Опасность TS-фоллбэка не в его математике (после `mat.reset()` Cauchy point по построению лежит на проекционном градиенте, поэтому конечный путь эквивалентен), а в **маскировке других багов**: знаковая ошибка в `solveM`, неверный знак в Cauchy `f'`-обновлении, неверная передача `θ` в `applyW` — всё это будет проявляться как `slope ≥ 0`, и фоллбэк тихо превратит метод в проекционный градиент, оставив видимость работоспособности.

### Изменения

Заменить блок fallback на (snapshot-буферы из commit 5 уже доступны):

```ts
if (!(slope < 0)) {
  if (retriedMemReset || mat.col === 0) {
    converged = false;
    break;
  }
  x.set(xSave); g.set(gSave); fCur = fSave;
  mat.reset();
  retriedMemReset = true;
  continue;
}
```

К этой точке итерации `(x, f, g)` ещё не мутированы (line search ещё не запускался), поэтому restore технически избыточен. Но он добавлен для единообразия с другими recovery-точками и устойчивости к будущим рефакторингам.

### Тесты

Новых не нужно. Существующие интеграционные тесты не должны были опираться на эту ветку (если опираются — это и есть тот замаскированный баг, который мы хотим обнаружить).

**Watchpoint:** в bounded-Rosenbrock-тесте после этого коммита может всплыть скрытый баг в Cauchy/subspace, ранее замаскированный fallback. Если тест краснеет — это **не** регрессия commit 6 как такового; это сигнал, что один из багов из remediation-серии (вероятно #1, исправляется в commit 8) уже играл роль и до этого. В таком случае оставлять тест красным до commit 8 нельзя — нужно дополнительно временно ослабить ассерты или пометить как `xit` с TODO-комментарием, ссылающимся на commit 8.

---

## Commit 7 — `stpMax ≤ 0` → recovery

### Bug

```ts
const stpMax = maxFeasibleStep(x, d, lower, upper, nbd);
if (!(stpMax > 0)) {
  converged = true;
  break;
}
```

Это самый опасный из найденных багов: при нулевом `stpMax` оптимизатор сообщает success в произвольной точке, даже если `pgNorm ≫ gradTol`. Fortran трактует ту же ситуацию через line search: если `stpmx = 0`, то `dcsrch` сразу возвращает ошибку (`stp ≤ stpmin`), и через `iback ≥ 20` идёт `RESTART_FROM_LNSRCH`.

### Изменения

```ts
const stpMax = maxFeasibleStep(x, d, lower, upper, nbd);
if (!(stpMax > 0)) {
  if (retriedMemReset || mat.col === 0) {
    converged = false;
    break;
  }
  x.set(xSave); g.set(gSave); fCur = fSave;
  mat.reset();
  retriedMemReset = true;
  continue;
}
```

После reset Cauchy point по построению лежит на проекционном градиенте, направление `d = xc − x` обычно не утыкается в границу, и `stpMax > 0`. Если и после reset `stpMax = 0` — это вырожденная задача (например, `x` в углу и градиент направлен в угол), и выход с `converged = false` корректен.

### Тесты

**N6 — добавить:** интеграционный тест с задачей, где на первой итерации `stpMax = 0` сохраняется после reset, и проверяется, что результат — `converged = false`, без ложного «успеха». В рамках текущего scope это **наиболее надёжный** вариант: он не зависит от того, найдёт ли Cauchy/subspace допустимое направление после `mat.reset()`, и проверяет именно главное обещание commit 7 — отсутствие false success.

Концепт целевой функции:

```ts
it('does not report converged when stpMax = 0 at unsolved point', () => {
  // Двумерная линейная задача f(x) = -x[0] - x[1] с границами [0, 1]^2,
  // старт в углу x₀ = (1, 1) (верхняя граница для обеих координат).
  //
  // - Истинный минимум также в углу (1, 1) — функция на единичном квадрате
  //   достигает минимума -2 именно там. Но pgnorm в (1, 1) = 0 (градиент (-1,-1)
  //   "толкает" вовне границы → проекция нулевая) → pre-loop check сразу
  //   вернёт converged. Это НЕ то, что нужно.
  //
  // Поэтому используем f(x) = x[0] + x[1] (минимум в (0, 0), но старт (1,1)),
  // ИЛИ запретим pre-loop convergence, выбрав границы так, чтобы на старте
  // не все компоненты были active с правильным знаком.
  //
  // Рабочий вариант: f(x) = -x[0] - x[1] с границами [-1, 1]^2 и x₀ = (1, 1).
  // - g(x₀) = (-1, -1), pgnorm в (1,1) на границе [-1,1] с верхней: проекция
  //   градиента отрицательна, значит pgnorm > 0, pre-loop НЕ выходит.
  // - Cauchy: ищет точку вдоль -g = (1, 1), но обе координаты уже на верхней
  //   границе → Cauchy point = x = (1, 1).
  // - Free set пуст → endpoint = xc = x → d = 0 → stpMax = 0.
  // - До commit 7: converged = true (false success в немимимальной точке (1,1),
  //   при том что истинный минимум этой задачи на квадрате [-1,1]^2 — в (1,1),
  //   так что example нужно ещё подкрутить).
  //
  // Корректный setup: f(x) = -x[0] + x[1]^2, границы x[0] ∈ [-1, 1],
  // x[1] ∈ [-1, 1], старт x₀ = (1, 0.5). Градиент = (-1, 2·0.5) = (-1, 1).
  // x[0] на верхней границе, проекция -g[0] обнуляется → free set = {1}.
  // Здесь будет нормальный шаг, не stpMax=0. — Этот вариант тоже не подходит.
  //
  // Конкретный setup откладывается до момента, когда commit 7 готов и можно
  // отлаживать на конкретных вызовах. Гарантия: воспроизвести stpMax = 0 в
  // точке с pgnorm > gradTol реально, но требует тонкой подгонки. Если не
  // получится подобрать целевую функцию за разумное время — fallback на
  // mock через exposed entry point (см. ниже).
});
```

**Fallback-вариант:** если интеграционный тест не получается скомпоновать, экспортировать `runSync` из модуля и собрать тест с искусственно подобранным `BFGSMat` + начальным состоянием, гарантированно дающим `stpMax = 0` на первой итерации. Это требует добавления небольшого `__test__`-only API (например, `runSyncWithState`), но компромисс приемлемый: bug #2 — самый опасный (false success), и тест на него **обязан** быть детерминированным, а не «работает у меня на конкретной задаче».

**Минимальная фиксация без рабочего N6:** даже если N6 откладывается, в этом коммите **нужно** добавить хотя бы **юнит-тест на отказ pre-loop convergence**: оптимизатор в неоптимальной точке с искусственно подобранной начальной конфигурацией не должен возвращать `converged = true` без доказательств. Это страховка против будущих регрессий аналогичного класса.

---

## Commit 8 — Morales–Nocedal 2011

### Bug

Текущий `subspace.ts` реализует «project + backtrack от `x_c` с тестом убывания модели», что не совпадает с патчем Morales–Nocedal 2011. Канонический патч (`subsm.f` v3.0, маркеры `c-jlm-jn`) выполняет:

1. Найти безусловный subspace minimizer `x̂ = x_c + Z · d̄ᵘ` (без изменений в текущем коде, шаги (1)–(6)).
2. Проектировать **один раз**: `x̄ = P_{[l,u]}(x̂)`.
3. Сформировать **направление от исходного итерата `x_k`**: `d_k = x̄ − x_k`.
4. Проверить **strong-descent / angle test** на истинном градиенте: `g_k^T d_k ≤ −η · ‖d_k‖ · ‖g_k‖`.
5. Если тест прошёл — endpoint = `x̄` (вернуть его в `ws.xHat`); никакого внутреннего бэктрекинга, line search в драйвере сделает свою работу позже.
6. Если не прошёл — **возврат к схеме 1997**: усечение по линии `x_k → x̂` с максимальным `α ∈ (0, 1]`, при котором `x_k + α(x̂ − x_k)` остаётся feasible.

Текущая реализация неправильна по всем четырём измерениям (опорная точка, критерий, бэктрекинг, fallback), что превращает subspace-фазу в no-op на сложных конфигурациях и оптимизатор вырождается в чистый Cauchy/проекционный градиент.

### Изменения

Шаги (1)–(6) в `subspaceMin` сохраняются дословно. Меняется заключительная часть (шаги (7)–(8) текущего кода).

Замена начиная с момента, где вычислен `du`:

```ts
// (7) Form the unconstrained Newton-like point xHatFull = xc + du.
//     On non-free coords du is zero, so xHatFull[i] = xc[i] there.
const xHatFull = ws.scratchN;     // переиспользуем буфер
for (let i = 0; i < n; i++) xHatFull[i] = xc[i] + du[i];

// (8) Project once onto [l, u].
const xBar = ws.delta;            // переиспользуем буфер размера n
project(xHatFull, lower, upper, nbd, xBar);

// (9) Form direction d_k = xBar - x_k.
//     Compute angle test: g · d ≤ -η · ‖d‖ · ‖g‖.
let gd = 0;
let dd = 0;
let gg = 0;
for (let i = 0; i < n; i++) {
  const di = xBar[i] - x[i];
  ws.bDelta[i] = di;              // храним d_k для (10) если понадобится
  gd += g[i] * di;
  dd += di * di;
  gg += g[i] * g[i];
}
const ETA = 1e-2;                 // порог из subsm.f; при необходимости калибровать
const dNorm = Math.sqrt(dd);
const gNorm = Math.sqrt(gg);

if (gd <= -ETA * dNorm * gNorm) {
  // Angle test passed → endpoint = xBar.
  for (let i = 0; i < n; i++) ws.xHat[i] = xBar[i];
  return {improved: true};
}

// (10) Fallback: 1997 truncation along x_k → xHatFull.
//      Compute α = max{α ∈ (0, 1] : x_k + α (xHatFull - x_k) feasible}.
let alpha = 1;
for (let i = 0; i < n; i++) {
  const stepDir = xHatFull[i] - x[i];
  if (stepDir > 0 && (nbd[i] === BOUND_UPPER || nbd[i] === BOUND_BOTH)) {
    const cap = (upper[i] - x[i]) / stepDir;
    if (cap < alpha) alpha = Math.max(0, cap);
  } else if (stepDir < 0 && (nbd[i] === BOUND_LOWER || nbd[i] === BOUND_BOTH)) {
    const cap = (lower[i] - x[i]) / stepDir;
    if (cap < alpha) alpha = Math.max(0, cap);
  }
}

if (alpha > 0) {
  for (let i = 0; i < n; i++)
    ws.xHat[i] = x[i] + alpha * (xHatFull[i] - x[i]);
  return {improved: true};
}

// (11) Truncation gave α = 0 → fall back to xc.
for (let i = 0; i < n; i++) ws.xHat[i] = xc[i];
return {improved: false};
```

Важные моменты:

- **Контракт `ws.xHat` сохраняется**: всегда содержит endpoint, который драйвер использует как `endpoint - x` для построения `d`. Driver не меняется.
- **`improved`** теперь означает «`xHat` отличается от `xc` по нетривиальной причине». В случае `alpha = 0` — fallback на `xc`, как в текущем коде.
- **Буферы `xHatFull`, `xBar`**: используем существующие scratch-поля workspace (`scratchN`, `delta` или `bz`). Если их не хватает — добавить два n-буфера в `SubspaceWorkspace`.
- **`η`** взят как `1e-2` по `subsm.f`. Калибровка по `bounded-benchmarks.md` (см. правила ниже), т. к. SciPy-parity infrastructure в текущем suite нет.
- **`project(xHatFull, ...)` пишет в `xBar`**: вызов `project(x, l, u, nbd, out)` уже есть в `bounds.ts` и не предполагает алиасинга, поэтому отдельный буфер обязателен.

#### Зачистка устаревшего кода

Удалить из `subspaceMin`:

- весь блок `// (8) Morales–Nocedal backtrack` (цикл по `lambda = 1, 0.5, 0.25, ...`);
- вызов `modelDelta` и подготовку `h = g + B·z` для него;
- параметр `maxBacktrack` функции (унести из сигнатуры).

Из `SubspaceWorkspace`:

- проверить, какие из полей `bz`, `h` после правок ещё используются. `modelDelta` больше не вызывается — поля для него можно удалить или переиспользовать под `xHatFull`/`xBar`.

Функцию `modelDelta` (helper в том же файле) — удалить целиком.

### Тесты

`lbfgs-b-subspace.test.ts`:

#### Существующие тесты

| Тест | Действие |
|---|---|
| `'col=0 → xHat = xᶜ, improved=false'` | оставить |
| `'freeCount=0 → xHat = xᶜ, improved=false'` | оставить |
| `'produces m_k(xHat) < m_k(xᶜ) on unconstrained problem'` | переименовать → `'unconstrained: xHat is the projected Newton step'`. Логика теста (Newton-step property) корректна для unconstrained случая, в котором проекция тождественна и angle test заведомо проходит. |
| `'xHat stays feasible even after backtrack'` | переименовать → `'xHat stays feasible after projection or truncation'`. |
| `'B·(xHat−xᶜ) ≈ −r̄ᶜ for t=n'` | оставить — справедливо в unconstrained случае. |
| `'projects and backtracks when unconstrained step is infeasible'` | **полностью переписать** на N3 + N4. |

#### Новые тесты

**N3 — angle-test path:**

```ts
it('uses projected endpoint when angle test passes', () => {
  // Сконструировать задачу, где du, спроецированное на бокс, даёт
  // направление d_k = P(xc + du) - x_k с углом < 90° к -g.
  const n = 2;
  const mat = new BFGSMat(n, 5);
  mat.update(
    new Float64Array([0.3, 0.2]),
    new Float64Array([0.5, 0.3]),
    0,
    1,  // commit 4
  );
  const x = new Float64Array([0.5, 0.5]);
  const g = new Float64Array([1, 1]);
  // … подобрать xc, c, freeSet так, чтобы projected endpoint был
  // в strong-descent направлении.
  // Ожидаем: ws.xHat ≠ xc и (ws.xHat - x_k) · g < -η ‖d‖ ‖g‖.
});
```

**N4 — truncation fallback path:**

```ts
it('falls back to 1997 truncation when projected endpoint fails angle test', () => {
  // Сконструировать du такое, что P(xc + du) - x_k почти ортогонально g.
  // Тогда angle test проваливается, должна сработать ветка усечения.
  // Ожидаем: ws.xHat = x + α (xc + du - x) для некоторого α ∈ (0, 1)
  // ИЛИ ws.xHat = xc (если α = 0).
});
```

Конкретные numerical setups для N3/N4 потребуют подбора `mat`, `xc`, `c`, `freeSet` — проще всего получить их прогоном `cauchyPoint` на специально подобранной задаче, как делают существующие тесты «model decrease» и «Newton-like».

#### Валидация

После 8-го коммита прогнать **полный** интеграционный suite. Bounded-Rosenbrock test (`lbfgs-b-integration.test.ts`, секция «minimize Rosenbrock 2D on [-1.5, 0.5]²») — главный канарейка: он вовлекает active-bounds + subspace, и старая (неправильная) реализация на нём могла случайно работать благодаря fallback на projected SD из commit 6 (теперь снятому). Если этот тест после коммита 8 ломается — значит, либо `η` подобрана неверно (см. таблицу калибровки в §«Возможные осложнения»), либо в truncation-цикле промах в индексах `nbd`.

---

## Валидация после всех коммитов

### Локальный CI

После каждого коммита (не только финального):

```
npm run lint-fix && npm run build && npm test
```

Все существующие тесты должны быть зелёными. Если коммит N меняет тесты — изменения в тестовых файлах должны быть в **том же** коммите, что изменения в коде (атомарность).

### Бенчмарки

Сравнение с baseline-снапшотом из commit 0 — основной критерий регрессии. После каждого коммита (не только финального) прогонять **релевантные** runner-ы согласно таблице из §«Сводная таблица правок» и сравнивать с `/tmp/baseline-*.md`:

1. `unconstrained-benchmarks.ts` — ожидание: качество L-BFGS-B без bounds ≥ качества L-BFGS. Регрессия здесь маловероятна (subspace-фаза в unconstrained случае почти не меняется), но если есть — сигнал об ошибке в commit 4 (curvature gate) или 3 (test reorder).

2. `multistart-benchmarks.ts` — guards against x₀-sensitivity. Здесь ожидается **улучшение**: правки 1, 2, 7, 8 убирают пути, на которых старая реализация теряла итерации.

3. `bounded-benchmarks.ts` — главная цель. На задачах с активными bounds (bounded Rosenbrock, half-bounded, fixed-variable) ожидается значительное улучшение по числу итераций после commit 8 (правильный M&N 2011 = смысл всего апдейта 2011 года). Также вход для калибровки `η`.

**Правило интерпретации:** «улучшение» = средн. итераций ↓ или success rate ↑ при том же или лучшем feasibility violation. Просто меньше итераций при возросшем `feas_vio` — **не** улучшение.

### SciPy parity (если применимо)

Если в проекте есть инфраструктура `§12.3` (cross-validation против SciPy) — после commit 4 (curvature gate) и commit 3 (test reorder) она должна показать numerical parity (final f до 1e-8, final x до 1e-6) на стандартном наборе задач. До commit 4/3 parity не гарантируется из-за расхождения в `nskip` и порядке update/test.

---

## Возможные осложнения

### Калибровка `η` в commit 8

Значение `η = 1e-2` взято из `subsm.f`. Решающее правило для калибровки (привязка к `baseline-bounded.md` из commit 0):

| Метрика на bounded-benchmarks | Действие |
|---|---|
| Все 7 задач сходятся, средн. итераций ≤ baseline ÷ 1.5 | **OK**, оставить `η = 1e-2`. |
| Хотя бы 1 задача не сходится (`feas_vio < 1e-6` нарушено или `f − f*` > tol) | Попробовать `η = 1e-6` — слабее порог, ближе к «direction is not actively ascent». |
| Все сходятся, но средн. итераций ≥ baseline (т. е. правка не дала ускорения, хоть и не регрессия) | Попробовать `η = epsmch^(1/3) ≈ 6e-6`; если не помогло — оставить `1e-2`, искать причину в реализации (проверить порядок проекции, знак `du`). |
| Регрессия на одной задаче, ускорение на остальных | Залогировать, оставить `η = 1e-2`, отдельным follow-up-PR разобрать конкретную задачу. |

Любую калибровку вне `1e-2` отметить в TSDoc к константе с комментарием «calibrated against bounded-benchmarks at <commit hash>».

### Удаление `modelDelta` ломает что-нибудь?

`modelDelta` сейчас используется только внутри Morales–Nocedal backtrack-цикла и нигде больше. После commit 8 функция полностью не нужна. Если в кодовой базе появятся другие потребители — оставить, но в текущей реализации можно удалить.

### Связь с finiteDifferences fallback

Все правки в драйвере не затрагивают путь `evalGrad` (который при отсутствии `gradFn` использует FD). Тестов с FD в существующем suite нет (все интеграционные тесты передают `gradFn` либо неявно через сходимость на простых задачах). Это пробел тестового покрытия, но он не относится к текущему набору правок и решается отдельно.
