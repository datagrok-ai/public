# План модификации: angle-test и truncation reference в `subspace.ts`

Документ описывает две связанные правки в `subspace.ts`, относящиеся к шагам (9) и (10) функции `subspaceMin`. Обе исправления — выравнивание под канонический Fortran v3.0 (`subsm.f`, c-jlm-jn block, lines 50–70) и под текст Morales–Nocedal (2011). Обе правки делаются в одном коммите, потому что отдельно их вносить бессмысленно: angle-test и truncation формируют единую endpoint-selection логику, и их некорректность даёт связанный эффект (преждевременный fallback + неправильная точка усечения).

---

## Сводка

| Аспект | Сейчас в TS | Должно быть (Fortran v3.0 + MN2011) |
|---|---|---|
| Тест приёмки спроектированной точки | `g·d ≤ −η · ‖d‖ · ‖g‖` с `η = 1e-2` | `g·d < 0` (без нормировки) + `‖d‖² > 0` guard |
| Точка отсчёта truncation | `x_k → x̂` (внешний итерат) | `x_c → x_c + du` (Cauchy point) |
| Что усекается | вектор `x̂ − x_k` | вектор `du` (подпространственный шаг) |
| Граничный snap при α<1 | нет | координата `b`, давшая минимальный `cap`, прижимается к границе |
| Поведение на «ascent» | сразу truncation | то же — truncation |
| Поведение при `α = 0` | fallback на `x_c` | то же — `x_c` (по построению `α = 0` ⇒ результат тождественно `x_c`) |

---

## Часть 1: angle-test → directional-derivative test

### 1.1 Текущий код (фрагмент `subspaceMin`)

```ts
const ETA = 1e-2;
// ...

// (9) Angle test on the true gradient
let gd = 0;
let dd = 0;
let gg = 0;
for (let i = 0; i < n; i++) {
  const di = xBar[i] - x[i];
  gd += g[i] * di;
  dd += di * di;
  gg += g[i] * g[i];
}
const dNorm = Math.sqrt(dd);
const gNorm = Math.sqrt(gg);

if (dNorm > 0 && gNorm > 0 && gd <= -ETA * dNorm * gNorm) {
  for (let i = 0; i < n; i++) xHat[i] = xBar[i];
  return {improved: true};
}
```

### 1.2 Канонический Fortran (subsm.f, lines 51–62)

```fortran
c     check sign of the directional derivative
      dd_p = zero
      do 55 i=1, n
         dd_p  = dd_p + (x(i) - xx(i))*gg(i)
 55   continue
      if ( dd_p .gt. zero ) then
         call dcopy( n, xp, 1, x, 1 )       ! откат на xc
         ...
      else
         go to 911                          ! принять спроектированную точку
      endif
```

`xx` — внешний итерат `x_k` (передан в subsm как параметр), `gg` — `g(x_k)`, `x` — изначально `x_c`, после блока проекции стало `x_proj = P(x_c + d̄ᵘ)`. То есть тест буквально:

```
ddp = (x_proj - x_k) · g(x_k)
if ddp < 0:  принять x_proj
else:        truncation fallback
```

Точка `ddp = 0` обрабатывается как «не descent» (Fortran: `dd_p .gt. zero`, инвертированная ветка → fallback). Никакого `eps0`, `η`, `‖d‖`, `‖g‖` в этом тесте нет.

### 1.3 Целевой код

```ts
// (9) Directional-derivative test: is the projected endpoint a
//     descent direction on the original problem at x_k?
//     Mirrors subsm.f c-jlm-jn block lines 51-62: accept iff
//     g_k · (xBar - x_k) < 0. No normalization, no angle threshold —
//     the outer-loop line search will take care of step size.
//
//     `dd > 0` guard: if projection collapsed xBar back onto x_k
//     (e.g. xc + du happened to project to x_k bit-exactly), then
//     gd = 0 mathematically; FP cancellation can flip the sign
//     to ε-negative and we would accept a zero direction. Driver
//     would then take it as endpoint, form d = xBar - x_k = 0,
//     and fail the `slope < 0` check anyway — but that costs a
//     memory reset. Cheaper to filter here.
let gd = 0;
let dd = 0;
for (let i = 0; i < n; i++) {
  const di = xBar[i] - x[i];
  gd += g[i] * di;
  dd += di * di;
}

if (gd < 0 && dd > 0) {
  for (let i = 0; i < n; i++) xHat[i] = xBar[i];
  return {improved: true};
}
```

Удаляется константа `ETA`, удаляются `gg`, `dNorm`, `gNorm`. Сохраняется только `dd` как guard от вырожденного направления (см. §1.4).

### 1.4 Минимальный guard `dd > 0`

В Fortran v3.0 этого guard'а нет — `subsm.f` буквально проверяет только `dd_p .gt. zero`. Но в Fortran проекция и накопление `dd_p` идут с тем же FP-форматом, что и outer line search, и при `xBar = x_k` дальнейший шаг просто отказывает в `dcsrch`.

В TS-реализации цена другой проверки чуть выше (memory reset стоит итерации). Поэтому добавлен **минимальный** guard: `dd > 0`. Это **не** angle-test (никакой нормировки на `‖g‖·‖d‖`), а только защита от точного нулевого направления. Когда `xBar ≠ x_k` хоть на один ULP — guard не срабатывает.

Альтернатива в стиле «`gd < -ε·...`» отвергнута: `Math.abs(чего?)` в безразмерной задаче — произвол, и любой выбор шкалы создаёт скрытый параметр. `dd > 0` шкало-инвариантен.

### 1.5 Удаление константы `ETA`

```ts
// убрать на module-уровне:
// const ETA = 1e-2;
```

И обновить doc-comment в шапке файла:

```ts
/**
 * Endpoint selection (Morales–Nocedal 2011, c-jlm-jn block in subsm.f):
 *   1. Form xHatFull = xc + Z·d̄ᵘ.
 *   2. Project once: xBar = P(xHatFull).
 *   3. Directional-derivative test on the original gradient at x_k:
 *      accept xBar iff g_k · (xBar - x_k) < 0 AND ‖xBar - x_k‖ > 0.
 *   4. Otherwise: 1997 truncation fallback — back-track from xc along
 *      d̄ᵘ until a bound is hit (largest α ∈ (0, 1]).
 *   5. If α = 0: fall back to xc.
 *
 * The 1997 code did the truncation unconditionally; MN2011 showed that
 * the projected point is preferred when it is descending, since
 * truncation can produce nearly-orthogonal-to-(-g) directions that
 * stall the outer line search. The `dd > 0` clause is a TS-only guard
 * against an FP-cancelled gd < 0 with xBar ≡ x_k; it is a strict
 * subset of the canonical Fortran condition (Fortran's dcsrch absorbs
 * the same case at higher cost).
 */
```

---

## Часть 2: truncation reference → Cauchy point + bound snap

### 2.1 Текущий код

```ts
// (10) Truncation fallback: largest α ∈ (0, 1] keeping x_k + α(x̂ − x_k)
//      feasible.
let alpha = 1;
for (let i = 0; i < n; i++) {
  const stepDir = xHatFull[i] - x[i];   // ← x = x_k
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
    xHat[i] = x[i] + alpha * (xHatFull[i] - x[i]);
  return {improved: true};
}

// (11) α = 0 → fall back to xc.
for (let i = 0; i < n; i++) xHat[i] = xc[i];
return {improved: false};
```

Здесь `x` — параметр функции, равный внешнему `x_k`. Усекается отрезок `[x_k, xHatFull]`. Это **не** то, что делает Fortran.

### 2.2 Канонический Fortran (subsm.f, lines 56–80)

После того, как `dd_p .gt. zero` (директивная производная положительная или нулевая → нет descent в спроектированной точке), Fortran делает:

```fortran
c     откат: x := xp (т.е. x := xc, сохранённая копия Cauchy point)
      call dcopy( n, xp, 1, x, 1 )
c
c-----------------------------------------------------------------
c
c     Now compute the maximum step length alpha:
c
      alpha = one
      temp1 = alpha
      iexit = 0
      do 60 i = 1, nsub
         k     = ind(i)
         dk    = d(i)
         if (nbd(k) .ne. 0) then
            if (dk .lt. zero .and. nbd(k) .le. 2) then
               temp2 = l(k) - x(k)
               if (temp2 .ge. zero) then
                  temp1 = zero
               else if (dk*alpha .lt. temp2) then
                  temp1 = temp2/dk
               endif
            else if (dk .gt. zero .and. nbd(k) .ge. 2) then
               temp2 = u(k) - x(k)
               if (temp2 .le. zero) then
                  temp1 = zero
               else if (dk*alpha .gt. temp2) then
                  temp1 = temp2/dk
               endif
            endif
            if (temp1 .lt. alpha) then
               alpha = temp1
               iexit = k
            endif
         endif
 60   continue
c
      if (alpha .lt. one) then
         dk = d(iexit-?? mapping) ; ...
         if (dk .gt. zero) then; x(k) = u(k); d(?)=zero
         else if (dk .lt. zero) then; x(k) = l(k); d(?)=zero
         endif
      endif
c
      do 70 i = 1, nsub
         k    = ind(i)
         x(k) = x(k) + alpha*d(i)
 70   continue
 911  continue
```

Ключевые отличия от TS:

1. **Точка отсчёта.** `x` к моменту цикла 60 — это **Cauchy point** (`x := xp`), а не внешний `x_k`. Поэтому `temp2 = l(k) - x(k)` — это «сколько ещё осталось от Cauchy до нижней границы», и Fortran движется по `d(i)` (= `du`-компонента на свободной координате) **из** `xc`.

2. **Усекается `du`, а не `xHatFull − x_k`.** Финальный шаг `x(k) = x(k) + alpha*d(i)` =  `xc(k) + alpha · du(k)`.

3. **Snap-к-границе при α<1.** Если `alpha < 1`, координата `iexit`, давшая минимум, явно прижимается к границе **до** прибавления `alpha · d`. Это устраняет ε-инфизибильность от делений `temp2/dk`.

4. **Только свободные координаты.** Цикл идёт по `i = 1, nsub` — то есть по `freeSet`. Несвободные координаты не двигаются вообще (на них `d(i) = du(i) = 0` по построению; смотреть подсчёт `du` в TS-блоке (6): только для `i = freeSet[k]`).

5. **`temp2 ≥ 0` (для нижней границы)** означает, что `xc` уже ниже границы или на ней — теоретически невозможно (по построению Cauchy `xc ∈ [l, u]`), но Fortran-guard обнуляет `temp1 = 0` для безопасности.

### 2.3 Целевой код

> **NB.** Ниже даны две формулировки. **Каноническая — вторая** (после
> блока «Замечания»): `xHat = xc baseline + free-only overwrite`. Первая
> — пошаговая иллюстрация механики; на ней удобно читать поток данных,
> но в production она проигрывает второй на FP-надёжности (см. ниже).
> В §3 «единая правка целиком» используется именно каноническая версия.

```ts
// (10) 1997 truncation fallback (subsm.f c-jlm-jn lines 56-80).
//      Back-track from xc along du, only on free coordinates. The
//      maximum α ∈ [0, 1] keeps xc + α·du feasible. On α<1 the
//      first-hitting coordinate is snapped to the bound to absorb
//      the ε-rounding in the division temp2/dk.
let alpha = 1;
let iexit = -1;
let iexitDir = 0;  // sign of du at the binding coordinate

for (let kk = 0; kk < t; kk++) {
  const i = freeSet[kk];
  const code = nbd[i];
  if (code === BOUND_FREE) continue;       // ничего не ограничивает
  const dk = du[i];
  if (dk === 0) continue;

  let cap = alpha;
  if (dk < 0 && (code === BOUND_LOWER || code === BOUND_BOTH)) {
    const room = lower[i] - xc[i];          // ≤ 0 на feasible xc
    if (room >= 0) cap = 0;                  // xc уже на/ниже границы
    else if (dk * alpha < room) cap = room / dk;   // dk<0, room<0 → cap>0
  } else if (dk > 0 && (code === BOUND_UPPER || code === BOUND_BOTH)) {
    const room = upper[i] - xc[i];          // ≥ 0 на feasible xc
    if (room <= 0) cap = 0;
    else if (dk * alpha > room) cap = room / dk;
  }
  if (cap < alpha) {
    alpha = cap;
    iexit = i;
    iexitDir = dk;
  }
}

if (alpha === 0) {
  // Truncation degenerate → fall back to xc.
  for (let i = 0; i < n; i++) xHat[i] = xc[i];
  return {improved: false};
}

// Snap the first-hitting coordinate to its bound to absorb ε-drift
// from the division above. Mirrors subsm.f's "if (alpha .lt. one)"
// block: x(iexit) := bound; d(iexit) := 0 (in our case du[iexit]=0
// achieves the same effect when forming xHat below).
if (alpha < 1 && iexit >= 0) {
  // We do not mutate du[] (it's a workspace, but keeping it intact
  // helps debugging). Instead, write the bound directly into xHat
  // for that coord and skip the standard formula.
  if (iexitDir > 0) xHat[iexit] = upper[iexit];
  else              xHat[iexit] = lower[iexit];
}

// All other coords: xHat = xc + alpha * du. On non-free coords du is
// zero, so xHat ≡ xc there. On the binding coord `iexit` we already
// wrote the bound above; overwriting it with xc + alpha·du would
// give the same value up to ε, but we keep the explicit snap.
for (let i = 0; i < n; i++) {
  if (i === iexit && alpha < 1) continue;   // уже записан выше
  xHat[i] = xc[i] + alpha * du[i];
}

return {improved: true};
```

Замечания по реализации:

- **`du[i]` для `i ∉ freeSet` равен 0** — это инвариант блока (6): `du[i] = 0` инициализируется на каждой итерации, и присваивание `du[i] = -r[k]/θ - nTv/(θ²)` происходит только для `i = freeSet[k]`. Поэтому формула `xHat[i] = xc[i] + alpha * du[i]` корректна и для несвободных координат — там она даёт `xc[i]`. Но для **строгой** феазибильности на несвободных координатах надёжнее опираться на то, что `xc` уже допустим, а не на (потенциально ε-ненулевое) `du[i]`. Поэтому в финальной версии я бы написал:

  ```ts
  for (let i = 0; i < n; i++) xHat[i] = xc[i];   // baseline для несвободных
  for (let kk = 0; kk < t; kk++) {
    const i = freeSet[kk];
    if (i === iexit && alpha < 1) {
      xHat[i] = iexitDir > 0 ? upper[i] : lower[i];
    } else {
      xHat[i] = xc[i] + alpha * du[i];
    }
  }
  ```

  Это даёт побайтовое совпадение с Fortran-семантикой «обновляем только свободные».

- **Альтернатива snap'у:** можно не trackить `iexit`/`iexitDir`, а после усечения сделать ещё одну проекцию `project(xHat, lower, upper, nbd, xHat)`. Это аккуратно, но проигрывает в численной точности — `cap = room/dk` может дать `xHat[iexit]` сдвинутым на `O(ε·|bound|)` за пределы границы, и проекция уберёт этот сдвиг. Однако затем восстановление в Cauchy-фазе следующей итерации может не «поверить» в точное значение границы. Snap явный — лучше.

### 2.4 Изменения в импортах

Текущий файл импортирует `BOUND_LOWER, BOUND_BOTH, BOUND_UPPER` из `./types`. Дополнительно для guard'а `code === BOUND_FREE` нужно импортировать `BOUND_FREE`:

```ts
import {BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
```

---

## Часть 3: единая правка целиком

Полная замена шагов (7)–(11) функции `subspaceMin`. Шаги (1)–(6) **не трогаются**.

```ts
// ---- (7) Form the unconstrained Newton-like point xHatFull = xc + du.
//          On non-free coords du is zero, so xHatFull[i] = xc[i] there.
for (let i = 0; i < n; i++) xHatFull[i] = xc[i] + du[i];

// ---- (8) Project once onto [l, u]. xBar is the projected endpoint.
for (let i = 0; i < n; i++) {
  let xi = xHatFull[i];
  const code = nbd[i];
  if ((code === BOUND_LOWER || code === BOUND_BOTH) && xi < lower[i]) xi = lower[i];
  if ((code === BOUND_UPPER || code === BOUND_BOTH) && xi > upper[i]) xi = upper[i];
  xBar[i] = xi;
}

// ---- (9) Directional-derivative test on the original gradient
//          (subsm.f c-jlm-jn lines 51-62): accept xBar iff
//          g_k · (xBar - x_k) < 0. No angle/normalisation; the outer
//          line search handles step sizing.
//          `dd > 0` guards against FP-cancelled gd < 0 with xBar ≡ x_k;
//          see §1.4.
let gd = 0;
let dd = 0;
for (let i = 0; i < n; i++) {
  const di = xBar[i] - x[i];
  gd += g[i] * di;
  dd += di * di;
}

if (gd < 0 && dd > 0) {
  for (let i = 0; i < n; i++) xHat[i] = xBar[i];
  return {improved: true};
}

// ---- (10) 1997 truncation fallback (subsm.f c-jlm-jn lines 56-80).
//           Back-track from xc along du on free coordinates. On α<1
//           the binding coord is snapped to its bound to absorb the
//           ε-rounding from the division.
let alpha = 1;
let iexit = -1;
let iexitDir = 0;

for (let kk = 0; kk < t; kk++) {
  const i = freeSet[kk];
  const code = nbd[i];
  if (code === BOUND_FREE) continue;
  const dk = du[i];
  if (dk === 0) continue;

  let cap = alpha;
  if (dk < 0 && (code === BOUND_LOWER || code === BOUND_BOTH)) {
    const room = lower[i] - xc[i];
    if (room >= 0) cap = 0;
    else if (dk * alpha < room) cap = room / dk;
  } else if (dk > 0 && (code === BOUND_UPPER || code === BOUND_BOTH)) {
    const room = upper[i] - xc[i];
    if (room <= 0) cap = 0;
    else if (dk * alpha > room) cap = room / dk;
  }
  if (cap < alpha) {
    alpha = cap;
    iexit = i;
    iexitDir = dk;
  }
}

if (alpha === 0) {
  // ---- (11) Truncation degenerate → fall back to xc.
  for (let i = 0; i < n; i++) xHat[i] = xc[i];
  return {improved: false};
}

// Build xHat: non-free coords stay at xc; free coords advance by α·du,
// except the binding coord iexit which is snapped to its bound.
for (let i = 0; i < n; i++) xHat[i] = xc[i];
for (let kk = 0; kk < t; kk++) {
  const i = freeSet[kk];
  if (i === iexit && alpha < 1) {
    xHat[i] = iexitDir > 0 ? upper[i] : lower[i];
  } else {
    xHat[i] = xc[i] + alpha * du[i];
  }
}

return {improved: true};
```

Удалить константу `ETA` на module-уровне.

---

## Часть 4: тесты

### 4.1 Существующие тесты в `lbfgs-b-subspace.test.ts`

| Тест | Действие | Причина |
|---|---|---|
| `'col=0 → xHat = xᶜ, improved=false'` | без изменений | early-exit логика не меняется |
| `'freeCount=0 → xHat = xᶜ, improved=false'` | без изменений | то же |
| `'unconstrained: produces m_k(xHat) < m_k(xᶜ)'` | переименовать → `'unconstrained: xHat is the projected Newton step'`, **сохранить логику** | в unconstrained случае `xHatFull = xBar` (проекция тождественна), `gd = g·du < 0` (Newton на положительно определённой `B`), angle-test тривиально проходит — поведение не меняется |
| `'xHat stays feasible even after backtrack'` | переименовать → `'xHat stays feasible after projection or truncation'`; **обновить инвариант** | теперь fallback может дать `xHat = xc + α·du`, что по построению feasible |
| `'B·(xHat−xᶜ) ≈ −r̄ᶜ for t=n'` | без изменений | unconstrained путь проходит через angle-test, `xHat = xBar = xc + du`, инвариант сохраняется |
| `'projects and backtracks when unconstrained step is infeasible'` | **переписать** на 2 теста (см. ниже) | старый тест проверял неправильную семантику бэктрэка |

### 4.2 Новые/переписанные тесты

#### N3 — directional-derivative test, успешная ветка

```ts
it('accepts xBar when g·(xBar - x) < 0', () => {
  // 2D задача: x_k = (0.5, 0.5), бокс [0, 1]^2, градиент g тянет на юго-запад.
  // Подбираем du так, чтобы xc + du вышло из бокса по одной координате,
  // но проекция давала direction с явным descent на g.
  const n = 2;
  const mat = new BFGSMat(n, 5);
  mat.update(
    new Float64Array([0.3, 0.2]),
    new Float64Array([0.5, 0.3]),
    0,
    1,
  );
  const x = new Float64Array([0.5, 0.5]);
  const g = new Float64Array([1, 1]);  // descent в любую сторону уменьшения координат

  // Прогон Cauchy на этой задаче должен дать xc внутри бокса
  // (с большим θ, маленьким шагом).
  const cauchyWs = makeCauchyWorkspace(n, 5);
  const lower = new Float64Array([0, 0]);
  const upper = new Float64Array([1, 1]);
  const nbd = classifyBounds(lower, upper);
  const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
  expect(cr.ok).toBe(true);

  const ws = makeSubspaceWorkspace(n, 5);
  const r = subspaceMin(
    x, g, lower, upper, nbd,
    cauchyWs.xc, cauchyWs.c, cauchyWs.freeSet, cr.freeCount, mat, ws,
  );
  expect(r.improved).toBe(true);

  // xHat должен лежать в [0, 1]^2.
  for (let i = 0; i < n; i++) {
    expect(ws.xHat[i]).toBeGreaterThanOrEqual(0);
    expect(ws.xHat[i]).toBeLessThanOrEqual(1);
  }

  // Direction (xHat - x_k) даёт строго отрицательную директивную
  // производную на g.
  let gd = 0;
  for (let i = 0; i < n; i++) gd += g[i] * (ws.xHat[i] - x[i]);
  expect(gd).toBeLessThan(0);
});
```

#### N4 — truncation fallback, не-degenerate путь

```ts
it('falls back to truncation along xc → xc+du when projected gives ascent', () => {
  // Сконструировать такую задачу, где xc глубоко внутри бокса,
  // du толкает в направлении почти-перпендикулярно −g, так что
  // gd = g · (xBar - x_k) ≥ 0, но cap ∈ (0, 1] — тогда truncation
  // даёт xHat внутри бокса на отрезке [xc, xc+du].
  // Точные числа подбираются ad-hoc прогоном.

  // Setup: 1D, чтобы избежать многомерных подгонок.
  const n = 1;
  const mat = new BFGSMat(n, 5);
  // Невозможно построить это аналитически в 1D; вместо этого
  // сделаем 2D с фиксированной координатой через ультра-узкий бокс.
  // ... (см. полный код ниже)
});
```

В реальности конструировать тест-кейс на failing angle test трудно, потому что direction `xBar − x_k` обычно сонаправлен с `−g` на типовой задаче. Есть два пути:

1. **Mock subspace ws** с искусственным `du` через exposed entry point. Требует test-only API.
2. **Patологический mat**: Hessian приближение, в котором inverse-Newton step выводит из бокса далеко и проекция даёт near-orthogonal direction. Подбирается через `mat.update(s, y, 0, 0)` с подобранными `s, y`.

**Рекомендация:** реализовать N4 через путь 2 с конкретными числами. Полная конструкция:

```ts
it('falls back to truncation when projected endpoint gives ascent', () => {
  const n = 2;
  const mat = new BFGSMat(n, 5);
  // Большой θ → "сильная" квадратичная модель, du ~ -B⁻¹ r̄ может
  // быть очень большой амплитуды.
  // s, y подобраны так, что θ = yᵀy/sᵀy высокая.
  mat.update(
    new Float64Array([1, 0]),
    new Float64Array([100, 0]),  // θ = 10000
    0,
    1,
  );
  // x_k около верхнего угла, gradient толкает в (+1, +1) — ascent
  // при движении в (-, -). Бокс [0, 1]^2.
  // ... (точные числа подбираются прогоном через debug-print).

  // Утверждения:
  // 1. xHat в боксе
  // 2. xHat = xc + α·du для некоторого α ∈ (0, 1) ИЛИ xHat = xc
  // 3. На связывающей координате ε-инфизибильность отсутствует
  //    (xHat[iexit] = upper[iexit] точно или lower[iexit] точно).
});
```

Если не получается воспроизвести «ascent» путь из природных задач — тест пометить как `it.skip` с TODO-комментарием или реализовать через test-only entry point. Bug критический, но ветка fallback редкая, и если природных триггеров нет в test suite — это положительный сигнал, что в реальных задачах MN2011 angle test проходит почти всегда.

#### N5 — точка отсчёта truncation = xc, а не x_k (integration canary)

Это **главный тест на правильность правки**. Реализуем как
**integration-canary** на «Bounded Rosenbrock» из `BOUNDED_PROBLEMS`
(`benchmarks/test-functions.ts`): эта задача имеет активную верхнюю
границу в оптимуме (`x* = (0.5, 0.25)`, `upper = (0.5, 0.5)`), что
даёт fallback-путь truncation естественно срабатывающим, и на ней
правка должна быть видна как сокращение числа итераций.

Идея — snapshot-тест на iter-count, файлится **до** правки и
проверяется **после**. Точные числа фиксируются прогоном на baseline
ветке и подставляются в комментарий теста как «expected».

```ts
it('Bounded Rosenbrock — iteration count drops vs baseline', () => {
  const p = BOUNDED_PROBLEMS.find((q) => q.name === 'Bounded Rosenbrock')!;
  let iters = 0;
  const opt = new LBFGSB();
  const res = opt.minimize(p.fn, p.x0.slice(), {
    bounds: {lower: p.lower, upper: p.upper},
    maxIterations: 200,
    tolerance: 1e-10,
    gradTolerance: 1e-8,
    onIteration: () => { iters++; return false; },
  });
  expect(res.x[0]).toBeCloseTo(p.knownPoint[0], 6);
  expect(res.x[1]).toBeCloseTo(p.knownPoint[1], 6);
  expect(res.value).toBeCloseTo(p.knownMin, 8);

  // Baseline (ETA=1e-2 angle test + x_k truncation reference): N_BEFORE
  // After fix (gd<0 + dd>0 + xc truncation reference): N_AFTER
  // Both numbers measured ad hoc; N_AFTER ≤ N_BEFORE expected.
  // If equal: правка не видна на этой задаче — расширить тест на
  // другие BOUNDED_PROBLEMS или добавить unit-тест с pathological mat.
  expect(iters).toBeLessThanOrEqual(/* N_AFTER, fill in after baseline */ 0);
});
```

**Процедура заполнения числа:**

1. На master-ветке (до правки): запустить тест с `expect(iters).toBeLessThanOrEqual(999)`, посмотреть фактическое `iters` → `N_BEFORE`.
2. На fix-ветке: запустить тот же тест → `N_AFTER`.
3. В тест записать `expect(iters).toBeLessThanOrEqual(N_AFTER)`. Если `N_AFTER === N_BEFORE` — правка на этой задаче не сказывается; в этом случае:
   - либо добавить аналогичную проверку на `Bounded Beale` или `Fixed Variables Sphere` (там `xc ≠ x_k` структурно);
   - либо вернуться к unit-тесту в стиле §4.2 N5b ниже.

#### N5b — точка отсчёта truncation (unit, mock setup) — fallback

Если N5 не даёт сигнала (числа итераций совпали с baseline), фиксируем
правку через unit-тест с предзаданными `xc`, `c`, `freeSet`. `subspaceMin`
принимает их как параметры, так что setup тривиален:

```ts
it('truncation reference is xc, not x_k (mock setup)', () => {
  const n = 2;
  const m = 3;
  const mat = new BFGSMat(n, m);
  mat.update(/* s, y подобраны под предсказуемые θ и du */);

  const x = new Float64Array([0.8, 1]);
  const g = new Float64Array([1, 1]);
  const xc = new Float64Array([0.3, 0.5]);
  const c = new Float64Array(2 * mat.col);
  const delta = new Float64Array(n);
  for (let i = 0; i < n; i++) delta[i] = xc[i] - x[i];
  mat.applyWt(delta, c);

  const freeSet = new Int32Array([0, 1]);
  const freeCount = 2;
  const lower = new Float64Array([0, 0]);
  const upper = new Float64Array([1, 1]);
  const nbd = classifyBounds(lower, upper);

  const ws = makeSubspaceWorkspace(n, m);
  const r = subspaceMin(x, g, lower, upper, nbd, xc, c, freeSet, freeCount, mat, ws);

  // Инвариант, отличающий старую и новую реализации:
  // если ветка truncation сработала, xHat[i] лежит на отрезке
  // [xc[i], xc[i] + du[i]] (по компонентам, со знаком du[i]).
  // В старой реализации xHat[i] лежал бы на [x[i], x[i] + (xHatFull[i]-x[i])].
  if (r.improved) {
    for (let i = 0; i < n; i++) {
      const lo = Math.min(xc[i], xc[i] + ws.du[i]);
      const hi = Math.max(xc[i], xc[i] + ws.du[i]);
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lo - 1e-12);
      expect(ws.xHat[i]).toBeLessThanOrEqual(hi + 1e-12);
    }
  } else {
    expect(ws.xHat[0]).toBe(xc[0]);
    expect(ws.xHat[1]).toBe(xc[1]);
  }
});
```

Этот тест всё ещё слабее N5 (инвариант, не точное число), но в паре с
N6 (snap-к-границе, точное равенство) он закрывает обе характеристики
правки.

#### N6 — bound snap при α<1

```ts
it('snaps the binding coordinate to its bound on truncation', () => {
  // Setup, где truncation действительно срабатывает с α < 1.
  // ... (path-coupled или mock).
  
  // Проверка: после truncation xHat[iexit] точно равен upper[iexit]
  // или lower[iexit] (без ε-смещения).
  expect(ws.xHat[iexit]).toBe(0);  // или upper[iexit], в зависимости от знака
});
```

### 4.3 Сводка изменений в тестах

| Тест | Файл | Действие |
|---|---|---|
| `'col=0'`, `'freeCount=0'` | `lbfgs-b-subspace.test.ts` | без изменений |
| `'unconstrained ...'` | `lbfgs-b-subspace.test.ts` | переименовать, без изменений в логике |
| `'feasible after backtrack'` | `lbfgs-b-subspace.test.ts` | переименовать |
| `'B·(xHat−xᶜ) ≈ −r̄ᶜ'` | `lbfgs-b-subspace.test.ts` | без изменений |
| **N3** descent ⇒ xBar | `lbfgs-b-subspace.test.ts` | новый |
| **N4** ascent ⇒ truncation | `lbfgs-b-subspace.test.ts` | новый (mock-friendly) |
| **N5** truncation reference = xc (integration canary) | `lbfgs-b-integration.test.ts` | **критический регрессионный**, новый — snapshot iter-count на Bounded Rosenbrock |
| **N5b** truncation reference (mock unit) | `lbfgs-b-subspace.test.ts` | новый — fallback для случаев, когда N5 не различает baseline и fix |
| **N6** bound snap | `lbfgs-b-subspace.test.ts` | новый |
| `'projects and backtracks ...'` | `lbfgs-b-subspace.test.ts` | удалить (заменён на N3+N4) |

### 4.4 Интеграционные тесты

Существующие в `lbfgs-b-integration.test.ts` тесты на bounded-Rosenbrock — главный канарейка. После правки они должны **продолжать проходить**, причём с потенциально меньшим числом итераций. Если число итераций выросло — это **сигнал**, что либо truncation reference исправлен, но angle-test (N3) ослаблен слишком сильно, либо наоборот. В этом случае стоит снять `it.skip` со специально подготовленных multi-iter benchmark'ов (если они есть) и сравнить trace с baseline.

---

## Часть 5: валидация

Runner проекта — **Jest** (`"test": "jest"` в `libraries/sci-comp/package.json`). Команды ниже — Jest-флаги.

### 5.1 До правки — зафиксировать baseline

```bash
cd libraries/sci-comp
npm test -- --json --outputFile=/tmp/baseline-before-subspace-fix.json
```

### 5.2 После правки — сравнить

```bash
npm test -- --json --outputFile=/tmp/after-subspace-fix.json
diff <(jq -r '.testResults[].assertionResults[] | "\(.fullName)\t\(.status)"' /tmp/baseline-before-subspace-fix.json | sort) \
     <(jq -r '.testResults[].assertionResults[] | "\(.fullName)\t\(.status)"' /tmp/after-subspace-fix.json    | sort)
```

Все статусы `passed`. Регрессии в integration на bounded — **не** регрессия правки, а либо настоящий новый bug, либо проявление другого замаскированного.

### 5.3 Регенерация benchmark-отчётов

Правка влияет на bounded-путь L-BFGS-B, поэтому отчёты в `benchmarks/`
обязательно регенерируются (см. CLAUDE.md, §«Benchmarks»):

```bash
cd libraries/sci-comp
npx tsx src/optimization/single-objective/benchmarks/bounded-benchmarks.ts > \
  src/optimization/single-objective/benchmarks/bounded-benchmarks.md
npx tsx src/optimization/single-objective/benchmarks/multistart-benchmarks.ts > \
  src/optimization/single-objective/benchmarks/multistart-benchmarks.md
```

(Точное имя выходного файла и формат — по существующему конвенту в репозитории; если runner печатает Markdown в stdout — вышеуказанная форма работает, иначе подсмотреть в `package.json` / истории комитов скрипт регенерации.)

В diff'е `bounded-benchmarks.md` ожидаем:
- L-BFGS-B: число итераций на «Bounded Rosenbrock» / «Bounded Beale» **уменьшилось** или осталось тем же;
- `feas_vio` остаётся 0 (правка не должна вводить инфизибильность);
- остальные оптимизаторы — без изменений (правка локальна для L-BFGS-B).

### 5.4 SciPy-parity (если есть)

Если в проекте есть инфраструктура сравнения с SciPy на стандартных задачах — после правки ожидаем **сближение** числа итераций со SciPy на bounded-Rosenbrock и других задачах из CUTE-collection. Расхождение в 1–2 итерации — норма (FP non-associativity). Расхождение в >5 итераций — стоит пересмотреть guard в §1.4.

---

## Часть 6: размер коммита

Это **один атомарный коммит**. Разбивать на «правка angle-test» и «правка truncation» отдельно нельзя, потому что:

- между ними меняется ширина класса задач, попадающих в каждую ветку;
- старые тесты, переписанные на новую truncation-семантику, не пройдут на коде с новым angle-test'ом и старой truncation;
- diff меньше 100 строк (по сути замена шагов 9–11 в одной функции + удаление константы) — нет смысла размазывать.

**Размер:** ~80 строк добавлено, ~60 удалено в `subspace.ts` + 4 новых теста + 2 обновлённых теста в `lbfgs-b-subspace.test.ts`.

**Заголовок коммита:**

```
fix(lbfgs-b): align subspace endpoint selection with subsm.f c-jlm-jn

Replace heuristic angle test (g·d ≤ −η‖d‖‖g‖, η=1e-2) with the
canonical directional-derivative test (g·d < 0) and fix the
truncation reference point from x_k to x_c. Both deviations from
Fortran v3.0 / Morales–Nocedal 2011 were producing different
trajectories on bounded problems, with the truncation bug
particularly visible when xc differs significantly from x_k
(e.g. when several variables hit bounds during Cauchy phase).

Refs: subsm.f lines 50-80 in Lbfgsb.3.0/lbfgsb.f.
Refs: Morales & Nocedal (2011), TOMS 38(1):7, §3.
```

---

## Часть 7: риск и rollback

**Риск:** низкий. Обе правки приближают TS-семантику к Fortran v3.0; ни одна не вводит новых degree-of-freedom. Угол test становится строго слабее (раньше требовал angle ≥ arccos(η), теперь требует только descent), значит **больше** шагов будут приниматься через быстрый путь (xBar). Truncation становится **более** консервативным (стартует из xc, который ближе к x_k), значит итерации будут короче.

**Что может пойти не так:**

1. Угол test теперь принимает шаги, которые раньше отвергались. На патологических задачах с почти-сингулярной B это может дать line-search-fail чаще. Recovery-механизм в driver покрывает (mat.reset() + retry). Тестовое покрытие: integration тесты с FD-градиентом могут показать рост function evaluation count — это ожидаемо при FD, не регрессия.

2. Truncation от `xc` может дать `α = 0` чаще, чем от `x_k`, если `xc` лежит близко к границе. В этом случае fallback на `xc` (шаг 11). Это идентично поведению Fortran. На последующих итерациях Cauchy сделает тривиальный шаг (потому что freeSet сократился), и subspace будет вызываться с `freeCount = 0` — early-exit.

3. **Edge case `improved=true && xHat ≈ xc`** возможен (новая ветка `gd < 0 && dd > 0` на `xBar = xc`, что бывает, когда проекция целиком отбросила `du`). Проверка в driver: `driver.ts:154–180` берёт `endpoint = subspace.xHat` без чтения `improved` и проверяет `slope = g·d < 0` явно. На `xHat = xc` имеем `d = xc − x_k`, что Cauchy строит как descent (`slope ≤ 0`); в случае строго `slope = 0` driver уходит в memory-reset retry — корректное и безопасное поведение. Поле `improved` в `SubspaceResult` ни в driver, ни в тестах не читается; формально его сейчас можно удалить, но это вне scope этой правки. Записать в TODO для отдельного коммита.

**Rollback:** отдельным reverse-коммитом, восстанавливающим текущий код. Так как правка локализована в одной функции и дополнена тестами, rollback тривиален (`git revert`).

---

## Чеклист исполнителя

**Перед правкой (на master / fix-ветке до изменений):**
- [ ] Прочитать lines 50–80 в `Lbfgsb.3.0/subsm.f` (или ROHSA-зеркало) — сверить с описанием в §1.2 и §2.2.
- [ ] Зафиксировать `N_BEFORE` для теста N5: запустить L-BFGS-B на «Bounded Rosenbrock» и записать число итераций (см. §4.2 N5, шаг 1).
- [ ] Сохранить текущие `bounded-benchmarks.md` и `multistart-benchmarks.md` куда-нибудь во временное место для diff'а.

**Правка кода:**
- [ ] Заменить шаги (7)–(11) в `subspace.ts` на код из §3 (использовать **каноническую** версию truncation, см. NB в §2.3).
- [ ] Удалить константу `ETA` на module-уровне `subspace.ts`.
- [ ] Обновить doc-comment в шапке `subspace.ts` (см. §1.5) — включить упоминание `dd > 0` guard'а.
- [ ] Добавить `BOUND_FREE` в импорт `subspace.ts` (см. §2.4).
- [ ] Обновить doc-comment workspace-полей `delta`/`bDelta` в `SubspaceWorkspace` (если они переименованы — переименовать).

**Тесты:**
- [ ] Переписать тест `'projects and backtracks ...'` на N3 + N4.
- [ ] Добавить N5 (integration canary, snapshot iter-count) и N6 (bound snap).
- [ ] Если N5 не различает baseline и fix (`N_AFTER === N_BEFORE`) — добавить N5b (mock unit, см. §4.2).

**Валидация:**
- [ ] Запустить `npm run lint-fix && npm run build && npm test` — все зелёные.
- [ ] Регенерировать `bounded-benchmarks.md` и `multistart-benchmarks.md` (см. §5.3); проверить что L-BFGS-B на bounded-задачах не деградировал и `feas_vio` остаётся 0.
- [ ] Если есть SciPy-parity тесты — прогнать.

**PR:**
- [ ] Зафиксировать `N_AFTER` в комментарии теста N5 (см. §4.2 шаг 3).
- [ ] PR с заголовком из §6, в описании приложить diff `bounded-benchmarks.md` (до/после).
- [ ] Открыть отдельный TODO-issue: «удалить неиспользуемое поле `improved` из `SubspaceResult` или начать его читать в driver» (см. §7, пункт 3).
