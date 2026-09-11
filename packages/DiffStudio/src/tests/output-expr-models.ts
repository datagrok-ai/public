// Test-only IVP models with #output expressions in #loop / #update models (GROK-20866).

/** Cyclic model whose #output includes an expression using math + #constants. */
export const LOOP_OUTPUT_EXPRESSIONS = `#name: Loop output expressions
#equations:
  dx/dt = -k * x + weight
  dy/dt = k * x

#expressions:
  weight = amp * sin(PI * t)
  energy = k * exp(-t) * (x * x + y * y)

#constants:
  k = 0.5
  amp = 2

#loop:
  count = 3
  x += dose

#argument: t
  start = 0
  finish = 2
  step = 0.1

#inits:
  x = 1
  y = 0

#parameters:
  dose = 0.5

#output:
  t
  x
  energy {caption: Energy}`;

/** Multistage model whose #output includes an expression using math + #constants. */
export const UPDATE_OUTPUT_EXPRESSIONS = `#name: Update output expressions
#equations:
  dx/dt = -k * x + weight
  dy/dt = k * x

#expressions:
  weight = amp * sin(PI * t)
  energy = k * exp(-t) * (x * x + y * y)

#constants:
  k = 0.5
  amp = 2

#argument: t, 1-st stage
  start = 0
  finish = 2
  step = 0.1

#update: 2-nd stage
  duration = extra
  x += boost

#inits:
  x = 1
  y = 0

#parameters:
  extra = 2
  boost = 1

#output:
  t
  x
  energy {caption: Energy}`;
