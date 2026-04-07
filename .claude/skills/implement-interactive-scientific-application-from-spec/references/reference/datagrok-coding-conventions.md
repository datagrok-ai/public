# Coding Standard: Datagrok Interactive Scientific Applications

> **End-to-end example:** application code in `../examples/lotka-volterra-spec/code/src/lotka-volterra/`.

Extracted from the `DiffStudio` package (`packages/DiffStudio`).

## 1. File Structure

### 1.1. Module Organization

Each file has a single responsibility:

- `package.ts` — entry point, package function registration via decorators.
- `app.ts` — main UI application class.
- `constants.ts` — domain constants (not UI).
- `ui-constants.ts` — UI constants: tooltips, titles, errors, links, timeouts.
- `error-utils.ts` — custom errors and error display utilities.
- `utils.ts` — general utilities.
- `model.ts` — types and model class.
- `solver-tools.ts` — wrappers over external library (core).
- `callbacks/` — pattern: base class + concrete implementations in separate files.
- `demo/` — standalone demo models in separate files.
- `tests/` — tests, separated by categories.

### 1.2. Comment at the Beginning of a File

Each file starts with a single-line comment describing the module's purpose:

```typescript
// Solver of initial value problem
```

For files with a detailed description, a block comment is used:

```typescript
/* Scripting tools for the Initial Value Problem (IVP) solver:
     - parser of formulas defining IVP;
     - JS-script generator.
*/
```

### 1.3. Import Order

Imports follow a fixed order:

1. Datagrok API (always three lines):

```typescript
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
```

2. External libraries (`diff-grok`, `codemirror`, `dayjs`, etc.).
3. Internal package modules (`./solver-tools`, `./constants`, etc.).
4. CSS styles (`'../css/app-styles.css'`).

Groups are separated by a blank line.

## 2. Constants

### 2.1. All String Literals Go Into Constants

String literals are not used inline in the code. All messages, titles, tooltips, and keywords are extracted into enums or consts.

### 2.2. Constants Are Grouped Into Enums by Purpose

```typescript
/** Tooltips messages */
export enum HINT {
  HELP = 'Open help in a new tab',
  OPEN = 'Open model',
  // ...
}; // HINT

/** UI titles */
export enum TITLE {
  LOAD = 'Load...',
  IMPORT = 'Import...',
  // ...
}; // TITLE

/** Error messages */
export enum ERROR_MSG {
  SOLVING_FAILS = 'Solving fails',
  // ...
};
```

### 2.3. Enum Naming

- Enum name — `UPPER_CASE` (e.g., `HINT`, `TITLE`, `ERROR_MSG`, `UI_TIME`).
- Enum values — `UPPER_CASE` (e.g., `HINT.OPEN`, `UI_TIME.DOCK_EDITOR_TIMEOUT`).

### 2.4. Composite Constants Are Built from Base Ones

```typescript
const META = `${CONTROL_TAG}meta`;

export enum CONTROL_EXPR {
  NAME = `${CONTROL_TAG}name`,
  SOLVER = `${META}.solver`,
  // ...
};
```

### 2.5. Maps for Linking Constants

`Map` is used for mapping between sets of constants:

```typescript
export const MODEL_HINT = new Map([
  [TITLE.BASIC, HINT.BASIC],
  [TITLE.ADV, HINT.ADV],
  // ...
]);
```

### 2.6. Separation of Domain and UI Constants

- `constants.ts` — domain: parser formulas, solver settings, column names.
- `ui-constants.ts` — UI: tooltips, titles, errors, links, timeouts, dock ratios.

## 3. Types

### 3.1. Local Types Are Defined Near Their Usage

```typescript
/** Numerical input specification */
export type Input = {
  value: number,
  annot: string | null,
};

/** Argument of IVP specification */
type Arg = {
  name: string,
  initial: Input,
  final: Input,
  step: Input,
};
```

### 3.2. Types Used in Multiple Modules Are Exported

Types needed in other files are exported from the defining file and imported at the point of use.

## 4. Classes

### 4.1. Access Modifiers

Fields and methods are explicitly marked with `private` or `public`:

```typescript
export class ModelError extends Error {
  private helpUrl: string;
  private toHighlight: string = undefined;

  public getHelpUrl() { return this.helpUrl; }
  public getToHighlight() { return this.toHighlight; }
}
```

### 4.2. Closing Comment for a Class

A comment with the class name is placed after the closing brace:

```typescript
}; // ModelError
```

```typescript
}; // Model
```

Similarly for large enums and functions:

```typescript
}; // HINT

} // error

} // showModelErrorHint
```

### 4.3. Inheritance Pattern (Callbacks)

The base class is in a separate file, concrete implementations are each in their own file:

```
callbacks/
  callback-base.ts       — base class Callback
  callback-tools.ts      — factory function getCallback
  iter-checker-callback.ts  — IterCheckerCallback extends Callback
  time-checker-callback.ts  — TimeCheckerCallback extends Callback
```

## 5. Functions

### 5.1. JSDoc Comment Before Each Function

Every function (exported and internal) has a single-line JSDoc comment:

```typescript
/** Return solution as a dataframe */
function getSolutionDF(odes: ODEs, solutionArrs: Float64Array[]): DG.DataFrame {

/** Default solver of initial value problem. */
export function solveDefault(odes: ODEs): DG.DataFrame {

/** Return unused IVP-file name */
export function unusedFileName(name: string, files: string[]): string {
```

### 5.2. Arrow Functions for Short Utilities

```typescript
const getMethod = (options?: Partial<SolverOptions>) => {
  // ...
};

const strToVal = (s: string) => {
  const num = Number(s);
  return !isNaN(num) ? num : s === 'true' ? true : s === 'false' ? false : s;
};
```

### 5.3. Comments Inside Functions Mark Logical Steps

```typescript
// Get numerical solution
const approxSolution = method(corProb.odes);

// Compute error
for (let i = 0; i < pointsCount; ++i) {
```

```typescript
// extract function values
const vx = _y[2];

// evaluate expressions
const v = Math.PI * dB ** 3 / 6;

// compute output
_output[0] = vx;
```

## 6. Error Handling

### 6.1. Custom Error Class

A separate class extending `Error` is created for domain errors:

```typescript
export class ModelError extends Error {
  private helpUrl: string;
  constructor(message: string, helpUrl: string) {
    super(message);
    this.helpUrl = helpUrl;
  }
}
```

### 6.2. Factory Functions for Common Errors

```typescript
/** Return ModelError corresponding to ".. is not defined" */
export function getIsNotDefined(msg: string): ModelError {
  // ...
}
```

### 6.3. User Notification via Datagrok Utilities

- `grok.shell.warning(...)` — warnings.
- `grok.shell.error(...)` — errors.
- `grok.shell.info(...)` — informational messages.

## 7. Package Function Registration

Datagrok supports two function registration approaches:

1. **JSDoc-style comments** (traditional) — `//name:`, `//tags: app`, `//input:`, `//output:`. Processed by `grok api` and `grok check`.
2. **Decorators** `@grok.decorators.*` (modern) — type-safe alternative used in newer packages.

Both approaches are valid. The example below uses decorators (as in the DiffStudio package):

```typescript
export class PackageFunctions {
  @grok.decorators.app({
    name: 'Diff Studio',
    description: 'Solver of ordinary differential equations systems',
    browsePath: 'Compute',
  })
  static async runDiffStudio(): Promise<DG.ViewBase> {
    // ...
  }

  @grok.decorators.func({})
  static solve(@grok.decorators.param({type: 'object'}) problem: ODEs): DG.DataFrame {
    return solveDefault(problem);
  }

  @grok.decorators.model({
    name: 'Ball flight',
    description: 'Ball flight simulation',
    // ...
  })
  static ballFlight(/* params */) {
    // ...
  }
}
```

## 8. Formatting (ESLint)

The configuration extends the `google` style guide.

- **Indentation**: 2 spaces.
- **Maximum line length**: 120 characters.
- **Curly braces**: `multi-or-nest` — single-line block without braces, multi-line — with braces.
- **Brace style**: `1tbs` with `allowSingleLine: true`.
- **Unused variables**: `warn`, exceptions — `_`, `ui`, `grok`, `DG`.
- **JSDoc**: `require-jsdoc: off`, `valid-jsdoc: off` — single-line `/** */` comments are used instead.
- **Line break**: `linebreak-style: off`.

## 9. Testing

### 9.1. Test Organization

Tests are separated by categories in individual files:

- `numerical-methods-tests.ts` — solver correctness and performance.
- `features-tests.ts` — IVP format features.
- `platform-funcs-tests.ts` — platform integration.
- `pipeline-tests.ts` — end-to-end pipelines.
- `test-utils.ts` — test utilities.

### 9.2. Framework

`category`, `test`, `expect` from `@datagrok-libraries/test` are used:

```typescript
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category(`Correctness: ${name}`, () => {
  corrProbs.forEach((problem) => test(problem.odes.name, async () => {
    const error = getError(method, problem);
    expect(
      error < TINY,
      true,
      `The ${name} method failed to solve "${problem.odes.name}", too big error: ${error}`,
    );
  }, {timeout: TIMEOUT}));
});
```

### 9.3. Test Options

- `{timeout: N}` — for tests with a time limit.
- `{benchmark: true}` — for performance tests.

## 10. Naming

### 10.1. Variables and Functions

- `camelCase`: `solveDefault`, `getScriptLines`, `showModelErrorHint`.

### 10.2. Classes and Types

- `PascalCase`: `DiffStudio`, `ModelError`, `CallbackAction`, `ModelInfo`.

### 10.3. Enums

- Enum name: `UPPER_CASE` (`HINT`, `TITLE`, `ERROR_MSG`).
- Values: `UPPER_CASE` (`HINT.OPEN`, `UI_TIME.DOCK_EDITOR_TIMEOUT`).

### 10.4. Constants Outside Enums

- `UPPER_CASE` for primitives: `const TINY = 0.0001;`
- `camelCase` for complex objects: `const completions = [...]`, `const modelImageLink = new Map(...)`.

### 10.5. Files

- `kebab-case`: `solver-tools.ts`, `error-utils.ts`, `ui-constants.ts`, `ball-flight.ts`.

## 11. CSS

Styles are placed in a separate CSS file (`css/app-styles.css`) and imported in the modules that use them:

```typescript
import '../css/app-styles.css';
```

CSS classes are named with a package prefix: `diff-studio-hint-btns-div`, `diff-studio-highlight-text`.
