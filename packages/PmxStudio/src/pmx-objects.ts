// (a leaf: ui.ts binds a figure, and the handlers' imports reach back to ui.ts)
import type {StudyInputs as StudyNames} from './explore';

/** An animal of a study. `twin`: Analyze draws it against its untreated twin; `details`: what its click shows. */
export class PmxAnimal {
  constructor(readonly names: StudyNames | null, readonly animal: string,
    readonly o: {color?: string, twin?: () => void, details?: () => HTMLElement} = {}) {}
  toString(): string {
    return this.animal;
  }
}

/** A group of a study, drawn in `color`. `details`: what its click shows. */
export class PmxGroup {
  constructor(readonly names: StudyNames | null, readonly group: string, readonly color: string,
    readonly o: {details?: () => HTMLElement} = {}) {}
  toString(): string {
    return this.group;
  }
}

/** A figure: a chart (its echarts host, named after its panel) or the result's figures, which pin to compare. */
export class PmxFigure {
  constructor(readonly name: string, readonly o: {chart?: HTMLElement, pin?: {on: () => boolean,
    set: (on: boolean) => void}} = {}) {}
  toString(): string {
    return this.name;
  }
}

/** A study of the user's, picked with its box in the Overview's list (a click on its row opens it). `details`: its
 *  panel; `open`: opens it in its stage (`stage`: that stage's name); `compare`: compares the picked studies, null when
 *  fewer than two are picked. */
export class PmxSavedStudy {
  constructor(readonly name: string, readonly o: {key: string, details: () => HTMLElement, stage: string,
    open: () => void,
    compare: () => (() => void) | null}) {}
  toString(): string {
    return this.name;
  }
}

/** A calendar day in the context panel (an element as the current object stays the first one set: the platform does
 *  not replace one element by another). `details`: its panel, actions included. */
export class PmxDay {
  constructor(readonly name: string, readonly details: () => HTMLElement) {}
  toString(): string {
    return this.name;
  }
}
