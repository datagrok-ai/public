import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {PropertyDesirability} from '../mpo';
import {MpoCategoricalEditor} from './mpo-categorical-editor';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';

export interface DesirabilityEditor<T = any> {
  supportsModeDialog: boolean;
  root: HTMLElement;
  onChanged: Subject<T>;
  redrawAll(): void;
  setColumn?(col: DG.Column | null): void;
}


export class DesirabilityEditorFactory {
  static create(
    prop: PropertyDesirability,
    width = 300,
    height = 80,
  ): DesirabilityEditor {
    if (prop.line)
      return new MpoDesirabilityLineEditor(prop, width, height);

    return new MpoCategoricalEditor(prop);
  }
}

