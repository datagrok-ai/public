import {Subject} from 'rxjs';

import {PropertyDesirability} from '../mpo';
import {MpoCategoricalEditor} from './mpo-categorical-editor';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';

export interface DesirabilityEditor<T = any> {
  root: HTMLElement;
  onChanged: Subject<T>;
  redrawAll(): void;
}


export class DesirabilityEditorFactory {
  static create(
    prop: PropertyDesirability,
    width = 300,
    height = 80,
  ): DesirabilityEditor {
    switch (prop.mode) {
    case 'categorical':
      return new MpoCategoricalEditor(prop);
    default:
      return new MpoDesirabilityLineEditor(prop, width, height);
    }
  }
}
