import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';

import {PropertyDesirability} from '../mpo-types';
import {applyDesirabilityPatch} from '../mpo';
import {DesirabilityEditor, DesirabilityEditorFactory} from './desirability-editor-factory';
import {DesirabilityModeDialog} from '../dialogs/desirability-mode-dialog';

import '../../../css/styles.css';

const DEFAULT_WIDTH = 300;
const DEFAULT_HEIGHT = 80;

export interface MpoRowOptions {
  name?: () => string;
  column?: () => DG.Column | null;
  width?: number;
  height?: number;
  design?: boolean;
  onChanged?: () => void;
  onReplaced?: (prop: PropertyDesirability) => void;
}

export interface MpoRow {
  root: HTMLElement;
  propertyCell: HTMLElement;
  weightCell: HTMLElement;
  editorHost: HTMLElement;
  editor: DesirabilityEditor;
  weightInput: DG.InputBase<number | null>;
  applyPatch(patch: Partial<PropertyDesirability>): void;
  sub: Subscription;
}

export function createMpoRow(prop: PropertyDesirability, options: MpoRowOptions = {}): MpoRow {
  const editor = DesirabilityEditorFactory.create(prop, options.width ?? DEFAULT_WIDTH,
    options.height ?? DEFAULT_HEIGHT, options.design ?? false);
  editor.root.classList.add('statistics-mpo-editor-fill');
  const column = options.column?.() ?? null;
  if (column)
    editor.setColumn?.(column);
  const sub = editor.onChanged.subscribe(() => options.onChanged?.());

  const weightInput = ui.input.float('', {value: prop.weight, min: 0, max: 1, format: '#0.000',
    onValueChanged: (v) => {
      prop.weight = Math.max(0, Math.min(1, v ?? 0));
      options.onChanged?.();
    },
  });
  weightInput.root.classList.add('statistics-mpo-weight-input');
  weightInput.captionLabel.classList.add('statistics-mpo-hidden');

  const applyPatch = (patch: Partial<PropertyDesirability>): void => {
    applyDesirabilityPatch(prop, patch);
    if ('weight' in patch) {
      weightInput.notify = false;
      weightInput.value = prop.weight;
      weightInput.notify = true;
    }
    editor.redrawAll(false);
    options.onChanged?.();
  };

  const editorHost = ui.divH([editor.root]);
  if (options.design) {
    editorHost.append(ui.icons.settings(() => new DesirabilityModeDialog(options.name?.() ?? '', prop,
      applyPatch, options.onReplaced, options.column?.() ?? null).show(), 'Settings'));
  }

  const propertyCell = ui.divV([], 'statistics-mpo-property-cell');
  const weightCell = ui.divH([weightInput.root], 'statistics-mpo-weight-cell');

  return {
    root: ui.divH([propertyCell, weightCell, editorHost], 'statistics-mpo-row'),
    propertyCell, weightCell, editorHost, editor, weightInput, applyPatch, sub,
  };
}
