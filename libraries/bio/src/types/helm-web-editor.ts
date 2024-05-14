import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface IEditorOptions {
  width: number;
  height: number;
}

export interface IEditor {
  get options(): IEditorOptions;

  get div(): HTMLDivElement;
  get m(): IEditorMol;

  resize(width: number, height: number): void;
  /** Clear all contents */ clear(redraw: boolean, fireevents: boolean): void;
  /** Resets and clears undo and redo buffers */ reset(): void;
  setData(data: string, format: string): void;
}

export interface IEditorMol {
  get atoms(): IEditorMolAtom[];

  clone(selectedOnly: boolean): IEditorMol;
}

export interface IEditorMolAtom {
  get p(): IEditorPoint;

  get elem(): string;
}

export interface IEditorPoint {
  get x(): number;

  get y(): number;
}
