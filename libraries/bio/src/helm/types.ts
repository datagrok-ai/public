import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export interface IWebEditorSizes {
  get rightwidth(): number;
  get topheight(): number;
  get bottomheight(): number;
}

export interface IWebEditorHelm {
  setSequence(seq: string, format: string,
    sugar?: string, linker?: string, append?: boolean, separator?: string): boolean;
}

export interface IWebEditorCanvas {
  get helm(): IWebEditorHelm;
  getHelm(ret: boolean): string;

  resize(width: number, height: number): void;
}

export interface IWebEditorApp {
  get toolbarheight(): number;
  get canvas(): IWebEditorCanvas;
  get structureview(): any;
  get mex(): any;

  calculateSizes(): IWebEditorSizes;
}

export interface IHelmWebEditorEditor {
  setHelm(helm: string): void;
}

export interface IHelmWebEditor {
  get editor(): IHelmWebEditorEditor;
  get host(): HTMLDivElement;

  resizeEditor(width: number, height: number): void;
}
