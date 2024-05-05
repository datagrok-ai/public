declare module 'JSDraw2' {
  export interface IEditorOptions {
    width: number;
    height: number;
    skin: string;
    viewonly: boolean;
  }

  export interface IEditorPoint {
    get x(): number;

    get y(): number;
  }

  export interface IAtomBio {
    type: 'HELM_BASE' | 'HELM_SUGAR' | 'HELM_LINKER' | 'HELM_AA' | 'HELM_CHEM' | 'HELM_BLOB' |
      /* type as in Pistoia.HELM-uncompressed.js */'HELM_NUCLETIDE';

  }

  export interface IEditorMolAtom {
    get p(): IEditorPoint;

    get elem(): string;

    get bio(): IAtomBio;
  }

  export interface IEditorMolBond {
    get a1(): IEditorMolAtom;
    get a2(): IEditorMolAtom;

    get type(): number;
  }

  export interface IEditorMol {
    get atoms(): IEditorMolAtom[];
    get bonds(): IEditorMolBond[];

    clone(selectedOnly: boolean): IEditorMol;
  }

  export class Editor {
    get options(): IEditorOptions;

    get div(): HTMLDivElement;

    get m(): IEditorMol;

    resize(width: number, height: number): void;

    /** Clear all contents */ clear(redraw: boolean, fireevents: boolean): void;

    /** Resets and clears undo and redo buffers */ reset(): void;

    setData(data: string, format: string): void;

    setHelm(helm: string): void;

    setSize(w: number, h: number): void;

    getMolfile(): string

    constructor(host: HTMLElement, options?: Partial<IEditorOptions>);
  }
}
