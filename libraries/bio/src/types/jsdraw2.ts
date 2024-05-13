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

  export interface IPistoiaBase {
    get T(): string;
  }

  export const enum BioTypes {
    AA = 'AA',
    //BASE: 'BASE',
    ANTIBODY = 'ANTIBODY',
    PROTEIN = 'PROTEIN',
    GENE = 'GENE',
    DNA = 'DNA',
    RNA = 'RNA',
    BASE_DNA = 'BASEDNA',
    BASE_RNA = 'BASERNA'
  }

  export type BioType = `${BioTypes}`;

  export interface IBio<TBio> {
    id?: number | null;
    type: TBio;
    ambiguity?: string | null;
    annotation?: string | null;
    annotationshowright?: string | null;
  }

  export abstract class JsAtom<TBio> implements IPistoiaBase {
    get T(): 'ATOM';

    get p(): Point;

    get elem(): string;

    get bio(): IBio<TBio>

    // [propName: string]: any;

    // get charge(): number;
    //
    // get isotope(): string | null;
    //
    // get radical(): string | null;
    //
    // get group(): string | null;
    //
    // get alias(): string | null;
    //
    // get superatom(): string | null;
    //
    // get attachpoints(): [];
    //
    // get rgroup(): string | null;
    //
    // get locked(): boolean;
    //
    // get hidden(): boolean | null;
    //
    // private get _rect(): object | null;
    //
    // get color(): string | null;
    //
    // get hcount(): number | null;
    //
    // get selected(): boolean;
    //
    // get f(): number | null;
    //
    // get bonds(): [] | null;
    //
    // get id(): string | null;
    //
    // get atommapid(): string | null;
  }

  export class Atom<TBio> extends JsAtom<TBio> {
    biotype(): TBio;

    constructor(p: Point, elem: string, bio: IBio<TBio>);
  }

  export class Point {
    public x: number;
    public y: number;

    constructor(x: number, y: number);
  }
}
