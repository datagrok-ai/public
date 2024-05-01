/* Do not import anything, otherwise get error TS2664: Invalid module name in augmentation */
declare module 'scil' {
  function apply<T>(dest: T, atts: Partial<T>, defaults?: Partial<T>): void;

  namespace Utils {
    export function alert(s: string): void;

    export function isNullOrEmpty(s: string): boolean;

    export function endswith(s: string, token: string, casesensitive?: boolean): boolean;
  }
}

declare module 'org' {
  namespace helm {
    export interface IAtom {

    }

    export interface IMonomer {

    }

    export interface IAppOptions {
      showabout: boolean;
      mexfontsize: string;
      mexrnapinontab: boolean;
      topmargin: number;
      mexmonomerstab: boolean;
      sequenceviewonly: boolean;
      mexfavoritefirst: boolean;
      mexfilter: boolean;
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

    export interface IWebEditorSizes {
      get rightwidth(): number;
      get topheight(): number;
      get bottomheight(): number;
    }

    export interface IWebEditorAppProperties {
      parent: HTMLElement;
    }

    export interface IWebEditorApp {
      notation: HTMLElement;
      sequence: HTMLElement;
      properties: IWebEditorAppProperties;
      get toolbarheight(): number;
      get canvas(): IWebEditorCanvas;
      get structureview(): any;
      get mex(): any;

      calculateSizes(): IWebEditorSizes;
    }

    export interface IApp {
      new(host: HTMLDivElement, options: IAppOptions): IWebEditorApp;
    }

    export interface IMolViewer {
      molscale: number;
    }

    export interface IMonomers {
      addOneMonomer(monomer: IMonomer): void;
      getMonomer(biotype: string, elem: string): IAtom;
      getMonomerSet(biotype: string): any;
      clear(): void;
    }

    export interface IOrgHelmWebEditor {
      App: IApp;
      Monomers: IMonomers;
      MolViewer: IMolViewer;
      kCaseSensitive: boolean;

      monomerTypeList(): string[];
    }

    export const webeditor: IOrgHelmWebEditor;
  }
}
