// hwe migration: the HELM value vocabulary and the editor graph/view shapes are
// owned by `@datagrok-libraries/hwe` (the standalone HELM editor). This barrel
// re-exports them under the legacy type names the platform still uses, and keeps
// the Datagrok-flavored monomer interfaces defined locally (rebuilt on hwe's
// value types). As a result `@datagrok-libraries/bio` no longer depends on the
// legacy `@datagrok-libraries/js-draw-lite` / `@datagrok-libraries/helm-web-editor`.
import type {
  HelmType, PolymerType, MonomerType,
  HelmAtomLike, HelmBondLike, HelmMolLike, HelmAtomBioLike,
  AppLike, JSDraw2EditorLike, MonomerExplorerLike, MonomersFuncsLike, DrawOptionsLike,
} from '@datagrok-libraries/hwe';

export type {HelmType, PolymerType, MonomerType};

/** 2-D point (the legacy JSDraw2 `Point` shape). */
export interface Point {
  x: number;
  y: number;
}

// ---- Legacy JSDraw2 / HELM graph shapes → hwe view types -------------------
// The platform refers to these JSDraw2-era names; at runtime every value is the
// hwe immutable model exposed through the adapter's legacy-shaped view
// (`buildHelmMolView`, `LegacyEditorWrapper`, `LegacyAppWrapper`). Phantom
// generic parameters preserve the legacy call sites (e.g.
// `Mol<HelmType, IHelmBio>`, `Editor<HelmType, IHelmBio, IHelmEditorOptions>`).
export type IBio<_TBio = unknown> = HelmAtomBioLike;
export type IHelmBio = HelmAtomBioLike;
export type Atom<_T = unknown, _TBio = unknown> = HelmAtomLike;
export type Bond<_T = unknown, _TBio = unknown> = HelmBondLike;
export type Mol<_T = unknown, _TBio = unknown> = HelmMolLike;
export type Editor<_T = unknown, _TBio = unknown, _O = unknown> = JSDraw2EditorLike;
export type IJsAtom = HelmAtomLike;

export type HelmAtom = HelmAtomLike;
export type HelmBond = HelmBondLike;
export type HelmMol = HelmMolLike;
export type HelmEditor = JSDraw2EditorLike;
export type HelmString = string;

export type App = AppLike;
export type MonomerExplorer = MonomerExplorerLike;
export type MonomersFuncs = MonomersFuncsLike;

// ---- Datagrok-flavored monomer interfaces (Pistoia HELMmonomer-shaped) ------
export type WebEditorRGroups = { [group: string]: string };

export interface IMonomerColors {
  linecolor: string;
  backgroundcolor: string;
  textcolor: string;
  nature?: string;
}

export interface IMonomer {
  /** symbol */ id?: string;
  /** name */ n?: string;
  /** natural analog */ na?: string;
  /** polymer type */ type?: PolymerType;
  /** monomer type */ mt?: MonomerType;
  /** molfile */ m?: string;
  /** molfile compressed */ mz?: any;
  /** substituents */ at?: WebEditorRGroups;
  /** number of substituents */ rs?: number;
  issmiles?: boolean;
  smiles?: string;
  name?: string;
  oldname?: string;
  /** Used by Formula */ stats?: any;
}

export interface IWebEditorMonomer extends IMonomer, Partial<IMonomerColors> {}

export type MonomerSetType = { [symbol: string]: IMonomer };

export type GetMonomerResType = IWebEditorMonomer | null;
export type GetMonomerFunc = (a: HelmAtom | HelmType, name?: string) => GetMonomerResType;

export interface IHelmDrawOptions extends DrawOptionsLike {
  getMonomer?: GetMonomerFunc;
}

export interface IHelmEditorOptions {
  width?: number;
  height?: number;
  viewonly?: boolean;
}

export interface ISeqMonomer {
  position: number;
  biotype: HelmType;
  symbol: string;
}

export interface IHelmWebEditor {
  get editor(): HelmEditor;
  get host(): HTMLDivElement;

  resizeEditor(width: number, height: number): void;
}
