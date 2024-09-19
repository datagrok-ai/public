export interface ISubstruct {
  atoms?: number[],
  bonds?: number[],
  highlightAtomColors?: { [key: number]: number[] | null },
  highlightBondColors?: { [key: number]: number[] | null }
}

export type MonomerHoverData = {
  dataFrameId: string,
  gridRowIdx: number,
  seqColName: string,
  seqPosition: number
  gridCell: any | null,
  /** Contains color of the monomer, empty lists on monomer that does not exist in molecule. */
  getSubstruct(): ISubstruct | undefined;
}

type MonomerHoverWindow = Window & {
  $monomerHover: MonomerHoverData | null;
}

declare const window: MonomerHoverWindow;

/** Return global monomer hover object. null - no monomer hover, negative seqPosition - hovered on not in a*/
export function getMonomerHover(): MonomerHoverData | null {
  return window.$monomerHover ?? null;
}

export function setMonomerHover(value: MonomerHoverData | null): void {
  window.$monomerHover = value;
}
