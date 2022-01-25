export interface RDModule {
  get_mol(molString: string, options?: string): RDMol;
  get_qmol(smarts: string): RDMol;
  get_inchikey_for_inchi(input: string): string;
  version(): string;
  prefer_coordgen(prefer: boolean): void;
}

export interface RDMol {
  d_defaultWidth: number;
  d_default_Height: number;

  get_smiles(): string;
  get_cxsmiles(): string;
  get_smarts(): string;
  get_cxsmarts(): string;
  get_molblock(): string;
  get_v3Kmolblock(): string;
  get_inchi(): string;
  get_json(): string;
  get_svg(width?: number, height?: number): string;
  get_svg_with_highlightes(details: string): string;
  get_substruct_match(qmol: RDMol) : string;
  get_substruct_matches(qmol: RDMol): string;
  get_descriptors(): string;

  get_morgan_fp(radius?: number, len?: number): string;
  get_morgan_fp_as_uint8array(radius?: number, len?: number): Uint8Array;
  get_pattern_fp(len?: number): string;
  get_pattern_fp_as_uint8array(len?: number): Uint8Array;

  condense_abbreviations(maxCoverage?: number, useLinkers?: boolean): string;
  condense_abbreviations_from_defs(definitions: string, maxCoverage: number, areLinkers: boolean): string;
  generate_aligned_coords(
      template: RDMol, useCoordGen?: boolean, allowOptionalAttachments?: boolean, acceptFailure?: boolean): string;

  draw_to_canvas_with_offset(
      canvas: HTMLCanvasElement, offsetX: number, offsetY: number, width: number, height: number): string;
  draw_to_canvas(canvas: HTMLCanvasElement, width: number, height: number): string;
  draw_to_canvas_with_highlights(canvas: HTMLCanvasElement, details: string): string;
  get_morgan_fp_as_uint8array(radius: number, fplen: number): Uint8Array;

  is_valid(): boolean;

  get_stereo_tags(): string;
  get_aromatic_form(): string;
  get_kekule_from(): string;
  get_new_coords(useCoordGen?: boolean): string;
  remove_hs(): string;
  add_hs(): string;

  normalize_2d_molblock(): string;
  straighten_2d_layout(): void;

  /** Reclaims the memory used for that molecule. */
  delete(): void;
}

export interface SubstructLibrary {
  d_num_bits: number;
  d_defaultNumBits: number;
  d_defaultUseChirality: boolean;
  d_defaultNumThreads: number;
  d_defaultMaxResults: number;

  add_mol(mol: RDMol): number;
  add_smiles(smiles: string): number;
  add_trusted_smiles(smiles: string): number;
  get_mol(i: number): RDMol;
  get_matches(qmol: RDMol, useChirality?: boolean, numThreads?: number, maxResults?: number): string;
  count_matches(qmol: RDMol, useChirality?: boolean, numThreads?: number): number;
}

