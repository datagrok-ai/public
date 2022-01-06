export interface RDModule {
  get_mol(molString: string): RDMol;
  get_qmol(smarts: string): RDMol;
}

export interface RDMol {

  get_smiles(): string;
  get_cxsmiles(): string;
  get_smarts(): string;
  get_cxsmarts(): string;
  get_molblock(): string;
  get_v3Kmolblock(): string;
  get_inchi(): string;
  get_json(): string;
  get_svg(width: number, height: number): string;

  is_valid(): boolean;

  generate_aligned_coords(template: RDMol, useCoordGen: boolean,
                          allowOptionalAttachments: boolean, acceptFailure: boolean): string;

  get_new_coords(useCoordGen: boolean): string;

  normalize_2d_molblock(): void;
  straighten_2d_layout(): void;

  /** Reclaims the memory used for that molecule. */
  delete(): void;
}