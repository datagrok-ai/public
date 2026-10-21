/**
 * What Chem calls in the one Crux file (vendor/chem-wasm) beyond RDKit.js's API: the search Collection (the API of the
 * search-only build it replaces), descriptors and the scaffold tree. Rendering uses the RDKit.js-shaped calls
 * (`RDModule`, chem-meta).
 */

export interface CruxCollection {
  /** The matching molecules' input positions (0 = no limit). Throws on a SMARTS crux cannot read. */
  substructureSearch(smarts: string, limit: number): Uint32Array;
  free(): void;
}

export interface CruxCollectionBuilder {
  /** Read what RDKit's default read rejects the way RDKit reads it without Kekulize (Chem's getMolSafe retry). */
  setLenient(lenient: boolean): void;
  /** The number of molecules added. */
  addMany(smiles: string[]): number;
  finish(): CruxCollection;
}

/** A descriptor's values as the chem-chem container returns them: Python's type of its first value, nulls kept. */
export interface DescriptorValues {
  type: 'int' | 'double';
  value: (number | null)[];
}

export interface CruxWorkerApi {
  CollectionBuilder: new (buildIndex: boolean) => CruxCollectionBuilder;
  /** The container's `GET /chem/descriptors/tree` body (JSON text): groups by RDKit module, the 3D group too. */
  descriptorTree(): string;
  /**
   * The container's `POST /chem/descriptors` body (JSON text, `{name: DescriptorValues}`): molecules as chem.py reads
   * them (null or unreadable: nulls), descriptor or group names, in request order.
   */
  calcDescriptors(molecules: (string | null)[], names: string[]): string;
  /** scripts/scaffold_tree.py's JSON text (ScaffoldGraph's Scaffold Tree); throws where the script fails. */
  scaffoldTree(smiles: string[], ringCutoff: number, dischargeAndDeradicalize: boolean): string;
}
