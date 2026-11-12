import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: ANARCI numbering (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy) of antibody and T-cell receptor chains, computed in the browser. Also returns chain, species, V/J germlines and E-values. ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.
//input: dataframe df 
//input: column seqCol { semType: Macromolecule }
//input: string scheme = 'imgt' { choices: ["imgt","kabat","chothia","martin","aho","wolfguy"] }
//input: string species = 'human and mouse' { optional: true; caption: Species; choices: ["human and mouse","all","human","mouse","rat","rabbit","rhesus","pig","alpaca","cow"]; description: Preferred species of the germline models (ANARCI default: human and mouse, other species only when they score clearly better) }
//input: double bitScoreThreshold = 80 { optional: true; caption: Bit score threshold; description: Minimum HMMER bit score of a domain (ANARCI default 80) }
//input: string domain = 'first' { optional: true; caption: Domain; choices: ["first","heavy","light"]; description: Which domain of multi-domain chains (scFv, diabodies) to number: the N-terminal one, the first heavy (TCR beta/delta) or the first light (TCR alpha/gamma) domain }
//output: dataframe result
//meta.role: antibodyNumbering
//friendlyName: ANARCI (HMMER)
export async function anarciNumbering(df: DG.DataFrame, seqCol: DG.Column<any>, scheme: string, species: string, bitScoreThreshold: number, domain: string) : Promise<any> {
  return await PackageFunctions.anarciNumbering(df, seqCol, scheme, species, bitScoreThreshold, domain);
}

//description: Adds chain type, species, closest V and J germline genes with identities, and the HMMER E-value of antibody and T-cell receptor chains (ANARCI method, computed in the browser). ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.
//input: dataframe table { caption: Table; nullable: false }
//input: column sequence { semType: Macromolecule; caption: Sequence; nullable: false }
//input: string species = 'human and mouse' { caption: Species; choices: ["human and mouse","all","human","mouse","rat","rabbit","rhesus","pig","alpaca","cow"]; description: Preferred species (ANARCI default: human and mouse) }
//friendlyName: Antibody Germlines (ANARCI)
//top-menu: Bio | Annotate | Germlines and Species (ANARCI)...
export async function antibodyGermlines(table: DG.DataFrame, sequence: DG.Column<any>, species: string) : Promise<void> {
  await PackageFunctions.antibodyGermlines(table, sequence, species);
}

//description: Searches a sequence column with a profile HMM (Pfam family, Pfam accession or HMM file) using hmmsearch in the browser; adds score, E-value and domain columns and annotates the domains. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//friendlyName: Profile HMM Search
//top-menu: Bio | Search | Profile HMM Search...
//meta.requiresSemType: Macromolecule
export async function profileHmmSearchDialog() : Promise<void> {
  await PackageFunctions.profileHmmSearchDialog();
}

//description: hmmsearch of a sequence column with one profile HMM. Model: HMMER text, Pfam accession(s) (PF07686) or builtin:<index>. Returns score, E-value, domains, significant and regions columns. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: dataframe table { caption: Table }
//input: column sequence { semType: Macromolecule; caption: Sequence }
//input: string model { caption: Model }
//input: double evalue = 10 { caption: E-value }
//input: bool gathering = false { caption: Gathering cutoff }
//output: dataframe result
export async function searchWithHmm(table: DG.DataFrame, sequence: DG.Column<any>, model: string, evalue: number, gathering: boolean) : Promise<any> {
  return await PackageFunctions.searchWithHmm(table, sequence, model, evalue, gathering);
}

//description: Aligns sequences to one profile HMM with hmmalign (HMMER 3.4, in the browser): the best-matching Pfam family of the built-in library, a Pfam family or accession, the best antibody/TCR germline model (ANARCI) or an HMM file. Columns are model positions and insertions between them. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: column sequences { semType: Macromolecule }
//input: string profile = 'Best Pfam family (auto)' { optional: true; caption: Profile; choices: ["Best Pfam family (auto)","Pfam family or accession","Antibody/TCR germline (auto)","HMM file"]; description: The profile HMM to align to }
//input: string family { optional: true; nullable: true; caption: Family; description: Pfam family of the built-in library (V-set) or any Pfam accession (PF07686), for Profile = Pfam family or accession }
//input: file hmmFile { optional: true; nullable: true; caption: HMM file; description: A file with one profile HMM (HMMER2/3 text or .h3m), for Profile = HMM file }
//input: bool trim = false { optional: true; caption: Trim; description: Trim residues outside the model's alignment (hmmalign --trim) }
//output: column result
//meta.role: sequenceMSA
//meta.acceptsCanonical: true
//friendlyName: Profile HMM (hmmalign)
export async function alignToProfile(sequences: DG.Column<any>, profile: string, family: string, hmmFile: any, trim: boolean) : Promise<any> {
  return await PackageFunctions.alignToProfile(sequences, profile, family, hmmFile, trim);
}

//description: Profiles for Bio's WebLogo (its Profile property): hmmlogo heights in bits of a Pfam family of the built-in library (V-set), a Pfam accession (PF07686), an ANARCI germline model (germline:human_H) or a profile built in the table. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: dataframe table { nullable: true }
//input: string column { nullable: true }
//input: string profile 
//output: object result
//meta.role: webLogoProfile
//friendlyName: Profile HMM Logo Heights
export async function hmmWebLogoProfile(table: any, column: string, profile: string) : Promise<any> {
  return await PackageFunctions.hmmWebLogoProfile(table, column, profile);
}

//description: Adds a WebLogo of a profile HMM (hmmlogo heights) to the current table view. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//friendlyName: Profile HMM Logo
//top-menu: Bio | Analyze | Profile HMM Logo...
//meta.requiresSemType: Macromolecule
export async function profileLogoDialog() : Promise<void> {
  await PackageFunctions.profileLogoDialog();
}

//description: Builds a profile HMM (hmmbuild, in the browser) from aligned rows of a column, then searches the table with it, annotates the hits and shows its logo. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//friendlyName: Build Profile HMM
//top-menu: Bio | Search | Build Profile HMM...
//meta.requiresSemType: Macromolecule
export async function buildProfileDialog() : Promise<void> {
  await PackageFunctions.buildProfileDialog();
}

//description: hmmbuild of the rows of an aligned sequence column (all, or the selected ones). Returns the HMMER3 model text, identical to C hmmbuild -n name; search with it via searchWithHmm. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: dataframe table { caption: Table }
//input: column alignment { semType: Macromolecule; caption: Alignment }
//input: string name { caption: Name }
//input: bool selectedOnly = false { optional: true; caption: Selected rows only }
//output: string result
export async function buildProfile(table: DG.DataFrame, alignment: DG.Column<any>, name: string, selectedOnly: boolean) : Promise<string> {
  return await PackageFunctions.buildProfile(table, alignment, name, selectedOnly);
}

//description: Searches a protein sequence column with one query sequence (a row or typed) using phmmer in the browser; adds bit score, E-value and domain columns and annotates the matches. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//friendlyName: Homology Search (phmmer)
//top-menu: Bio | Search | Homology Search (phmmer)...
//meta.requiresSemType: Macromolecule
export async function homologySearchDialog() : Promise<void> {
  await PackageFunctions.homologySearchDialog();
}

//description: phmmer of one query protein sequence against a sequence column (Z = non-empty rows). Returns score, E-value, domains, significant and regions columns. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: dataframe table { caption: Table }
//input: column sequence { semType: Macromolecule; caption: Sequence }
//input: string query { caption: Query }
//input: double evalue = 10 { optional: true; caption: E-value }
//input: string matrix = 'BLOSUM62' { optional: true; caption: Matrix; choices: ["BLOSUM62","BLOSUM45","BLOSUM50","BLOSUM80","BLOSUM90","PAM30","PAM70","PAM120","PAM240"] }
//output: dataframe result
export async function homologySearch(table: DG.DataFrame, sequence: DG.Column<any>, query: string, evalue: number, matrix: string) : Promise<any> {
  return await PackageFunctions.homologySearch(table, sequence, query, evalue, matrix);
}

//description: Annotates protein domains (built-in Pfam biologics library, Pfam accessions or an HMM file) on a sequence column using hmmscan in the browser. Overlapping matches of one Pfam clan are resolved as Pfam does (this package's implementation, checked against pfam_scan.pl). HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//friendlyName: Find Domains (HMMER)
//top-menu: Bio | Annotate | Find Domains (HMMER)...
//meta.requiresSemType: Macromolecule
export async function findDomainsDialog() : Promise<void> {
  await PackageFunctions.findDomainsDialog();
}

//description: hmmscan of a sequence column against a model library: builtin, Pfam accession(s) or HMMER text. Returns one row per domain (row, model, accession, from, to, score, evalue, clan); overlapping domains of a Pfam clan are resolved as Pfam does. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//input: dataframe table { caption: Table }
//input: column sequence { semType: Macromolecule; caption: Sequence }
//input: string library = 'builtin' { caption: Library }
//input: double evalue = 0.01 { caption: E-value }
//input: bool gathering = true { caption: Gathering cutoffs }
//input: bool clanCompetition = true { optional: true; caption: Clan competition; description: Of overlapping domains from the same Pfam clan keep the most significant (Pfam families only) }
//output: dataframe result
export async function findDomains(table: DG.DataFrame, sequence: DG.Column<any>, library: string, evalue: number, gathering: boolean, clanCompetition: boolean) : Promise<any> {
  return await PackageFunctions.findDomains(table, sequence, library, evalue, gathering, clanCompetition);
}

//name: Bioinformatics | ANARCI
//description: Chain, species, germlines and IMGT CDRs of an antibody or TCR sequence (ANARCI method). ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule; units: fasta }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function anarciPanel(sequence: DG.SemanticValue) : any {
  return PackageFunctions.anarciPanel(sequence);
}

//name: Bioinformatics | Pfam Domains
//description: Pfam domains of a protein sequence (built-in library, gathering cutoffs; hmmscan). HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule; units: fasta }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function pfamDomainsPanel(sequence: DG.SemanticValue) : any {
  return PackageFunctions.pfamDomainsPanel(sequence);
}

//description: Engine status (Rusty-HMMER, a Rust/WebAssembly reimplementation of HMMER 3.4): worker count and loaded model databases.
//output: object result
export function hmmerEngineInfo() : any {
  return PackageFunctions.hmmerEngineInfo();
}
