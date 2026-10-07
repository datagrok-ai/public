import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: ANARCI numbering (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy) of antibody and T-cell receptor chains, computed in the browser by an exact WebAssembly port of HMMER 3.4. Also returns chain, species, V/J germlines and E-values.
//input: dataframe df 
//input: column seqCol { semType: Macromolecule }
//input: string scheme = 'imgt' { choices: ["imgt","kabat","chothia","martin","aho","wolfguy"] }
//output: dataframe result
//meta.role: antibodyNumbering
//friendlyName: ANARCI (HMMER)
export async function anarciNumbering(df: DG.DataFrame, seqCol: DG.Column<any>, scheme: string) : Promise<any> {
  return await PackageFunctions.anarciNumbering(df, seqCol, scheme);
}

//description: Adds chain type, species, closest V and J germline genes with identities, and the HMMER E-value of antibody and T-cell receptor chains (ANARCI method, computed in the browser).
//input: dataframe table { caption: Table; nullable: false }
//input: column sequence { semType: Macromolecule; caption: Sequence; nullable: false }
//input: string species = 'human, mouse' { caption: Species; choices: ["human, mouse","all","human","mouse","rat","rabbit","rhesus","pig","alpaca","cow"]; description: Preferred species (ANARCI default: human and mouse) }
//friendlyName: Antibody Germlines (ANARCI)
//top-menu: Bio | Annotate | Germlines and Species (ANARCI)...
export async function antibodyGermlines(table: DG.DataFrame, sequence: DG.Column<any>, species: string) : Promise<void> {
  await PackageFunctions.antibodyGermlines(table, sequence, species);
}

//description: Searches a sequence column with a profile HMM (Pfam family, Pfam accession or HMM file) using hmmsearch in the browser; adds score, E-value and domain columns and annotates the domains.
//friendlyName: Profile HMM Search
//top-menu: Bio | Search | Profile HMM Search...
export async function profileHmmSearchDialog() : Promise<void> {
  await PackageFunctions.profileHmmSearchDialog();
}

//description: hmmsearch of a sequence column with one profile HMM. Model: HMMER text, Pfam accession(s) (PF07686) or builtin:<index>. Returns score, E-value, domains, significant and regions columns.
//input: dataframe table { caption: Table }
//input: column sequence { semType: Macromolecule; caption: Sequence }
//input: string model { caption: Model }
//input: double evalue = 10 { caption: E-value }
//input: bool gathering = false { caption: Gathering cutoff }
//output: dataframe result
export async function searchWithHmm(table: DG.DataFrame, sequence: DG.Column<any>, model: string, evalue: number, gathering: boolean) : Promise<any> {
  return await PackageFunctions.searchWithHmm(table, sequence, model, evalue, gathering);
}

//description: Annotates protein domains (built-in Pfam biologics library, Pfam accessions or an HMM file) on a sequence column using hmmscan in the browser.
//friendlyName: Find Domains (HMMER)
//top-menu: Bio | Annotate | Find Domains (HMMER)...
export async function findDomainsDialog() : Promise<void> {
  await PackageFunctions.findDomainsDialog();
}

//description: hmmscan of a sequence column against a model library: builtin, Pfam accession(s) or HMMER text. Returns one row per domain (row, model, accession, from, to, score, evalue).
//input: dataframe table { caption: Table }
//input: column sequence { semType: Macromolecule; caption: Sequence }
//input: string library = 'builtin' { caption: Library }
//input: double evalue = 0.01 { caption: E-value }
//input: bool gathering = true { caption: Gathering cutoffs }
//output: dataframe result
export async function findDomains(table: DG.DataFrame, sequence: DG.Column<any>, library: string, evalue: number, gathering: boolean) : Promise<any> {
  return await PackageFunctions.findDomains(table, sequence, library, evalue, gathering);
}

//name: Bioinformatics | ANARCI
//description: Chain, species, germlines and IMGT CDRs of an antibody or TCR sequence (ANARCI on HMMER)
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function anarciPanel(sequence: DG.SemanticValue) : any {
  return PackageFunctions.anarciPanel(sequence);
}

//name: Bioinformatics | Pfam Domains
//description: Pfam domains of a protein sequence (built-in library, gathering cutoffs; HMMER hmmscan)
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function pfamDomainsPanel(sequence: DG.SemanticValue) : any {
  return PackageFunctions.pfamDomainsPanel(sequence);
}

//description: HMMER engine status: worker count and loaded model databases.
//output: object result
export function hmmerEngineInfo() : any {
  return PackageFunctions.hmmerEngineInfo();
}
