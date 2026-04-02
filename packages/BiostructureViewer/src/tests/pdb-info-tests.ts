import * as grok from 'datagrok-api/grok';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {parsePdbHeaders, PdbHeaderInfo} from '../utils/pdb-helper';
import {RcsbGraphQLAdapter} from '../utils/rcsb-gql-adapter';
import {_package} from '../package-test';


// -- Sample PDB text for testing (RCSB-style with full headers) --
const RCSB_PDB_SNIPPET = `\
HEADER    HYDROLASE (ACID PROTEINASE)             25-APR-97   1QBS
TITLE     HIV-1 PROTEASE INHIBITORS WITH LOW NANOMOLAR POTENCY
COMPND    MOL_ID: 1;
COMPND   2 MOLECULE: HIV-1 PROTEASE;
COMPND   3 CHAIN: A, B;
COMPND   4 EC: 3.4.23.16;
COMPND   5 ENGINEERED: YES
SOURCE    MOL_ID: 1;
SOURCE   2 ORGANISM_SCIENTIFIC: HUMAN IMMUNODEFICIENCY VIRUS 1;
SOURCE   3 EXPRESSION_SYSTEM: ESCHERICHIA COLI;
KEYWDS    HYDROLASE (ACID PROTEINASE)
EXPDTA    X-RAY DIFFRACTION
REMARK   2
REMARK   2 RESOLUTION.    1.80 ANGSTROMS.
REMARK   3   PROGRAM     : X-PLOR
REMARK   3   R VALUE            (WORKING SET) : 0.189
JRNL        TITL   CYCLIC HIV PROTEASE INHIBITORS
JRNL        REF    J.MED.CHEM.                   V.  39  3514 1996
JRNL        DOI    10.1021/JM9602571
JRNL        PMID   8784449
HET    DMP  A 323      42
HETNAM     DMP SOME INHIBITOR COMPOUND
FORMUL   3  DMP    C28 H32 N2 O6
MODRES 1QBS CSO A   67  CYS  S-HYDROXYCYSTEINE
SITE     1 AC1 21 ASP A  25  GLY A  27  ALA A  28  ASP A  29
SITE     2 AC1 21 ASP A  30  ILE A  47  GLY A  48  GLY A  49
SEQADV 1QBS TRP A  243  UNP  A0A7Z6MKR LEU   243 ENGINEERED MUTATION
SSBOND   1 CYS A   36    CYS A   49
SSBOND   2 CYS A  181    CYS A  187
HELIX    1   1 ARG A   87  ILE A   93  1
HELIX    2   2 ARG B   87  LEU B   90  1
SHEET    1   A 5 GLN A   2  THR A   4  0
SHEET    2   A 5 LEU A  10  ILE A  15  0
REMARK 465
REMARK 465 MISSING RESIDUES
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
REMARK 465
REMARK 465   M RES C SSSEQI
REMARK 465     ARG A    23
REMARK 465     SER A    24
REMARK 465     PRO A    25
ATOM      1  N   PRO A   1       2.000   3.000   4.000  1.00 20.00           N
ATOM      2  CA  PRO A   1       3.000   4.000   5.000  1.00 20.00           C
ATOM      3  N   GLN A   2       4.000   5.000   6.000  1.00 20.00           N
ATOM      4  N   PRO B   1       5.000   6.000   7.000  1.00 20.00           N
END
`;

// -- Minimal PDB (MOE/Schrödinger-style, no headers) --
const MINIMAL_PDB_SNIPPET = `\
HEADER    RNA BINDING PROTEIN                                 denovo_1
REMARK  99
REMARK  99 MOE v2022.02 (Chemical Computing Group ULC) Wed Dec 04 02:59:53 2024
SSBOND   1 CYS A   36    CYS A   49
ATOM      1  N   PRO A   1       2.000   3.000   4.000  1.00 20.00           N
ATOM      2  CA  PRO A   1       3.000   4.000   5.000  1.00 20.00           C
ATOM      3  N   GLN A   2       4.000   5.000   6.000  1.00 20.00           N
ATOM      4  N   PRO B   1       5.000   6.000   7.000  1.00 20.00           N
END
`;


category('PDB Header Parsing', () => {
  // -- Basic fields --
  test('parsePdbHeaders: PDB ID', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.pdbId, '1QBS');
  });

  test('parsePdbHeaders: title', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.title, 'HIV-1 PROTEASE INHIBITORS WITH LOW NANOMOLAR POTENCY');
  });

  test('parsePdbHeaders: classification', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.classification, 'HYDROLASE (ACID PROTEINASE)');
  });

  test('parsePdbHeaders: deposit date', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.depositDate, '25-APR-97');
  });

  test('parsePdbHeaders: method', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.method, 'X-RAY DIFFRACTION');
  });

  test('parsePdbHeaders: resolution', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.resolution, 1.80);
  });

  test('parsePdbHeaders: R-factor', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.rFactor, 0.189);
  });

  // -- Organism & expression system --
  test('parsePdbHeaders: organism', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.organism, 'HUMAN IMMUNODEFICIENCY VIRUS 1');
  });

  test('parsePdbHeaders: expression system', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.expressionSystem, 'ESCHERICHIA COLI');
  });

  // -- Entities --
  test('parsePdbHeaders: entities', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.entities != null, true);
    expect(info.entities!.length, 1);
    expect(info.entities![0].molId, '1');
    expect(info.entities![0].molecule, 'HIV-1 PROTEASE');
    expect(info.entities![0].chains!.join(','), 'A,B');
    expect(info.entities![0].engineered === true, true);
  });

  // -- Citation --
  test('parsePdbHeaders: citation', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.citation != null, true);
    expect(info.citation!.title, 'CYCLIC HIV PROTEASE INHIBITORS');
    expect(info.citation!.doi, '10.1021/JM9602571');
    expect(info.citation!.pmid, 8784449);
  });

  // -- Ligands --
  test('parsePdbHeaders: ligands (excludes modified residues)', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.ligands != null, true);
    expect(info.ligands!.length, 1);
    expect(info.ligands![0].id, 'DMP');
    expect(info.ligands![0].name, 'SOME INHIBITOR COMPOUND');
  });

  // -- Modified residues --
  test('parsePdbHeaders: modified residues', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.modifiedResidues != null, true);
    expect(info.modifiedResidues!.length, 1);
    expect(info.modifiedResidues![0].resName, 'CSO');
    expect(info.modifiedResidues![0].standardRes, 'CYS');
  });

  // -- Secondary structure --
  test('parsePdbHeaders: secondary structure', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.secondaryStructure != null, true);
    expect(info.secondaryStructure!.helices.length, 2);
    expect(info.secondaryStructure!.sheets.length, 2);
  });

  // -- Sites --
  test('parsePdbHeaders: binding sites', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.sites != null, true);
    expect(info.sites!.length, 1);
    expect(info.sites![0].name, 'AC1');
    expect(info.sites![0].residues.length, 8);
  });

  // -- Missing residues --
  test('parsePdbHeaders: missing residues', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.missingResidues != null, true);
    expect(info.missingResidues!.length, 3);
    expect(info.missingResidues![0].resName, 'ARG');
    expect(info.missingResidues![0].resSeq, 23);
  });

  // -- Mutations --
  test('parsePdbHeaders: mutations', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.mutations != null, true);
    expect(info.mutations!.length, 1);
    expect(info.mutations![0].chain, 'A');
  });

  // -- Disulfide bonds --
  test('parsePdbHeaders: disulfide bonds', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.disulfideBondCount, 2);
  });

  // -- Atom count and chains --
  test('parsePdbHeaders: atom count and chains', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.atomCount, 4);
    expect(info.chains!.join(','), 'A,B');
  });

  // -- Residue count by chain --
  test('parsePdbHeaders: residue count by chain', async () => {
    const info = parsePdbHeaders(RCSB_PDB_SNIPPET);
    expect(info.residueCountByChain != null, true);
    expect(info.residueCountByChain!['A'], 2); // PRO 1, GLN 2
    expect(info.residueCountByChain!['B'], 1); // PRO 1
  });
});


category('PDB Header Parsing: Minimal PDB', () => {
  test('parsePdbHeaders: MOE file has no valid PDB ID', async () => {
    const info = parsePdbHeaders(MINIMAL_PDB_SNIPPET);
    expect(info.pdbId == null, true);
  });

  test('parsePdbHeaders: MOE file has classification', async () => {
    const info = parsePdbHeaders(MINIMAL_PDB_SNIPPET);
    expect(info.classification, 'RNA BINDING PROTEIN');
  });

  test('parsePdbHeaders: MOE file has software remarks', async () => {
    const info = parsePdbHeaders(MINIMAL_PDB_SNIPPET);
    expect(info.softwareRemarks != null, true);
    expect(info.softwareRemarks!.length > 0, true);
  });

  test('parsePdbHeaders: MOE file has no title/method/citation', async () => {
    const info = parsePdbHeaders(MINIMAL_PDB_SNIPPET);
    expect(info.title == null, true);
    expect(info.method == null, true);
    expect(info.citation == null, true);
  });

  test('parsePdbHeaders: MOE file still has atom count and chains', async () => {
    const info = parsePdbHeaders(MINIMAL_PDB_SNIPPET);
    expect(info.atomCount, 4);
    expect(info.chains!.join(','), 'A,B');
    expect(info.disulfideBondCount, 1);
  });
});


category('PDB Header Parsing: Sample file', () => {
  test('parsePdbHeaders: 1bdq.pdb from samples', async () => {
    const pdbText = await grok.dapi.files.readAsText(
      `System:AppData/${_package.name}/samples/1bdq.pdb`);
    const info = parsePdbHeaders(pdbText);
    // 1bdq.pdb is minimal - no HEADER line, but has HELIX/SHEET and ATOM records
    expect(info.atomCount! > 0, true);
    expect(info.chains != null, true);
    expect(info.secondaryStructure != null, true);
    expect(info.secondaryStructure!.helices.length, 2);
    expect(info.secondaryStructure!.sheets.length > 0, true);
  });
});


category('RCSB GraphQL Adapter', () => {
  test('getEntryInfo: fetches basic info for 1QBS', async () => {
    const info = await RcsbGraphQLAdapter.getEntryInfo('1QBS');
    expect(info.id, '1QBS');
    expect(info.title != null, true);
    expect(info.resolution != null, true);
    expect(info.experimentalMethod != null, true);
    expect(info.classification != null, true);
  });

  test('getEntryInfo: fetches polymer entities', async () => {
    const info = await RcsbGraphQLAdapter.getEntryInfo('1QBS');
    expect(info.polymerEntities != null, true);
    expect(info.polymerEntities!.length > 0, true);
    expect(info.polymerEntities![0].entityId != null, true);
    expect(info.polymerEntities![0].chains.length > 0, true);
  });

  test('getEntryInfo: fetches nonpolymer entities', async () => {
    const info = await RcsbGraphQLAdapter.getEntryInfo('4UJ1');
    expect(info.nonpolymerEntities != null, true);
    expect(info.nonpolymerEntities!.length > 0, true);
    const ne = info.nonpolymerEntities![0];
    expect(ne.compId != null, true);
    expect(ne.smiles != null, true);
  });

  test('getEntryInfo: fetches citation', async () => {
    const info = await RcsbGraphQLAdapter.getEntryInfo('4UJ1');
    expect(info.citation != null, true);
    expect(info.citation!.title != null, true);
    expect(info.citation!.pubmedId != null, true);
  });

  test('getEntryInfo: fetches binding affinity when available', async () => {
    // 3HVG has known binding affinity data in BindingDB
    const info = await RcsbGraphQLAdapter.getEntryInfo('3HVG');
    expect(info.bindingAffinities != null, true);
    expect(info.bindingAffinities!.length > 0, true);
    expect(info.bindingAffinities![0].value > 0, true);
  });
});
