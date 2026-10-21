// Node environment for the crux verification scripts: Chem's RDKit build, a crux wasm build, and Chem's own
// query code (utils/mol-creation_rdkit.ts) transpiled on the fly, so the scripts check the code the package runs, not a
// copy of it.
import {createRequire} from 'module';
import {existsSync, mkdtempSync, readFileSync, writeFileSync} from 'fs';
import {tmpdir} from 'os';
import {join} from 'path';
import {fileURLToPath, pathToFileURL} from 'url';

export const CHEM_DIR = fileURLToPath(new URL('../../', import.meta.url));
/** The reddata checkout's data submodule (demo datasets). */
export const DATA_DIR = join(CHEM_DIR, '..', '..', '..', 'data');
export const requireFromChem = createRequire(join(CHEM_DIR, 'package.json'));
const SRC = join(CHEM_DIR, 'src');
const tmp = mkdtempSync(join(tmpdir(), 'chem-crux-'));

export async function loadRdkit() {
  const version = requireFromChem(join(SRC, 'rdkit_lib_version.js'));
  const init = requireFromChem(join(SRC, 'RDKit_minimal.js'));
  const rdkit = await init({locateFile: () => join(SRC, `${version}.wasm`)});
  rdkit.use_legacy_stereo_perception(false);
  return rdkit;
}

/** A crux build (`--target web`): the out dir in CRUX_WASM (a chem-wasm or a crux-wasm build), or Chem's one Crux file
 * (vendor/chem-wasm). */
export async function loadCrux() {
  const dir = process.env.CRUX_WASM ?? join(CHEM_DIR, 'vendor', 'chem-wasm');
  const name = existsSync(join(dir, 'chem_wasm.js')) ? 'chem_wasm' : 'crux_wasm';
  // the glue is an ES module and Chem's package.json does not declare one, so it is imported from an .mjs copy
  const glue = join(tmp, `${name}.mjs`);
  writeFileSync(glue, readFileSync(join(dir, `${name}.js`)));
  const crux = await import(pathToFileURL(glue).href);
  const wasm = await crux.default({module_or_path: readFileSync(join(dir, `${name}_bg.wasm`))});
  return {crux, wasm, dir};
}

/** getMolSafe / getQueryMolSafe from utils/mol-creation_rdkit.ts (Chem's search builds the query with getQueryMolSafe,
 * on RDKit and on the Crux module). */
export async function loadChemQueryCode() {
  const chemConstants = readFileSync(join(SRC, 'utils', 'chem-constants.ts'), 'utf8');
  const maxSmilesLength = chemConstants.match(/MAX_SMILES_LENGTH = (\d+)/)[1];
  const molCreation = await importChemModule('utils/mol-creation_rdkit.ts', [
    [/import \{hasNewLines, isMolBlock\} from '\.\/chem-common';/,
      `const hasNewLines = (s) => s.includes('\\n') || s.includes('\\r');
       const isMolBlock = (s) => s.includes('M  END');`],
    [/import \{MAX_SMILES_LENGTH\} from '\.\/chem-constants';/, `const MAX_SMILES_LENGTH = ${maxSmilesLength};`],
    // only reached for molblocks by isFragment / _isSmarts, which the search does not call on molblocks
    [/import \{MolfileHandler\} from .*;/, 'const MolfileHandler = undefined;'],
  ]);
  return molCreation;
}

async function importChemModule(path, replacements) {
  const ts = requireFromChem('typescript');
  let code = readFileSync(join(SRC, path), 'utf8');
  for (const [from, to] of replacements) {
    if (!from.test(code))
      throw new Error(`${path} changed: no match for ${from} (update misc/crux/node-env.mjs)`);
    code = code.replace(from, to);
  }
  const js = ts.transpileModule(code,
    {compilerOptions: {module: ts.ModuleKind.ESNext, target: ts.ScriptTarget.ES2022}}).outputText;
  const file = join(tmp, path.replace(/[\\/]/g, '_').replace(/\.ts$/, '.mjs'));
  writeFileSync(file, js);
  return import(pathToFileURL(file).href);
}

/**
 * First column whose header mentions smiles (else the first one) of a CSV with one molecule per line, or the first
 * field of each line of a headerless .smi file (crux-bench's corpora).
 */
export function readSmilesCsv(path, limit = Infinity) {
  const lines = readFileSync(path, 'utf8').split(/\r?\n/);
  if (path.endsWith('.smi'))
    return lines.filter((l) => l.length).slice(0, limit).map((l) => l.split(/[\t ]/)[0]);
  const header = lines[0].split(',');
  const col = Math.max(0, header.findIndex((h) => /smiles/i.test(h)));
  return lines.slice(1).filter((l) => l.length).slice(0, limit).map((l) => l.split(',')[col]);
}
