import * as DG from 'datagrok-api/dg';
//@ts-ignore
import initBBWasm, {WasmBitBirch} from '../bbwasm.js';
import {_package} from '../package';
import {Fingerprint, defaultMorganFpLength} from '../utils/chem-common';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {getRDKitFpAsUint8Array} from '../chem-searches';

let bbwasmInitialized = false;

async function ensureBBWasmInit(): Promise<void> {
  if (!bbwasmInitialized) {
    await initBBWasm(`${_package.webRoot}/dist/bbwasm_bg.wasm`);
    bbwasmInitialized = true;
  }
}

/** Run BitBIRCH clustering on a molecule column. Returns an int column with cluster assignments. */
export async function bitbirchClustering(
  molCol: DG.Column,
  threshold: number = 0.65,
  fingerprintType: Fingerprint = Fingerprint.Morgan,
  nFeatures: number = defaultMorganFpLength,
  branchingFactor: number = 50,
): Promise<DG.Column> {
  await ensureBBWasmInit();

  const nBytesPerFp = Math.ceil(nFeatures / 8);
  const rowCount = molCol.length;
  const rdkit = getRdKitModule();

  // Compute fingerprints and pack into contiguous buffer
  const packedFps = new Uint8Array(rowCount * nBytesPerFp);
  let validCount = 0;
  // Map from position in valid-only array back to original row index
  const validIndices: number[] = [];

  for (let i = 0; i < rowCount; i++) {
    if (molCol.isNone(i)) continue;
    const molString = molCol.get(i);
    if (!molString) continue;

    const molSafe = getMolSafe(molString, {}, rdkit);
    if (!molSafe.mol) continue;
    try {
      const fp = getRDKitFpAsUint8Array(molSafe.mol, fingerprintType);
      packedFps.set(fp, validCount * nBytesPerFp);
      validIndices.push(i);
      validCount++;
    } finally {
      molSafe.mol.delete();
    }
  }

  if (validCount === 0) {
    const col = DG.Column.int('Cluster (BitBIRCH)', rowCount);
    col.init((_i) => DG.INT_NULL);
    return col;
  }

  // Run BitBIRCH on valid fingerprints only
  const validFps = packedFps.slice(0, validCount * nBytesPerFp);
  const bb = new WasmBitBirch(threshold, branchingFactor, nFeatures);
  try {
    bb.fit(validFps, validCount);
    const assignments: Int32Array = bb.get_assignments();

    // Map back to full column, leaving nulls for invalid rows
    const result = DG.Column.int('Cluster (BitBIRCH)', rowCount);
    result.init((_i) => DG.INT_NULL);
    for (let j = 0; j < validCount; j++)
      result.set(validIndices[j], assignments[j]);

    return result;
  } finally {
    bb.free();
  }
}
