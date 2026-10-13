import type {SearchOptions, SequenceInput} from '../hmmer/engine.ts';
import type {AnarciOptions} from '../hmmer/anarci/anarci.ts';
import type {Scheme} from '../hmmer/anarci/types.ts';
import type {DomainChoice} from '../numbering-rows';
import type {LogoMode, QueryModelOptions} from '../hmmer/tools.ts';

type Named = {name: string; residues: string};

/** Messages from the pool to a HMMER worker. */
export type WorkerRequest =
  | {op: 'init'; module: WebAssembly.Module; base: string}
  | {op: 'anarci'; sequences: [string, string][]; options: AnarciOptions}
  // ANARCI with the results reduced to numbering frame rows in the worker (numbering-rows.ts).
  | {op: 'anarciRows'; sequences: [string, string][]; options: AnarciOptions; scheme: Scheme;
      choice: DomainChoice}
  | {op: 'load'; key: string; bytes: ArrayBuffer}
  | {op: 'unload'; key: string}
  | {op: 'scan'; key: string; sequences: SequenceInput[]; options: SearchOptions}
  | {op: 'search'; key: string; model: number; sequences: SequenceInput[]; options: SearchOptions; offset: number}
  // The tools module (hmmer_tools.wasm), compiled by the pool on first use.
  | {op: 'queryModels'; tools: WebAssembly.Module; queries: Named[]; options: QueryModelOptions}
  | {op: 'alignTraces'; tools: WebAssembly.Module; model: Uint8Array; index: number; sequences: Named[]}
  | {op: 'alignRows'; tools: WebAssembly.Module; model: Uint8Array; index: number; sequences: Named[];
      chunks: Uint8Array[]; trim: boolean}
  | {op: 'logos'; tools: WebAssembly.Module; models: Uint8Array; mode: LogoMode}
  | {op: 'extract'; tools: WebAssembly.Module; models: Uint8Array; index: number}
  | {op: 'build'; tools: WebAssembly.Module; rows: {name: string; aligned: string}[]; name: string;
      alphabet: 'amino' | 'dna' | 'rna'};

export type WorkerResponse =
  | {id: number; ok: true; result: unknown}
  | {id: number; ok: false; error: string};
