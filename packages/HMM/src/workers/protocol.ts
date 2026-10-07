import type {SearchOptions, SequenceInput} from '../hmmer/engine.ts';
import type {AnarciOptions} from '../hmmer/anarci/anarci.ts';

/** Messages from the pool to a HMMER worker. */
export type WorkerRequest =
  | {op: 'init'; module: WebAssembly.Module; base: string}
  | {op: 'anarci'; sequences: [string, string][]; options: AnarciOptions}
  | {op: 'load'; key: string; bytes: ArrayBuffer}
  | {op: 'unload'; key: string}
  | {op: 'scan'; key: string; sequences: SequenceInput[]; options: SearchOptions}
  | {op: 'search'; key: string; model: number; sequences: SequenceInput[]; options: SearchOptions; offset: number};

export type WorkerResponse =
  | {id: number; ok: true; result: unknown}
  | {id: number; ok: false; error: string};
