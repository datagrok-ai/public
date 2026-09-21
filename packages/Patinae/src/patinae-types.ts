/** Hand-written subset of the @patinae/viewer API (web/ts/src/core/api.ts, types.ts upstream).
 * The vendored bundle in wasm/ is imported at runtime by URL, so no .d.ts is available. */

export interface CommandMessage {
  level: 'info' | 'warning' | 'error' | 'clear';
  text: string;
}

export interface CommandOutput {
  messages: CommandMessage[];
}

export interface ObjectInfo {
  name: string;
  object_type: 'molecule' | 'map' | 'measurement' | 'label';
  atom_count: number;
  enabled: boolean;
}

export type PanelName = 'repl' | 'objects' | 'sequence' | 'movie';
export type PanelSlot = 'top' | 'right' | 'bottom';

export interface PatinaeViewerOptions {
  layout?: {name: PanelName, slot: PanelSlot, collapsed?: boolean}[];
  slots?: Partial<Record<PanelSlot, HTMLElement>>;
  picking?: boolean;
  selectionOverlay?: boolean;
}

export interface PatinaeViewer {
  init(): Promise<void>;
  loadData(data: Uint8Array, name: string, format: string): void;
  execute(command: string): CommandOutput;
  /** `execute` for `fetch`/`load <url>` lines is fire-and-forget; this variant awaits them. */
  executeAsync(command: string): Promise<CommandOutput>;
  getObjectInfos(): ObjectInfo[];
  destroy(): void;
}

export interface PatinaeModule {
  PatinaeViewer: new (container: HTMLElement, options?: PatinaeViewerOptions) => PatinaeViewer;
}
