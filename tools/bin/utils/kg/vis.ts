/// The render tier of the graph browser (core/docs/knowledge-graph/vis/): one binary blob per generation with
/// every node and edge as typed arrays, the ids and names beside it, and the type tree the panels are built
/// from. Derived from the JSONL after the manifest, cached under `<gen>/vis/`, never part of the batch id.
import * as fs from 'fs';
import * as path from 'path';
import {TypeSystem, NodeType, EdgeType, nodeOrder, edgeOrder, graphLabel, firstSentence} from './types';
import {readJsonl} from './kuzu';

export const VIS_DIR = 'vis';
export const BLOB = 'graph.bin';
export const INDEX = 'index.json';
export const SCHEMA = 'schema.json';
const MAGIC = 'KGV1';
const ENUMS = ['layer', 'visibility', 'status', 'provenance', 'derivedBy'] as const;

export interface Section {
  dtype: 'u8' | 'u32' | 'i32';
  offset: number;
  length: number;
}

export interface BlobHeader {
  version: 1;
  batch: string;
  nodes: number;
  edges: number;
  /** Edges whose end is not a node of the generation, left out of the blob. */
  dropped: number;
  nodeTypes: string[];
  /** Edge type names and, for `ref` rows, the reference property name: the Kuzu rel table each maps to is in schema.json. */
  edgeKinds: string[];
  enums: Record<typeof ENUMS[number], string[]>;
  sections: Record<string, Section>;
}

export interface VisSummary {
  dir: string;
  nodes: number;
  edges: number;
  dropped: number;
  bytes: number;
}

export function visDir(genDir: string): string {
  return path.join(genDir, VIS_DIR);
}

/** Bumped when the exported files change shape; an older cache is exported again on the next serve. */
export const VIS_VERSION = 2;

export function hasVis(genDir: string): boolean {
  if (![BLOB, INDEX, SCHEMA].every((f) => fs.existsSync(path.join(visDir(genDir), f)))) return false;
  try {
    return JSON.parse(fs.readFileSync(path.join(visDir(genDir), SCHEMA), 'utf8')).version === VIS_VERSION;
  }
  catch {
    return false;
  }
}

/** Writes graph.bin, index.json and schema.json under `<genDir>/vis/` from the generation's JSONL. */
export function exportVis(genDir: string, system: TypeSystem, batch: string): VisSummary {
  const dir = visDir(genDir);
  fs.mkdirSync(dir, {recursive: true});
  const enums = new EnumTables();
  const ids: string[] = [];
  const names: string[] = [];
  const index = new Map<string, number>();
  const nodeTypes: string[] = [];
  const type: number[] = [], layer: number[] = [], visibility: number[] = [], status: number[] = [], provenance: number[] = [];
  for (const t of nodeOrder(system)) {
    const file = path.join(genDir, 'data', 'nodes', `${t.name}.jsonl`);
    if (t.abstract || !fs.existsSync(file)) continue;
    const typeIndex = nodeTypes.push(t.name) - 1;
    for (const row of readJsonl(file)) {
      const id = String(row.id);
      index.set(id, ids.length);
      ids.push(id);
      names.push(String(row.name ?? id.split('/').pop()));
      type.push(typeIndex);
      layer.push(enums.code('layer', row.source_layer));
      visibility.push(enums.code('visibility', row.visibility));
      status.push(enums.code('status', row.status));
      provenance.push(enums.code('provenance', row.provenance));
    }
  }
  const edgeKinds: string[] = [];
  const from: number[] = [], to: number[] = [], kind: number[] = [], confidence: number[] = [], derivedBy: number[] = [];
  const degree = new Uint32Array(ids.length);
  const parent = new Int32Array(ids.length).fill(-1);
  let dropped = 0;
  for (const file of fs.readdirSync(path.join(genDir, 'data', 'edges')).filter((f) => f.endsWith('.jsonl')).sort()) {
    const name = file.slice(0, -'.jsonl'.length);
    const kindIndex = edgeKinds.push(name) - 1;
    for (const row of readJsonl(path.join(genDir, 'data', 'edges', file))) {
      const a = index.get(String(row.from)), b = index.get(String(row.to));
      if (a === undefined || b === undefined) {
        dropped++;
        continue;
      }
      from.push(a);
      to.push(b);
      kind.push(kindIndex);
      confidence.push(Math.round(Math.min(1, Math.max(0, Number(row.confidence ?? 1))) * 100));
      derivedBy.push(enums.code('derivedBy', row.derived_by));
      degree[a]++;
      degree[b]++;
      if (name === 'part-of') parent[a] = b;
    }
  }
  const sections = new SectionWriter();
  sections.add('type', 'u8', type);
  sections.add('layer', 'u8', layer);
  sections.add('visibility', 'u8', visibility);
  sections.add('status', 'u8', status);
  sections.add('provenance', 'u8', provenance);
  sections.add('degree', 'u32', degree);
  sections.add('parent', 'i32', parent);
  sections.add('from', 'u32', from);
  sections.add('to', 'u32', to);
  sections.add('kind', 'u8', kind);
  sections.add('confidence', 'u8', confidence);
  sections.add('derivedBy', 'u8', derivedBy);
  const header: BlobHeader = {version: 1, batch, nodes: ids.length, edges: from.length, dropped, nodeTypes, edgeKinds,
    enums: enums.tables(), sections: sections.table()};
  const blob = sections.blob(header);
  fs.writeFileSync(path.join(dir, BLOB), blob);
  fs.writeFileSync(path.join(dir, INDEX), JSON.stringify({ids, names}));
  fs.writeFileSync(path.join(dir, SCHEMA), JSON.stringify(schemaJson(system), null, 1));
  return {dir, nodes: ids.length, edges: from.length, dropped, bytes: blob.length};
}

/** The header, then every section at the offset the header names; readers view the sections without copying. */
export function readBlobHeader(file: string): BlobHeader {
  const fd = fs.openSync(file, 'r');
  try {
    const head = Buffer.alloc(8);
    fs.readSync(fd, head, 0, 8, 0);
    if (head.toString('latin1', 0, 4) !== MAGIC) throw new Error(`${file}: not a graph blob`);
    const length = head.readUInt32LE(4);
    const json = Buffer.alloc(length);
    fs.readSync(fd, json, 0, length, 8);
    return JSON.parse(json.toString('utf8'));
  }
  finally {
    fs.closeSync(fd);
  }
}

export interface SchemaJson {
  version: number;
  roots: string[];
  /** Edge groups in folder order (conventions.md §4), with what each holds. */
  edgeGroups: string[];
  nodeTypes: {name: string, extends?: string, root?: string, abstract: boolean, prefix?: string, hierarchical: boolean, authored: boolean,
    description: string, members: string[]}[];
  edgeTypes: {name: string, label: string, group: string, extends?: string, abstract: boolean, from: string[], to: string[], description: string,
    properties: string[]}[];
  /** Reference properties, each a rel table of its own: which node types declare it and what it points at. */
  refs: {name: string, declaredBy: string[], to: string[]}[];
}

export function schemaJson(system: TypeSystem): SchemaJson {
  const refs = new Map<string, {name: string, declaredBy: string[], to: string[]}>();
  for (const t of system.nodes.values())
    for (const m of Object.values(t.own)) {
      if (m.kind !== 'ref') continue;
      let ref = refs.get(m.name);
      if (!ref) refs.set(m.name, ref = {name: m.name, declaredBy: [], to: []});
      ref.declaredBy.push(t.name);
      for (const r of m.refs ?? [])
        if (!ref.to.includes(r)) ref.to.push(r);
    }
  const edgeGroups: string[] = [];
  for (const e of edgeOrder(system))
    if (e.group && !edgeGroups.includes(e.group)) edgeGroups.push(e.group);
  return {
    version: VIS_VERSION,
    roots: system.roots,
    edgeGroups,
    nodeTypes: nodeOrder(system).map((t: NodeType) => ({name: t.name, extends: t.extends, root: t.root, abstract: t.abstract, prefix: t.prefix,
      hierarchical: t.hierarchical, authored: t.authored, description: firstSentence(t.description), members: Object.keys(t.members)})),
    edgeTypes: edgeOrder(system).map((e: EdgeType) => ({name: e.name, label: graphLabel(e.name), group: e.group, extends: e.extends, abstract: e.abstract,
      from: e.from, to: e.to, description: firstSentence(e.description), properties: Object.keys(e.properties)})),
    refs: [...refs.values()].sort((a, b) => a.name < b.name ? -1 : 1),
  };
}

class EnumTables {
  private values = new Map<string, string[]>(ENUMS.map((e) => [e, []]));

  code(table: typeof ENUMS[number], value: unknown): number {
    const list = this.values.get(table)!;
    const text = value === undefined || value === null ? '' : String(value);
    let i = list.indexOf(text);
    if (i < 0) i = list.push(text) - 1;
    if (i > 255) throw new Error(`enum ${table} has more than 256 values`);
    return i;
  }

  tables(): Record<typeof ENUMS[number], string[]> {
    return Object.fromEntries(this.values) as Record<typeof ENUMS[number], string[]>;
  }
}

/** Sections are 4-byte aligned so a `Uint32Array` view can start at any of them. */
class SectionWriter {
  private parts: {name: string, dtype: Section['dtype'], data: Buffer}[] = [];

  add(name: string, dtype: Section['dtype'], values: ArrayLike<number>): void {
    const array = dtype === 'u8' ? Uint8Array.from(values) : dtype === 'u32' ? Uint32Array.from(values) : Int32Array.from(values);
    this.parts.push({name, dtype, data: Buffer.from(array.buffer, array.byteOffset, array.byteLength)});
  }

  table(): Record<string, Section> {
    const sections: Record<string, Section> = {};
    let offset = 0;
    for (const p of this.parts) {
      sections[p.name] = {dtype: p.dtype, offset, length: p.dtype === 'u8' ? p.data.length : p.data.length / 4};
      offset += align(p.data.length);
    }
    return sections;
  }

  blob(header: BlobHeader): Buffer {
    const json = Buffer.from(JSON.stringify(header), 'utf8');
    const head = Buffer.alloc(align(8 + json.length));
    head.write(MAGIC, 0, 'latin1');
    head.writeUInt32LE(json.length, 4);
    json.copy(head, 8);
    const body: Buffer[] = [];
    for (const p of this.parts) {
      body.push(p.data);
      const pad = align(p.data.length) - p.data.length;
      if (pad) body.push(Buffer.alloc(pad));
    }
    return Buffer.concat([head, ...body]);
  }
}

function align(n: number): number {
  return (n + 3) & ~3;
}
