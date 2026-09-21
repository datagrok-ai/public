/// Media (conventions.md §5.7, notation.md §3.7): every media file under the media roots and every file a page
/// embeds as a `media` node with its blob, bytes, delivery url, thumbnail and record fields; every hosted video a
/// page embeds or a record describes; one `embeds` edge per occurrence the docs extractor collected; `illustrates`
/// from the records, at annotation 1.0 when reviewed and as an llm proposal at 0.6 otherwise.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../emitter';
import {Row, compare} from '../../normalize';
import {BuildContext, Extractor} from '../context';
import {HOME_IGNORE, HomeSet} from '../../homes';
import {isMedia, formatOf, blobIds, hashBlob, MediaRecord} from '../../media';
import {Embed, targetKey} from '../../embeds';
import {mediaId, hostedMediaId, docId, sourceLayerOf} from '../../ids';
import {Roots, LANDING_PREFIX, localPath, existsAt, deliveryUrl} from '../../roots';
import {homesOf} from './markers';

/** Where every media file is a node, embedded or not; a file elsewhere is one only when a page embeds or a record describes it. */
export const MEDIA_ROOTS = ['public/help', 'public/docusaurus/static', 'core/docs'];
/** The same, in the site's checkout. */
const LANDING_ROOT = 'web';
const PROPOSAL_CONFIDENCE = 0.6;

export const mediaExtractor: Extractor = {
  name: 'media',
  describes: {media: 'media files, hosted videos, the pages that show them and the records that describe them'},
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    new MediaLayer(ctx, emitter, homesOf(ctx)).run();
  },
};

class MediaLayer {
  private broken = 0;
  private roots: Roots;

  constructor(private ctx: BuildContext, private emitter: Emitter, private homes: HomeSet) {
    this.roots = {repoRoot: ctx.repoRoot, landingDir: ctx.landingDir};
  }

  run(): void {
    const records = this.homes.media;
    const files = new Set<string>();
    for (const root of MEDIA_ROOTS)
      for (const f of globSync(`${root}/**/*`, {cwd: this.ctx.repoRoot, ignore: HOME_IGNORE, nodir: true, posix: true}))
        if (isMedia(f)) files.add(f);
    if (this.ctx.landingDir !== undefined)
      for (const f of globSync(`${LANDING_ROOT}/**/*`, {cwd: this.ctx.landingDir, ignore: HOME_IGNORE, nodir: true, posix: true}))
        if (isMedia(f)) files.add(`${LANDING_PREFIX}${f}`);
    for (const r of records.values())
      if (r.path) files.add(r.path);
    const shown: {page: string, embed: Embed}[] = [];
    const posters = new Map<string, string>();
    for (const {page, embed} of this.emitter.embeds) {
      if (embed.target.kind === 'file') {
        // a site path with no site given is neither shown nor broken: the source is missing
        if (embed.target.path.startsWith(LANDING_PREFIX) && this.ctx.landingDir === undefined) continue;
        if (!this.exists(embed.target.path)) {
          this.broken++;
          this.emitter.problem('broken_embeds', `${page}:${embed.line}: ${embed.raw}`);
          continue;
        }
        files.add(embed.target.path);
      }
      if (embed.poster && this.exists(embed.poster)) {
        files.add(embed.poster);
        if (embed.target.kind === 'hosted' && !posters.has(targetKey(embed.target))) posters.set(targetKey(embed.target), embed.poster);
      }
      shown.push({page, embed});
    }
    const blobs = blobIds(this.roots);
    const embedded = new Set(shown.filter((s) => s.embed.target.kind === 'file').map((s) => targetKey(s.embed.target)));
    const admitted = new Set<string>();
    for (const file of [...files].sort(compare)) {
      const id = mediaId(file);
      const record = records.get(id);
      let blob = blobs.get(file);
      // a git tree lists what it tracks; a local addition no page shows is reported, not indexed (a tree with no git tracks nothing)
      if (blob === undefined && blobs.size) {
        this.emitter.problem('untracked_media', file);
        if (!embedded.has(file) && !record) continue;
      }
      const local = localPath(this.roots, file)!;
      blob ??= hashBlob(fs.readFileSync(local));
      const thumb = /\.gif$/i.test(file) ? `${file.slice(0, -4)}-thumb.png` : undefined;
      const row: Row = {type: 'media', id, name: path.posix.basename(file), path: file, format: formatOf(file), blob,
        bytes: fs.statSync(local).size, url: deliveryUrl(file), thumbnail: thumb && files.has(thumb) ? mediaId(thumb) : undefined,
        provenance: record ? 'annotation' : 'filesystem', source_layer: sourceLayerOf(file), ...record?.data};
      if (!this.emitter.node(row).accepted) continue;
      admitted.add(id);
      if (record) this.illustrates(record);
    }
    const hosted = new Map<string, {provider: string, externalId: string}>();
    for (const {embed} of shown)
      if (embed.target.kind === 'hosted') hosted.set(hostedMediaId(embed.target.provider, embed.target.id), {provider: embed.target.provider, externalId: embed.target.id});
    for (const r of records.values())
      if (r.provider && r.externalId) hosted.set(r.id, {provider: r.provider, externalId: r.externalId});
    for (const [id, {provider, externalId}] of [...hosted].sort(([a], [b]) => compare(a, b))) {
      const record = records.get(id);
      const {name, ...data} = record?.data ?? {};
      const poster = posters.get(`${provider}:${externalId}`);
      const row: Row = {type: 'media', id, name: name ?? `${provider} ${externalId}`, url: `https://www.youtube.com/watch?v=${externalId}`, format: 'youtube',
        provider, external_id: externalId, thumbnail: poster && files.has(poster) ? mediaId(poster) : undefined,
        provenance: record ? 'annotation' : 'ast', source_layer: 'public', ...data};
      if (!this.emitter.node(row).accepted) continue;
      admitted.add(id);
      if (record) this.illustrates(record);
    }
    for (const {page, embed} of shown) {
      const to = embed.target.kind === 'file' ? mediaId(embed.target.path) : hostedMediaId(embed.target.provider, embed.target.id);
      if (!admitted.has(to)) continue;
      this.emitter.edge({type: 'embeds', from: docId(page), to, derived_by: 'ast', confidence: 1, evidence: [page], position: embed.position, line: embed.line,
        form: embed.form, anchor: embed.anchor, alt: embed.alt, title: embed.title, caption: embed.caption, start_seconds: embed.start_seconds});
    }
    this.emitter.source('media', this.broken ? 'partial' : 'ok');
  }

  /** What the record asserts: a reviewed record is authored fact, an unreviewed one a proposal below the review line (conventions.md §8). */
  private illustrates(record: MediaRecord): void {
    const reviewed = record.data.reviewed === true;
    for (const target of record.illustrates) {
      const home = this.homes.index.byId.get(target) ?? this.homes.index.byAlias.get(target);
      this.emitter.edge({type: 'illustrates', from: record.id, to: home?.id ?? target, derived_by: reviewed ? 'annotation' : 'llm',
        confidence: reviewed ? 1 : PROPOSAL_CONFIDENCE, evidence: [record.file]});
    }
  }

  private exists(file: string): boolean {
    return existsAt(this.roots, file);
  }
}
