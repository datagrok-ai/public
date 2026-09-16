/// The process layer (build-plan.md WO-5): tickets from the backlog snapshot, the releases recorded
/// under core/docs/release with the commits they picked, and the people all three name, resolved
/// through the person homes and the autofix roster.
import * as fs from 'fs';
import * as path from 'path';
import {spawnSync} from 'child_process';
import {Emitter} from '../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../context';
import {splitFrontmatter, parseYamlDocument} from '../../frontmatter';
import {HomeSet} from '../../homes';
import {ticketId, relId, custId, commitId} from '../../ids';
import {homesOf, idTokens, ticketTokens, ticketStub, resolveMention} from './markers';

/** The snapshot repo beside the monorepo, then the place it is cloned to on the dev boxes (build-plan.md WO-5). */
const RELEASE_DIR = 'core/docs/release';
const OWNERS = 'autofix/cfg/owners.json';
/** A fix version that names a release rather than a bucket such as `v1` or `Next patch version`. */
const VERSION = /^\d+(\.\d+)*$/;
const GENERATED = /generated\s+(\d{4}-\d{2}-\d{2})/i;
const IDENTITY_KEYS = ['github', 'bitbucket', 'jira'];
/** Separates the git --format fields, so a subject or a message body cannot be mistaken for one. */
const UNIT = '\u001f';
/** A reconstructed record is evidence of the reconstruction, not of the release; every edge it produces says so. */
const DRY_RUN_CONFIDENCE = 0.7;
/** What turns a commit into a claim of resolution, in its own message (conventions.md §6, edges/work/resolves.yaml). */
const RESOLUTION_CLAIM = /\b(?:fix|fixes|fixed|resolve|resolves|resolved|close|closes|closed)\b[\s:#-]*(GROK-\d+|#\d+)/gi;

export const processExtractor: Extractor = {
  name: 'process',
  describes: {backlog: 'tickets, their state and who they are assigned to', git: 'commits and what a release includes',
    people: 'people, teams and customers', process: 'releases, scenarios and tutorials', releases: 'release records'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    new ProcessLayer(ctx, emitter).run();
  },
};

/** A picked commit, with what its own message claims. */
interface Commit {
  id: string;
  sha: string;
  message: string;
}

interface Person {
  name?: string;
  github?: string;
  bitbucket?: string;
  jira?: string;
  slack?: string;
  emails?: string[];
  departed?: boolean;
}

/** The people of this build, loaded once: the packages extractor resolves its authors through them too. */
export function peopleOf(ctx: BuildContext, emitter: Emitter): People {
  return ctx.people ??= new People(ctx, emitter);
}

class ProcessLayer {
  private people: People;
  private homes: HomeSet;
  private gitPartial = false;

  constructor(private ctx: BuildContext, private emitter: Emitter) {
    this.people = peopleOf(ctx, emitter);
    this.homes = homesOf(ctx);
  }

  run(): void {
    this.people.load();
    this.releases();
    this.tickets();
    this.people.finish();
  }

  /** The snapshot the CLI resolved (`backlogRoot`); none means no ticket layer. */
  private backlogDir(): string | undefined {
    const dir = this.ctx.backlogDir;
    return dir && fs.existsSync(path.join(dir, 'index.jsonl')) ? dir : undefined;
  }

  private tickets(): void {
    const dir = this.backlogDir();
    if (!dir) {
      this.emitter.source('backlog', 'missing');
      return;
    }
    let watermark = '';
    for (const line of fs.readFileSync(path.join(dir, 'index.jsonl'), 'utf8').split('\n')) {
      if (!line.trim()) continue;
      let row: any;
      try {
        row = JSON.parse(line);
      }
      catch (e: any) {
        this.emitter.problem('invalid_rows', `backlog/index.jsonl: ${e.message}`);
        continue;
      }
      this.ticket(row, dir);
      if (typeof row.updated === 'string' && row.updated > watermark) watermark = row.updated;
    }
    this.emitter.source('backlog', watermark ? `ok@${watermark}` : 'ok');
  }

  private ticket(r: any, dir: string): void {
    const key = String(r.id);
    const id = ticketId(key);
    const github = r.source === 'github';
    const folder = github ? 'github' : 'jira';
    const file = `backlog/${folder}/${key}.md`;
    const full = path.join(dir, folder, `${key}.md`);
    const opened = fs.existsSync(full) ? splitFrontmatter(fs.readFileSync(full, 'utf8')) : undefined;
    const tracker = opened?.data?.tracker as Record<string, unknown> | undefined;
    const admitted = this.emitter.node({
      type: 'ticket', id, name: String(r.title ?? id), tracker: github ? 'github' : 'jira', key: github ? `#${key.replace(/^public-/, '')}` : key,
      kind: r.kind ?? 'unknown', state: r.status, raw_status: tracker?.raw_status, priority: r.priority ?? undefined,
      components: some(r.components), labels: some(r.labels), resolution: r.resolution ?? undefined, closed: r.closed ?? undefined,
      area: r.area ?? undefined, url: r.url, created: r.created, updated: r.updated,
      assignee: this.people.id(r.assignee), reporter: this.people.id(r.reporter),
      provenance: 'external', source_layer: 'process',
    });
    if (!admitted.accepted) return;
    for (const version of some(r.fix_versions) ?? []) {
      if (!VERSION.test(String(version))) continue;
      this.releaseStub(String(version), 'external');
      this.emitter.edge({type: 'targets-release', from: id, to: relId(String(version)), kind: 'fix-version', derived_by: 'external', confidence: 1, evidence: ['backlog/index.jsonl']});
    }
    for (const customer of some(r.customer) ?? []) {
      const name = this.people.customer(String(customer));
      this.emitter.stub(custId(slug(name)), 'customer', name, 'external', {relationship: 'customer', jira_client: name});
      this.emitter.edge({type: 'requested-by', from: id, to: custId(slug(name)), derived_by: 'external', confidence: 1, evidence: ['backlog/index.jsonl']});
    }
    if (opened) this.body(id, opened.body, file);
  }

  /** The Feature field the snapshot has no column for: `~id` tokens in the ticket text, and the tickets it names. */
  private body(id: string, body: string, file: string): void {
    for (const token of idTokens(body).keys()) {
      const target = resolveMention(this.emitter, this.homes, token, file);
      if (!target) continue;
      this.emitter.edge({type: target.root === 'feature' ? 'affects' : 'mentions', from: id, to: target.id, derived_by: 'annotation', confidence: 1, evidence: [file]});
    }
    for (const [ticket, count] of ticketTokens(body)) {
      if (ticket === id) continue;
      ticketStub(this.emitter, ticket);
      this.emitter.edge({type: 'mentions', from: id, to: ticket, count, derived_by: 'annotation', confidence: 1, evidence: [file]});
    }
  }

  private releases(): void {
    const dir = path.join(this.ctx.repoRoot, ...RELEASE_DIR.split('/'));
    if (!fs.existsSync(dir)) {
      this.emitter.source('releases', 'missing');
      return;
    }
    for (const name of fs.readdirSync(dir).sort())
      if (/\.ya?ml$/.test(name)) this.release(path.join(dir, name), `${RELEASE_DIR}/${name}`);
    this.emitter.source('releases', 'ok');
    this.emitter.source('git', this.gitPartial ? 'partial' : 'ok');
  }

  private release(full: string, record: string): void {
    const text = fs.readFileSync(full, 'utf8');
    const data = parseYamlDocument(text).data;
    if (!data || data.version === undefined) {
      this.emitter.problem('invalid_rows', `${record}: no version`);
      return;
    }
    const version = String(data.version);
    const id = relId(version);
    const base = typeof data.base === 'string' ? data.base.replace(/^release\//, '') : undefined;
    const row: Row = {type: 'release', id, name: version, version, kind: data.kind ?? kindOf(version), state: data.state ?? 'development',
      branch: `release/${version}`, base: base ? relId(base) : undefined, record, released: data.released ?? undefined,
      provenance: 'annotation', source_layer: 'process'};
    if (data.dry_run === true) {
      row.status = 'proposed';
      row.dry_run = true;
      row.description = `Dry run record generated ${GENERATED.exec(text)?.[1] ?? String((data.checks as any)?.checked_at ?? 'without a date')}; not the live record.`;
    }
    if (!this.emitter.node(row).accepted) return;
    const confidence = data.dry_run === true ? DRY_RUN_CONFIDENCE : 1;
    if (base) this.releaseStub(base, 'annotation');
    for (const pick of Array.isArray(data.picks) ? data.picks : []) this.pick(id, pick as Record<string, unknown>, record, confidence);
    for (const group of Object.values(data.features ?? {} as Record<string, unknown>))
      for (const item of Array.isArray(group) ? group : []) {
        const ticket = (item as Record<string, unknown>)?.t;
        if (typeof ticket !== 'string') continue;
        ticketStub(this.emitter, ticketId(ticket));
        this.emitter.edge({type: 'targets-release', from: ticketId(ticket), to: id, kind: 'fix-version', derived_by: 'annotation', confidence, evidence: [record]});
      }
  }

  /**
   * One line of `picks:`: the commit it names, what the release includes, and how the record relates it to a ticket.
   * The association itself is the record's word, not git's, so it is an annotation; resolution is claimed only by the
   * commit message or by an explicit `claim: resolves` on the pick, and anything else is a mention (review 3 #5).
   */
  private pick(release: string, pick: Record<string, unknown>, record: string, confidence: number): void {
    const repo = pick.repo === 'public' ? 'public' : 'reddata';
    const commit = this.commit(repo, String(pick.commit ?? ''), record);
    if (!commit) return;
    this.emitter.edge({type: 'includes', from: release, to: commit.id, derived_by: 'annotation', confidence, evidence: [record]});
    if (typeof pick.ticket === 'string' && pick.ticket) {
      const ticket = ticketId(pick.ticket);
      ticketStub(this.emitter, ticket);
      const claimed = pick.claim === 'resolves' || claimsResolution(commit.message, ticket);
      const evidence = claimed && pick.claim !== 'resolves' ? [record, `${repo}:${commit.sha}`] : [record];
      this.emitter.edge({type: claimed ? 'resolves' : 'mentions', from: commit.id, to: ticket, derived_by: 'annotation', confidence, evidence});
      this.emitter.edge({type: 'targets-release', from: ticket, to: release, kind: 'picked', derived_by: 'annotation', confidence, evidence: [record]});
      return;
    }
    for (const [ticket, count] of ticketTokens(String(pick.subject ?? ''))) {
      ticketStub(this.emitter, ticket);
      this.emitter.edge({type: 'mentions', from: commit.id, to: ticket, count, derived_by: 'annotation', confidence, evidence: [record]});
    }
  }

  /** The full sha, subject, date and message of an abbreviated pick; a sha this checkout does not have leaves git partial. */
  private commit(repo: 'reddata' | 'public', short: string, record: string): Commit | undefined {
    const cwd = repo === 'public' ? path.join(this.ctx.repoRoot, 'public') : this.ctx.repoRoot;
    const shown = spawnSync('git', ['-C', cwd, 'show', '-s', `--format=%H${UNIT}%s${UNIT}%cI${UNIT}%B`, `${short}^{commit}`], {encoding: 'utf8'});
    const [sha, subject, date, message] = shown.status === 0 ? shown.stdout.split(UNIT) : [];
    if (!sha) {
      this.gitPartial = true;
      this.emitter.problem('unresolved_ids', `${record}: ${repo} has no commit ${short}`);
      return undefined;
    }
    const id = commitId(repo, sha);
    if (!this.emitter.node({type: 'commit', id, name: sha.slice(0, 10), repo, sha, subject, date, provenance: 'git', source_layer: 'process'}).accepted) return undefined;
    return {id, sha, message: message ?? ''};
  }

  private releaseStub(version: string, provenance: string): void {
    this.emitter.stub(relId(version), 'release', version, provenance, {version, kind: kindOf(version), state: 'development', branch: `release/${version}`});
  }
}

/**
 * A tracker handle or display name as a person node: a home that declares the handle, a home whose cross-system
 * identifiers match the autofix roster's, then the roster itself as a stub. What none of the three know is listed
 * in `reports/unresolved-people.json` and leaves the people source partial.
 */
export class People {
  private roster: Record<string, Person> = {};
  private customers: Record<string, string> = {};
  private byHandle = new Map<string, string>();
  private byIdentity = new Map<string, string>();
  private byName = new Map<string, string>();
  private byEmail = new Map<string, string>();
  private resolved = new Map<string, string | null>();
  private unresolved = new Map<string, number>();
  private rosterFile = true;
  private loaded = false;

  constructor(private ctx: BuildContext, private emitter: Emitter) {}

  load(): void {
    if (this.loaded) return;
    this.loaded = true;
    const owners = path.join(this.ctx.repoRoot, ...OWNERS.split('/'));
    this.rosterFile = fs.existsSync(owners);
    if (this.rosterFile) this.roster = JSON.parse(fs.readFileSync(owners, 'utf8')).people ?? {};
    for (const [key, person] of Object.entries(this.roster))
      if (person.name) this.byName.set(person.name.toLowerCase(), key);
    for (const home of homesOf(this.ctx).homes) {
      if (home.type.root !== 'actor') continue;
      const data = home.data;
      if (typeof data.handle === 'string') this.byHandle.set(data.handle, home.id);
      for (const key of IDENTITY_KEYS)
        if (typeof data[key] === 'string') this.byIdentity.set(`${key}:${data[key]}`, home.id);
      for (const email of Array.isArray(data.emails) ? data.emails : []) this.byIdentity.set(`email:${email}`, home.id);
    }
    this.customers = this.backlogTaxonomy();
  }

  /** The canonical Jira Client name of a backlog customer label (`taxonomy.yaml` customers map). */
  customer(label: string): string {
    return this.customers[label] ?? label;
  }

  /**
   * The person a `package.json` author names: the home or the roster entry that lists the email, else the local part
   * as the handle. Emitted as a stub either way, so the owner edge stands wherever the person record lands later.
   */
  byAuthorEmail(email: string, name?: string): string {
    this.load();
    const key = email.toLowerCase();
    const cached = this.byEmail.get(key);
    if (cached) return cached;
    const handle = key.split('@')[0];
    const known = this.byIdentity.get(`email:${key}`) ?? this.fromRoster(key);
    const id = known ?? `P:${handle}`;
    this.emitter.stub(id, 'person', name || handle, 'registry', known ? {emails: [email]} : {handle, emails: [email]});
    this.byEmail.set(key, id);
    return id;
  }

  id(key: unknown): string | undefined {
    if (typeof key !== 'string' || !key) return undefined;
    const cached = this.resolved.get(key);
    if (cached !== undefined) {
      if (cached === null) this.unresolved.set(key, (this.unresolved.get(key) ?? 0) + 1);
      return cached ?? undefined;
    }
    const id = this.lookup(key);
    this.resolved.set(key, id ?? null);
    if (!id) this.unresolved.set(key, (this.unresolved.get(key) ?? 0) + 1);
    return id;
  }

  finish(): void {
    this.emitter.source('people', !this.rosterFile ? 'missing' : this.unresolved.size ? 'partial' : 'ok');
    if (this.unresolved.size)
      this.emitter.report('unresolved-people', [...this.unresolved].sort(([a, x], [b, y]) => y - x || (a < b ? -1 : 1)).map(([key, count]) => ({key, count})));
  }

  private lookup(key: string): string | undefined {
    const home = this.byHandle.get(key);
    if (home) return home;
    const rosterKey = this.roster[key] ? key : this.byName.get(key.toLowerCase());
    if (!rosterKey) return undefined;
    const known = this.byHandle.get(rosterKey);
    if (known) return known;
    const person = this.roster[rosterKey];
    for (const identity of identities(person)) {
      const match = this.byIdentity.get(identity);
      if (match) return match;
    }
    const id = `P:${rosterKey}`;
    this.emitter.stub(id, 'person', person.name ?? rosterKey, 'external', {handle: rosterKey, github: person.github, bitbucket: person.bitbucket,
      jira: person.jira, slack: person.slack, emails: person.emails, departed: person.departed});
    return id;
  }

  /** The roster person who lists [email], resolved the way a tracker handle is. */
  private fromRoster(email: string): string | undefined {
    for (const [key, person] of Object.entries(this.roster))
      if ((person.emails ?? []).some((e) => e.toLowerCase() === email)) return this.lookup(key);
    return undefined;
  }

  private backlogTaxonomy(): Record<string, string> {
    const file = this.ctx.backlogDir && path.join(this.ctx.backlogDir, 'taxonomy.yaml');
    if (file && fs.existsSync(file)) {
      const customers = parseYamlDocument(fs.readFileSync(file, 'utf8')).data?.customers;
      return customers && typeof customers === 'object' ? customers as Record<string, string> : {};
    }
    return {};
  }
}

/** Whether [message] itself claims to resolve [ticket]: `Fixes GROK-1`, `closes #7` and their variants. */
function claimsResolution(message: string, ticket: string): boolean {
  for (const m of message.matchAll(RESOLUTION_CLAIM))
    if (ticketId(m[1]) === ticket) return true;
  return false;
}

function identities(person: Person): string[] {
  return [...IDENTITY_KEYS.filter((k) => (person as any)[k]).map((k) => `${k}:${(person as any)[k]}`),
    ...(person.emails ?? []).map((e) => `email:${e}`)];
}

/** A patch release has a non-zero third segment; everything else is a minor. */
function kindOf(version: string): string {
  return Number(version.split('.')[2] ?? 0) > 0 ? 'patch' : 'minor';
}

/** `Insitro` -> `insitro`, `pharm-sphere` -> `pharm-sphere`, `JnJ` -> `jnj`. */
function slug(name: string): string {
  return name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
}

/** A list the snapshot leaves empty is left out of the row rather than written as `[]`. */
function some(value: unknown): unknown[] | undefined {
  return Array.isArray(value) && value.length ? value : undefined;
}
