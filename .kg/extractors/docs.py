"""Extract docs layer:

- DocPage + HelpAnchor + LinksTo from `help/**/*.md`
- DocPage from plugin `README.md`, `META.md`, `CLAUDE.md`
- `Documents` edges from `//help-url:` annotations in `package.g.ts`
- Tutorial + TutorialTrack from `packages/Tutorials/src/tracks/**/*.ts`
  (and any plugin's `tutorials/` directory using the same pattern)

URL normalization (so deduplication works):
  https://datagrok.ai/help/visualize/viewers/scatter-plot.md  ->  /help/visualize/viewers/scatter-plot
  /help/visualize/viewers/scatter-plot.md#fragment            ->  /help/visualize/viewers/scatter-plot + anchor 'fragment'
  /help/visualize/viewers/scatter-plot/                       ->  /help/visualize/viewers/scatter-plot
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import yaml

from schema.base import Provenance, SourceLayer
from schema.entities import DocPage, HelpAnchor, Tutorial, TutorialTrack
from schema.relations import Documents, HasAnchor, HasTutorial, LinksTo

from . import _common as c

EXTRACTOR_NAME = "docs"


# ---------------------------------------------------------------------------
# URL normalization
# ---------------------------------------------------------------------------

def normalize_help_url(raw: str) -> tuple[str, str | None] | None:
    """Return (canonical_url, anchor_or_None) or None if unparseable.

    canonical_url has no scheme/host, no trailing slash, no `.md` suffix.
    """
    if not raw:
        return None
    s = raw.strip()
    # Strip protocol+host
    s = re.sub(r"^https?://(?:www\.)?datagrok\.ai", "", s)
    s = re.sub(r"^https?://[^/]+", "", s)
    # Split anchor
    anchor: str | None = None
    if "#" in s:
        s, _, anchor = s.partition("#")
        anchor = anchor.strip() or None
    # Strip .md and trailing /
    s = re.sub(r"\.md$", "", s)
    s = s.rstrip("/")
    if not s.startswith("/"):
        s = "/" + s
    return s, anchor


def doc_id_for(canonical_url: str) -> str:
    return f"doc:{canonical_url}"


def anchor_id_for(canonical_url: str, anchor: str) -> str:
    return f"anchor:{canonical_url}#{anchor}"


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---\n", re.DOTALL)
LEGACY_TITLE_RE = re.compile(r"<!--\s*TITLE:\s*(.*?)\s*-->", re.IGNORECASE)
HEADING_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$", re.MULTILINE)
MD_LINK_RE = re.compile(r"\[(?:[^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
H1_RE = re.compile(r"^#\s+(.+?)\s*$", re.MULTILINE)


def _slugify(s: str) -> str:
    s = s.lower().strip()
    s = re.sub(r"[`*_]", "", s)
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")


def _parse_frontmatter(text: str) -> tuple[dict, str]:
    """Strip the YAML frontmatter (or legacy HTML title) and return
    (frontmatter_dict, remainder)."""
    fm: dict = {}
    m = FRONTMATTER_RE.match(text)
    if m:
        try:
            fm = yaml.safe_load(m.group(1)) or {}
        except yaml.YAMLError:
            fm = {}
        text = text[m.end():]
    elif (m := LEGACY_TITLE_RE.search(text)):
        fm["title"] = m.group(1).strip()
    return fm, text


def _classify_help_path(rel_path: str) -> tuple[bool, str | None, str | None]:
    """(is_public, audience, doc_kind) from the help/ relative path."""
    parts = rel_path.split("/")
    is_public = "_internal" not in parts and not any(p.startswith("_") for p in parts[:-1])
    top = parts[0] if parts else ""
    if top == "develop":
        return is_public, "developer", "guide"
    if top == "deploy":
        if "releases" in parts and "plugins" in parts:
            return is_public, "user", "release-note"
        return is_public, "admin", "guide"
    if top == "govern":
        return is_public, "admin", "guide"
    if top in {"explore", "transform", "visualize", "compute", "access",
                "collaborate", "learn", "datagrok", "domains"}:
        return is_public, "user", "guide"
    if top == "_internal":
        return False, "internal", "guide"
    return is_public, "user", "guide"


# ---------------------------------------------------------------------------
# help/ walker
# ---------------------------------------------------------------------------

def _walk_help(repo_root: Path, bundle, url_to_doc: dict[str, DocPage]) -> int:
    help_root = repo_root / "help"
    if not help_root.is_dir():
        return 0
    n_pages = 0
    for md in help_root.rglob("*.md"):
        try:
            text = md.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        rel_to_help = md.relative_to(help_root).as_posix()
        url = "/help/" + re.sub(r"\.md$", "", rel_to_help)
        url = url.rstrip("/")
        fm, body = _parse_frontmatter(text)
        title = fm.get("title")
        if not title:
            mh = H1_RE.search(body)
            title = mh.group(1).strip() if mh else md.stem
        is_public, audience, doc_kind = _classify_help_path(rel_to_help)
        if md.name.startswith("_") or md.parent.name.startswith("_"):
            is_public = False
        page = DocPage(
            id=doc_id_for(url),
            name=title or md.stem,
            source_layer=SourceLayer.DOCS,
            title=title,
            audience=audience,
            doc_kind=fm.get("doc_kind") or doc_kind,
            public=is_public,
            frontmatter=fm if isinstance(fm, dict) else {},
            paths=[c.file_ref(md, role="markdown")],
            extras={"url": url, "filename": md.name},
        )
        bundle.add(page)
        url_to_doc[url] = page
        n_pages += 1
        # Anchors
        for hm in HEADING_RE.finditer(body):
            level = len(hm.group(1))
            head = hm.group(2).strip()
            slug = _slugify(head)
            if not slug or level == 1:
                continue
            ha = HelpAnchor(
                id=anchor_id_for(url, slug),
                name=head, source_layer=SourceLayer.DOCS,
                doc_id=page.id, anchor=slug, level=level,
                paths=[c.file_ref(md, role="anchor")],
            )
            bundle.add(ha)
            bundle.add(HasAnchor(from_id=page.id, to_id=ha.id,
                                 derived_by=Provenance.AST))
        # Outgoing links
        for lm in MD_LINK_RE.finditer(body):
            href = lm.group(1)
            if href.startswith(("#", "mailto:", "http://", "javascript:")):
                continue
            if href.startswith("https://") and "datagrok.ai" not in href:
                continue
            target_anchor: str | None = None
            if "#" in href:
                href_path, _, target_anchor = href.partition("#")
                target_anchor = target_anchor or None
            else:
                href_path = href
            # Resolve relative path from this md file
            if href_path.startswith("/"):
                # absolute under site root; treat /help/... as canonical
                resolved = re.sub(r"\.md$", "", href_path).rstrip("/")
            elif href.startswith("https://"):
                norm = normalize_help_url(href)
                if not norm:
                    continue
                resolved, target_anchor = norm[0], norm[1] or target_anchor
            else:
                # relative path — resolve against this file's url dir
                base_dir = "/help/" + (rel_to_help.rsplit("/", 1)[0] if "/" in rel_to_help else "")
                # Use posix join semantics
                from pathlib import PurePosixPath
                joined = (PurePosixPath(base_dir) / href_path).as_posix()
                # Normalize ../ segments
                resolved = re.sub(r"//+", "/", joined)
                resolved = re.sub(r"\.md$", "", resolved).rstrip("/")
                # Manually collapse segments
                segs: list[str] = []
                for s in resolved.split("/"):
                    if s == "..":
                        if segs and segs[-1] not in ("", ".."):
                            segs.pop()
                    elif s and s != ".":
                        segs.append(s)
                resolved = "/" + "/".join(segs)
            if not resolved:
                continue
            tgt_doc_id = doc_id_for(resolved)
            tgt_id = anchor_id_for(resolved, target_anchor) if target_anchor else tgt_doc_id
            bundle.add(LinksTo(from_id=page.id, to_id=tgt_id,
                               derived_by=Provenance.AST,
                               extras={"href": href}))
    print(f"[{EXTRACTOR_NAME}] walked help/: {n_pages} pages")
    return n_pages


# ---------------------------------------------------------------------------
# package.g.ts: emit Documents edges from //help-url:
# ---------------------------------------------------------------------------

THIS_HELP_URL_RE = re.compile(
    r"this\.helpUrl\s*=\s*['\"](?P<url>[^'\"]+)['\"]"
)
CLASS_DECL_RE = re.compile(
    r"^\s*(?:export\s+)?(?:abstract\s+)?class\s+(?P<cls>[A-Z][\w]*)\s+",
    re.MULTILINE,
)


def _camel_split(name: str) -> list[str]:
    return re.findall(r"[A-Z][a-z0-9]*|[A-Z]+(?=[A-Z]|$)", name)


def _emit_help_url_edges(repo_root: Path, bundle,
                          url_to_doc: dict[str, DocPage],
                          package_filter: list[str] | None) -> int:
    """Emit Documents edges from three sources:

    1. `//help-url:` in `package.g.ts` annotation blocks  -> Function
    2. `//help-url:` in ApiSamples scripts                -> Package (proxy until ApiSample extractor exists)
    3. `this.helpUrl = '...'` assignments in viewer TS    -> Function (by class-name fuzzy match) or Package
    """
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return 0
    n = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        # 1) package.g.ts annotation //help-url:
        g_path = pkg_dir / "src" / "package.g.ts"
        if g_path.is_file():
            text = g_path.read_text(encoding="utf-8", errors="replace")
            for start, end, block_lines, anchor in c.find_annotation_blocks(text):
                ann = c.parse_annotation_block(block_lines)
                help_url = ann.get("help_url")
                name = ann.get("name")
                if not help_url or not name:
                    continue
                norm = normalize_help_url(help_url)
                if not norm:
                    continue
                url, anc = norm
                tgt = anchor_id_for(url, anc) if anc else doc_id_for(url)
                bundle.add(Documents(
                    from_id=tgt,
                    to_id=c.func_id(pkg_dir.name, name),
                    derivation="help_url_annotation",
                    derived_by=Provenance.ANNOTATION,
                    extras={"raw": help_url},
                ))
                n += 1

        # 2) Scripts in scripts/ folder with //help-url: header (e.g. ApiSamples)
        scripts_dir = pkg_dir / "scripts"
        if scripts_dir.is_dir():
            for sf in scripts_dir.rglob("*"):
                if not sf.is_file() or sf.suffix not in {".js", ".py", ".r", ".jl", ".m"}:
                    continue
                try:
                    head_lines: list[str] = []
                    with sf.open(encoding="utf-8", errors="replace") as f:
                        for line in f:
                            s = line.strip()
                            if s.startswith("//") or s.startswith("#"):
                                head_lines.append(line.rstrip("\n"))
                            elif not s:
                                continue
                            else:
                                break
                except OSError:
                    continue
                if not head_lines:
                    continue
                ann = c.parse_annotation_block(head_lines)
                help_url = ann.get("help_url")
                if not help_url:
                    continue
                norm = normalize_help_url(help_url)
                if not norm:
                    continue
                url, anc = norm
                tgt = anchor_id_for(url, anc) if anc else doc_id_for(url)
                # Without an ApiSample extractor yet, point at the Script we already extract
                script_name = ann.get("name") or sf.stem
                bundle.add(Documents(
                    from_id=tgt,
                    to_id=c.script_id(pkg_dir.name, script_name),
                    derivation="script_help_url",
                    derived_by=Provenance.ANNOTATION,
                    extras={"raw": help_url, "script_path": c.repo_rel(sf)},
                ))
                n += 1

        # 3) `this.helpUrl = '...'` in src/**/*.ts (skip .g.ts)
        src_dir = pkg_dir / "src"
        if src_dir.is_dir():
            # Pre-build a lookup of viewer-role functions in this package,
            # keyed by the canonical name -> id, plus a slug -> id index.
            # We scan package.g.ts again to gather names quickly.
            viewer_funcs: dict[str, str] = {}
            if g_path.is_file():
                gtext = g_path.read_text(encoding="utf-8", errors="replace")
                for s, e, blk, anc_line in c.find_annotation_blocks(gtext):
                    a = c.parse_annotation_block(blk)
                    if not a.get("name"):
                        continue
                    role = (a.get("meta", {}) or {}).get("role")
                    is_viewer = role == "viewer" or any(
                        (o.get("type") == "viewer") for o in a.get("outputs", []))
                    if not is_viewer:
                        continue
                    fid = c.func_id(pkg_dir.name, a["name"])
                    viewer_funcs[a["name"].lower().replace(" ", "")] = fid
                    viewer_funcs[_slugify(a["name"]).replace("-", "")] = fid

            for ts_file in src_dir.rglob("*.ts"):
                if ts_file.name.endswith(".g.ts") or "node_modules" in ts_file.parts:
                    continue
                try:
                    text = ts_file.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                if "this.helpUrl" not in text:
                    continue
                # Collect class boundaries to associate the assignment with a class
                class_positions = [(m.start(), m.group("cls")) for m in CLASS_DECL_RE.finditer(text)]
                for hm in THIS_HELP_URL_RE.finditer(text):
                    raw = hm.group("url")
                    norm = normalize_help_url(raw)
                    if not norm:
                        continue
                    url, anc = norm
                    tgt = anchor_id_for(url, anc) if anc else doc_id_for(url)
                    # Find enclosing class
                    cls = None
                    for pos, name in class_positions:
                        if pos < hm.start():
                            cls = name
                        else:
                            break
                    fid: str | None = None
                    if cls:
                        # Try fuzzy match: class FooBarViewer -> viewer named "Foo Bar"
                        bare = re.sub(r"Viewer$", "", cls).lower()
                        fid = viewer_funcs.get(bare) or viewer_funcs.get(cls.lower())
                    if fid:
                        bundle.add(Documents(
                            from_id=tgt, to_id=fid,
                            derivation="viewer_class_help_url",
                            derived_by=Provenance.AST,
                            extras={"raw": raw, "class": cls,
                                    "file": c.repo_rel(ts_file)},
                        ))
                    else:
                        bundle.add(Documents(
                            from_id=tgt, to_id=c.pkg_id(pkg_dir.name),
                            derivation="viewer_class_help_url_unmatched",
                            derived_by=Provenance.AST,
                            extras={"raw": raw, "class": cls,
                                    "file": c.repo_rel(ts_file)},
                        ))
                    n += 1
    print(f"[{EXTRACTOR_NAME}] {n} help-url Documents edges (annotation + script + viewer-class)")
    return n


# ---------------------------------------------------------------------------
# Plugin README walker
# ---------------------------------------------------------------------------

def _emit_plugin_readmes(repo_root: Path, bundle,
                          url_to_doc: dict[str, DocPage],
                          package_filter: list[str] | None) -> int:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return 0
    n = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        for fname, kind in (("README.md", "readme"),
                             ("META.md", "meta"),
                             ("CLAUDE.md", "claude-guide"),
                             ("CONTRIB.md", "contrib"),
                             ("CREDITS.md", "credits")):
            f = pkg_dir / fname
            if not f.is_file():
                continue
            try:
                text = f.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            fm, body = _parse_frontmatter(text)
            title = fm.get("title") or (
                H1_RE.search(body).group(1).strip() if H1_RE.search(body) else f"{pkg_dir.name} {kind}"
            )
            url = f"/packages/{pkg_dir.name}/{fname.lower()}"
            page = DocPage(
                id=doc_id_for(url),
                name=title, source_layer=SourceLayer.DOCS,
                title=title, audience="developer", doc_kind=kind,
                public=True, frontmatter=fm if isinstance(fm, dict) else {},
                paths=[c.file_ref(f, role=kind)],
                extras={"url": url, "package": pkg_dir.name},
            )
            bundle.add(page)
            url_to_doc[url] = page
            # Also emit a Documents edge: this README documents the Package
            bundle.add(Documents(
                from_id=page.id, to_id=c.pkg_id(pkg_dir.name),
                derivation="readme",
                derived_by=Provenance.ANNOTATION,
            ))
            n += 1
    print(f"[{EXTRACTOR_NAME}] {n} plugin doc files (README/META/CLAUDE/...)")
    return n


# ---------------------------------------------------------------------------
# Tutorial parser
# ---------------------------------------------------------------------------

CLASS_RE = re.compile(
    r"export\s+class\s+(?P<cls>[A-Z][\w]*)\s+extends\s+Tutorial\s*\{",
)
NAME_GETTER_RE = re.compile(
    r"get\s+name\s*\(\s*\)\s*\{[^}]*return\s+['\"](?P<name>[^'\"]+)['\"]",
    re.DOTALL,
)
DESC_GETTER_RE = re.compile(
    r"get\s+description\s*\(\s*\)\s*\{(?P<body>(?:[^{}]|\{[^{}]*\})*)\}",
    re.DOTALL,
)
STEPS_GETTER_RE = re.compile(
    r"get\s+steps\s*\(\s*\)\s*\{[^}]*return\s+(\d+)",
    re.DOTALL,
)
HELP_URL_FIELD_RE = re.compile(
    r"helpUrl\s*[:=]\s*['\"](?P<url>[^'\"]+)['\"]",
)
PREREQS_RE = re.compile(
    r"prerequisites\s*[:=]\s*[^=]*=?\s*\{\s*packages\s*:\s*\[(?P<list>[^\]]*)\]",
    re.DOTALL,
)
TRACK_CTOR_RE = re.compile(
    r"new\s+Track\s*\(\s*['\"](?P<name>[^'\"]+)['\"]\s*,\s*[^,]+,\s*['\"](?P<url>[^'\"]+)['\"]\s*\)",
    re.DOTALL,
)


def _parse_string_concat(body: str) -> str:
    """Crude: pull every string literal out of a getter body and concat."""
    parts = re.findall(r"['\"]([^'\"]+)['\"]", body)
    return " ".join(parts).strip()


def _parse_string_list(s: str) -> list[str]:
    return [x.strip().strip('"').strip("'") for x in s.split(",") if x.strip()]


def _emit_tutorial(ts_file: Path, pkg_name: str, bundle,
                    url_to_doc: dict[str, DocPage]) -> int:
    text = ts_file.read_text(encoding="utf-8", errors="replace")
    n = 0
    for cm in CLASS_RE.finditer(text):
        cls = cm.group("cls")
        # Limit search window to the file from the class onward
        scope = text[cm.start():]
        name = (NAME_GETTER_RE.search(scope) or {}).group("name") if NAME_GETTER_RE.search(scope) else cls
        desc_m = DESC_GETTER_RE.search(scope)
        desc = _parse_string_concat(desc_m.group("body")) if desc_m else None
        steps_m = STEPS_GETTER_RE.search(scope)
        steps = int(steps_m.group(1)) if steps_m else None
        url_m = HELP_URL_FIELD_RE.search(scope)
        help_url_raw = url_m.group("url") if url_m else None
        prereq_m = PREREQS_RE.search(scope)
        prereqs = _parse_string_list(prereq_m.group("list")) if prereq_m else []

        tid = f"tutorial:{pkg_name}:{cls}"
        t = Tutorial(
            id=tid, name=name, source_layer=SourceLayer.DOCS,
            package_id=c.pkg_id(pkg_name),
            class_name=cls,
            help_url=normalize_help_url(help_url_raw)[0] if help_url_raw and normalize_help_url(help_url_raw) else None,
            prerequisites=prereqs,
            step_count=steps,
            description=desc,
            description_provenance=Provenance.AST if desc else None,
            paths=[c.file_ref(ts_file, role="definition")],
            extras={"raw_help_url": help_url_raw},
        )
        bundle.add(t)
        # Documents edge: tutorial's helpUrl points at a DocPage/HelpAnchor
        if help_url_raw:
            norm = normalize_help_url(help_url_raw)
            if norm:
                url, anc = norm
                tgt = anchor_id_for(url, anc) if anc else doc_id_for(url)
                bundle.add(Documents(
                    from_id=tgt, to_id=tid,
                    derivation="tutorial_help_url",
                    derived_by=Provenance.AST,
                    extras={"raw": help_url_raw},
                ))
        n += 1
    return n


def _emit_tracks(ts_file: Path, pkg_name: str, bundle) -> int:
    text = ts_file.read_text(encoding="utf-8", errors="replace")
    n = 0
    for tm in TRACK_CTOR_RE.finditer(text):
        track_name = tm.group("name")
        track_url = tm.group("url")
        norm = normalize_help_url(track_url)
        tid = f"track:{pkg_name}:{c.lib_id(track_name).split(':',1)[1] if ':' in track_name else _slugify(track_name)}"
        tt = TutorialTrack(
            id=tid, name=track_name, source_layer=SourceLayer.DOCS,
            package_id=c.pkg_id(pkg_name),
            help_url=norm[0] if norm else None,
            paths=[c.file_ref(ts_file, role="definition")],
            extras={"raw_help_url": track_url},
        )
        bundle.add(tt)
        bundle.add(HasTutorial(
            from_id=c.pkg_id(pkg_name), to_id=tid,
            derived_by=Provenance.AST,
        ))
        n += 1
    return n


def _walk_tutorials(repo_root: Path, bundle,
                     url_to_doc: dict[str, DocPage],
                     package_filter: list[str] | None) -> tuple[int, int]:
    """Walk the Tutorials package + any plugin's `tutorials/` folder."""
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return 0, 0
    n_t = n_tk = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        # Tutorials package has tracks under src/tracks/
        # Other packages may have a `tutorials/` folder
        candidate_dirs = [pkg_dir / "src" / "tracks", pkg_dir / "tutorials",
                          pkg_dir / "src" / "tutorials"]
        for cand in candidate_dirs:
            if not cand.is_dir():
                continue
            for ts_file in cand.rglob("*.ts"):
                if "node_modules" in ts_file.parts:
                    continue
                try:
                    n_t += _emit_tutorial(ts_file, pkg_dir.name, bundle, url_to_doc)
                    n_tk += _emit_tracks(ts_file, pkg_dir.name, bundle)
                except Exception as e:
                    print(f"  ! tutorial parse {ts_file}: {e}")
    print(f"[{EXTRACTOR_NAME}] {n_t} tutorials, {n_tk} tracks")
    return n_t, n_tk


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    url_to_doc: dict[str, DocPage] = {}
    _walk_help(repo_root, bundle, url_to_doc)
    _emit_plugin_readmes(repo_root, bundle, url_to_doc, package_filter)
    _emit_help_url_edges(repo_root, bundle, url_to_doc, package_filter)
    _walk_tutorials(repo_root, bundle, url_to_doc, package_filter)
