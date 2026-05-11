"""Extract ChangelogEntry entities from `packages/*/CHANGELOG.md`.

Format (per root `CLAUDE.md`):
  ## v.next                              # unreleased
  * GROK-12345: Did something
  * Fixed Y

  ## 1.17.5 (2026-04-15)
  * GROK-23456: Description

Emits:
  - ChangelogEntry
  - HasChangelogEntry (Package -> ChangelogEntry)
  - JiraTicket  (placeholder node — only id+key, fleshed out by jira extractor)
  - MentionsTicket (ChangelogEntry -> JiraTicket)
"""

from __future__ import annotations

import re
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import ChangelogEntry, JiraTicket
from schema.relations import HasChangelogEntry, MentionsTicket

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_changelog"

VERSION_RE = re.compile(r"^##\s+(?P<ver>v\.next|\d+\.\d+\.\d+(?:-[\w.]+)?)\s*(?:\((?P<date>\d{4}-\d{2}-\d{2})\))?\s*$")
TICKET_RE = re.compile(r"GROK-\d+")
BULLET_RE = re.compile(r"^\s*(?:[-*+]\s+|\d+\.\s+)(?P<text>.+)$")


def _ticket_from(text: str) -> str | None:
    m = TICKET_RE.search(text)
    return m.group(0) if m else None


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return
    total = 0
    seen_tickets: set[str] = set()
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        cl = pkg_dir / "CHANGELOG.md"
        if not cl.is_file():
            continue
        try:
            text = cl.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

        cur_ver: str | None = None
        cur_date: str | None = None
        cur_idx = 0
        had_bullet_for_version = False
        for raw_line in text.splitlines():
            line = raw_line.rstrip()
            mv = VERSION_RE.match(line)
            if mv:
                cur_ver = mv.group("ver")
                cur_date = mv.group("date")
                cur_idx = 0
                had_bullet_for_version = False
                continue
            if not line.strip():
                continue
            if cur_ver is None:
                continue
            mb = BULLET_RE.match(line)
            if mb:
                bullet_text = mb.group("text").strip()
                had_bullet_for_version = True
            elif not had_bullet_for_version and not line.startswith(("#", "<!--")):
                # Paragraph-style entry: capture as a single bullet equivalent.
                # Only take the first paragraph line per version (avoids
                # multi-paragraph spam), then ignore until the next version.
                bullet_text = line.strip(" *")
                had_bullet_for_version = True
            else:
                continue
            ticket = _ticket_from(bullet_text)
            ent_id = c.clog_id(pkg_dir.name, cur_ver, cur_idx)
            entry = ChangelogEntry(
                id=ent_id, name=bullet_text[:80],
                source_layer=SourceLayer.PROCESS,
                package_id=c.pkg_id(pkg_dir.name),
                version=cur_ver,
                released_at=cur_date,
                ticket_id=ticket,
                description=bullet_text,
                description_provenance=Provenance.ANNOTATION,
                paths=[c.file_ref(cl, role="changelog")],
            )
            bundle.add(entry)
            bundle.add(HasChangelogEntry(
                from_id=c.pkg_id(pkg_dir.name), to_id=ent_id,
                derived_by=Provenance.ANNOTATION))
            if ticket:
                jid = c.jira_id(ticket)
                if ticket not in seen_tickets:
                    seen_tickets.add(ticket)
                    bundle.add(JiraTicket(
                        id=jid, name=ticket, source_layer=SourceLayer.PROCESS,
                        ticket_key=ticket,
                        description_provenance=Provenance.ANNOTATION,
                    ))
                bundle.add(MentionsTicket(
                    from_id=ent_id, to_id=jid,
                    derived_by=Provenance.ANNOTATION))
            cur_idx += 1
            total += 1
    print(f"[{EXTRACTOR_NAME}] {total} changelog entries, {len(seen_tickets)} tickets")
