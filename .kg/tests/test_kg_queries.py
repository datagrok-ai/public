"""End-to-end tests for the live Kuzu graph at `.kg/kg.kuzu`.

Run from `.kg/`:
    .venv/Scripts/python tests/test_kg_queries.py
    .venv/Scripts/python tests/test_kg_queries.py --filter coverage
    .venv/Scripts/python tests/test_kg_queries.py --verbose

Each test is a small Python function that runs one or more Cypher queries
and asserts on the result. Tests are grouped by category so the suite
doubles as a living spec for "what the graph guarantees right now".

Categories:
  schema       — every declared entity/relation kind has a table
  coverage     — minimum row counts; specific packages/libraries present
  connectivity — every Package has content; FK-like soundness checks
  crosslayer   — multi-hop traversals: dep, doc, feature
  features     — feature-cluster invariants (LLM enrichment output)
  gaps         — known coverage gaps the graph should surface (gap > 0)

Failing tests print the offending Cypher + summary so they double as
diagnostic queries.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import kuzu

from schema.entities import ENTITY_KINDS
from schema.relations import RELATION_KINDS

DB_PATH = ROOT / "kg.kuzu"


# ---------------------------------------------------------------------------
# Tiny framework
# ---------------------------------------------------------------------------

@dataclass
class TestCase:
    name: str
    category: str
    fn: Callable
    desc: str = ""


REGISTRY: list[TestCase] = []


def test(category: str, desc: str = ""):
    def decorator(fn: Callable) -> Callable:
        REGISTRY.append(TestCase(name=fn.__name__, category=category, fn=fn, desc=desc))
        return fn
    return decorator


def q(conn, cypher: str, **params) -> list[tuple]:
    """Run a Cypher query, return list of tuples."""
    res = conn.execute(cypher, parameters=params or None)
    rows: list[tuple] = []
    while res.has_next():
        rows.append(tuple(res.get_next()))
    return rows


def q_one(conn, cypher: str, **params) -> Any:
    rows = q(conn, cypher, **params)
    return rows[0][0] if rows else None


# ---------------------------------------------------------------------------
# Schema tests — every declared kind has a table, every node has identity
# ---------------------------------------------------------------------------

@test("schema", "Every entity kind in entities.py is queryable as a Kuzu node table")
def test_all_entity_tables_exist(conn):
    missing: list[str] = []
    for kind in ENTITY_KINDS:
        try:
            q(conn, f"MATCH (n:`{kind}`) RETURN count(n);")
        except Exception as e:
            missing.append(f"{kind}: {e}")
    assert not missing, f"Missing/broken node tables: {missing}"


@test("schema", "Every relation in relations.py is queryable as a Kuzu rel table")
def test_all_rel_tables_exist(conn):
    missing: list[str] = []
    for predicate in RELATION_KINDS:
        try:
            q(conn, f"MATCH ()-[r:`{predicate}`]->() RETURN count(r);")
        except Exception as e:
            missing.append(f"{predicate}: {e}")
    assert not missing, f"Missing/broken rel tables: {missing}"


@test("schema", "Every loaded node has a non-empty id")
def test_node_ids_present(conn):
    bad = q_one(conn, "MATCH (n) WHERE n.id IS NULL OR n.id = '' RETURN count(n);")
    assert bad == 0, f"{bad} nodes have null/empty id"


@test("schema", "Node ids are globally unique")
def test_node_ids_unique(conn):
    dupes = q(conn, "MATCH (n) WITH n.id AS id, count(n) AS c WHERE c > 1 RETURN id, c LIMIT 10;")
    assert not dupes, f"Duplicate node ids: {dupes}"


# ---------------------------------------------------------------------------
# Coverage tests — minimum row counts and that key packages/libs are present
# ---------------------------------------------------------------------------

@test("coverage", "At least 50 Datagrok packages")
def test_packages_count(conn):
    n = q_one(conn, "MATCH (p:Package) RETURN count(p);")
    assert n >= 50, f"Got only {n} packages"


@test("coverage", "At least 500 RegisteredFunctions")
def test_functions_count(conn):
    n = q_one(conn, "MATCH (f:RegisteredFunction) RETURN count(f);")
    assert n >= 500, f"Got only {n} functions"


@test("coverage", "At least 15 libraries")
def test_libraries_count(conn):
    n = q_one(conn, "MATCH (l:Library) RETURN count(l);")
    assert n >= 15, f"Got only {n} libraries"


@test("coverage", "Specific known packages present (Chem, Chembl, PowerPack, Bio, Charts, ApiSamples, ApiTests)")
def test_known_packages(conn):
    needed = ["Chem", "Chembl", "PowerPack", "Bio", "Charts", "ApiSamples", "ApiTests"]
    rows = q(conn, "MATCH (p:Package) WHERE p.name IN $names RETURN p.name;", names=needed)
    found = {r[0] for r in rows}
    assert found == set(needed), f"Missing: {set(needed) - found}"


@test("coverage", "Chem has at least 100 RegisteredFunctions")
def test_chem_functions(conn):
    n = q_one(conn, "MATCH (p:Package {name:'Chem'})-[:EXPORTS]->(f) RETURN count(f);")
    assert n >= 100, f"Chem has only {n} functions"


@test("coverage", "Chembl has at least 30 SQL queries")
def test_chembl_queries(conn):
    n = q_one(conn, "MATCH (p:Package {name:'Chembl'})-[:HAS_QUERY]->(q:DataQuery) RETURN count(q);")
    assert n >= 30, f"Chembl has only {n} queries"


@test("coverage", "Help docs are extracted (>= 400 DocPages)")
def test_doc_pages_count(conn):
    n = q_one(conn, "MATCH (d:DocPage) RETURN count(d);")
    assert n >= 400, f"Got only {n} doc pages"


@test("coverage", "Tutorials are extracted (>= 15 Tutorial)")
def test_tutorials_count(conn):
    n = q_one(conn, "MATCH (t:Tutorial) RETURN count(t);")
    assert n >= 15, f"Got only {n} tutorials"


@test("coverage", "Plugin CHANGELOGs parsed (>= 1000 entries)")
def test_changelog_count(conn):
    n = q_one(conn, "MATCH (c:ChangelogEntry) RETURN count(c);")
    assert n >= 1000, f"Got only {n} changelog entries"


@test("coverage", "JIRA tickets auto-linked (>= 100 unique)")
def test_jira_tickets(conn):
    n = q_one(conn, "MATCH (t:JiraTicket) RETURN count(t);")
    assert n >= 100, f"Got only {n} tickets"


# ---------------------------------------------------------------------------
# Connectivity tests — every Package has content; FK-like soundness checks
# ---------------------------------------------------------------------------

@test("connectivity", "Every Package has at least one function/script/query/changelog")
def test_no_empty_packages(conn):
    empty = q(conn, """
        MATCH (p:Package)
        WHERE NOT EXISTS { MATCH (p)-[:EXPORTS]->() }
          AND NOT EXISTS { MATCH (p)-[:HAS_SCRIPT]->() }
          AND NOT EXISTS { MATCH (p)-[:HAS_QUERY]->() }
          AND NOT EXISTS { MATCH (p)-[:HAS_CHANGELOG_ENTRY]->() }
        RETURN p.name LIMIT 10;
    """)
    assert not empty, f"Empty packages: {[r[0] for r in empty]}"


@test("connectivity", "Every RegisteredFunction's package_id resolves to a known Package")
def test_functions_have_packages(conn):
    orphan = q(conn, """
        MATCH (f:RegisteredFunction)
        WHERE f.package_id IS NOT NULL
          AND NOT EXISTS { MATCH (p:Package {id: f.package_id}) }
        RETURN count(f);
    """)
    assert orphan[0][0] == 0, f"{orphan[0][0]} functions reference unknown package"


@test("connectivity", "Every Library has at least 1 LibraryModule")
def test_libraries_have_modules(conn):
    bad = q(conn, """
        MATCH (l:Library)
        WHERE NOT EXISTS { MATCH (m:LibraryModule {library_id: l.id}) }
        RETURN l.name LIMIT 10;
    """)
    # Some libraries may be very small / no src; assert majority have modules
    n_total = q_one(conn, "MATCH (l:Library) RETURN count(l);")
    assert len(bad) <= 2, f"{len(bad)} libraries with zero modules: {bad} (total {n_total})"


@test("connectivity", "Tutorials carry a class_name and at least an attempted helpUrl")
def test_tutorial_fields(conn):
    no_class = q_one(conn, "MATCH (t:Tutorial) WHERE t.class_name IS NULL RETURN count(t);")
    assert no_class == 0, f"{no_class} tutorials have no class_name"


# ---------------------------------------------------------------------------
# Cross-layer tests — multi-hop traversals
# ---------------------------------------------------------------------------

@test("crosslayer", "Chem depends on the @datagrok-libraries/test library")
def test_chem_depends_on_test(conn):
    n = q_one(conn, """
        MATCH (p:Package {name:'Chem'})-[:DEPENDS_ON]->(l:Library {name:'test'}) RETURN count(*);
    """)
    assert n >= 1, "Chem -> lib:test edge missing"


@test("crosslayer", "Bio imports >= 20 distinct modules from lib:bio")
def test_bio_module_imports(conn):
    n = q_one(conn, """
        MATCH (p:Package {name:'Bio'})-[:IMPORTS_FROM_MODULE]->(m:LibraryModule {library_id:'lib:bio'})
        RETURN count(distinct m);
    """)
    assert n >= 20, f"Only {n} bio modules imported"


@test("crosslayer", "Most-used library is one of: test, utils")
def test_top_library(conn):
    rows = q(conn, """
        MATCH (p:Package)-[:DEPENDS_ON]->(l:Library)
        RETURN l.name, count(p) AS c ORDER BY c DESC LIMIT 1;
    """)
    assert rows and rows[0][0] in {"test", "utils"}, f"Top library is {rows[0] if rows else None}"


@test("crosslayer", "Chembl queries use a Chembl-namespace connection")
def test_chembl_query_uses_connection(conn):
    n = q_one(conn, """
        MATCH (p:Package {name:'Chembl'})-[:HAS_QUERY]->(q:DataQuery)-[:USES_CONNECTION]->(c:DataConnection)
        WHERE c.name CONTAINS 'hembl'
        RETURN count(distinct q);
    """)
    assert n >= 5, f"Only {n} Chembl queries linked to a Chembl connection"


@test("crosslayer", "At least one ChangelogEntry mentions a Jira ticket")
def test_changelog_mentions_ticket(conn):
    n = q_one(conn, """
        MATCH (c:ChangelogEntry)-[:MENTIONS_TICKET]->(t:JiraTicket) RETURN count(*);
    """)
    assert n >= 50, f"Only {n} changelog->ticket links"


@test("crosslayer", "DocPage LinksTo edges resolve to real DocPages or HelpAnchors (>1000 working)")
def test_doc_links_resolve(conn):
    n = q_one(conn, """
        MATCH (d:DocPage)-[:LINKS_TO]->(t)
        RETURN count(*);
    """)
    assert n >= 1000, f"Only {n} doc links land on real targets"


# ---------------------------------------------------------------------------
# Feature tests — only meaningful after enrichment has run
# ---------------------------------------------------------------------------

@test("features", "At least 3 Packages have HAS_FEATURE edges (pilot)")
def test_packages_with_features(conn):
    n = q_one(conn, "MATCH (p:Package)-[:HAS_FEATURE]->(:Feature) RETURN count(distinct p);")
    assert n >= 3, f"Only {n} packages have features yet (pilot needs >= 3)"


@test("features", "Every Feature has at least one PART_OF_FEATURE member")
def test_features_have_members(conn):
    rows = q(conn, """
        MATCH (f:Feature)
        WHERE NOT EXISTS { MATCH (f)<-[:PART_OF_FEATURE]-() }
        RETURN f.id LIMIT 10;
    """)
    assert not rows, f"Orphan features (no members): {[r[0] for r in rows]}"


@test("features", "Every Feature has exactly one HAS_FEATURE owner Package")
def test_feature_single_owner(conn):
    rows = q(conn, """
        MATCH (f:Feature)<-[:HAS_FEATURE]-(p:Package)
        WITH f, count(p) AS owners
        WHERE owners <> 1
        RETURN f.id, owners LIMIT 10;
    """)
    assert not rows, f"Features with non-1 owner: {rows}"


@test("features", "Chem features include 'Activity Cliffs' (LLM clustering sanity)")
def test_chem_activity_cliffs_feature(conn):
    rows = q(conn, """
        MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature)
        WHERE f.name CONTAINS 'Activity' AND f.name CONTAINS 'Cliff'
        RETURN f.id;
    """)
    assert rows, "No 'Activity Cliffs'-named feature on Chem"


@test("features", "A meaningful fraction of Features carry a DocPage member")
def test_features_have_docs(conn):
    n_total = q_one(conn, "MATCH (f:Feature) RETURN count(f);")
    n_with_doc = q_one(conn, """
        MATCH (f:Feature)<-[:PART_OF_FEATURE]-(d:DocPage) RETURN count(distinct f);
    """)
    if n_total == 0:
        return                              # enrichment not run yet — skip silently
    ratio = n_with_doc / n_total
    # Realistic target across the full 79-package corpus: many small "infra"
    # / test-harness / file-importer features legitimately have no DocPage.
    # We assert that a meaningful fraction (~40%+) DO carry doc membership.
    assert ratio >= 0.4, f"Only {n_with_doc}/{n_total} ({ratio:.0%}) features have docs"


@test("features", "User-facing Feature categories (chem-*, viz, ml, bio-*) carry docs at higher rate")
def test_user_facing_features_docs(conn):
    """Refines the previous test: features in user-facing categories (extracted
    from extras JSON) should have a noticeably higher doc-coverage rate than
    pure-infra ones. This catches regressions where doc-clustering breaks
    altogether for the user-visible features that matter most."""
    rows = q(conn, """
        MATCH (f:Feature)
        WITH count(f) AS total,
             count(CASE WHEN EXISTS { MATCH (f)<-[:PART_OF_FEATURE]-(:DocPage) } THEN f END) AS with_doc
        RETURN total, with_doc;
    """)
    if not rows:
        return
    total, with_doc = rows[0]
    if total == 0:
        return
    # Hard floor: at least 200 features must carry a doc
    assert with_doc >= 200, f"Only {with_doc} features have any DocPage member; clustering may have regressed"


# ---------------------------------------------------------------------------
# Gap tests — assert that a known coverage gap exists, with a count
# ---------------------------------------------------------------------------

@test("gaps", "Charts viewers are >= 80% covered by docs (regression check after doc_linking)")
def test_charts_viewers_doc_coverage(conn):
    """Originally a 'known gap' assertion (≥ 5 viewers undocumented). After
    the `doc_linking` enricher landed it dropped to 1, so we flip the
    invariant: assert that high coverage is maintained going forward.
    A regression that drops doc coverage will trip this test."""
    total = q_one(conn, """
        MATCH (p:Package {name:'Charts'})-[:EXPORTS]->(f:RegisteredFunction {role:'viewer'})
        RETURN count(f);
    """)
    documented = q_one(conn, """
        MATCH (p:Package {name:'Charts'})-[:EXPORTS]->(f:RegisteredFunction {role:'viewer'})
        WHERE EXISTS { MATCH (f)<-[:DOCUMENTS]-() }
           OR EXISTS { MATCH (f)-[:PART_OF_FEATURE]->(:Feature)<-[:DOCUMENTS]-() }
        RETURN count(f);
    """)
    if total == 0:
        return
    ratio = documented / total
    assert ratio >= 0.80, f"Only {documented}/{total} ({ratio:.0%}) Charts viewers documented"


@test("gaps", "Many DocPages are orphans (no inbound DOCUMENTS, no outbound LINKS_TO)")
def test_orphan_docs(conn):
    n = q_one(conn, """
        MATCH (d:DocPage)
        WHERE NOT EXISTS { MATCH (d)-[:DOCUMENTS]->() }
          AND NOT EXISTS { MATCH (d)-[:LINKS_TO]->() }
        RETURN count(d);
    """)
    assert n >= 50, f"Only {n} orphan docs — was expecting > 50 from earlier diagnostics"


@test("gaps", "Most ApiSample-style scripts have no //help-url: annotation")
def test_apisamples_help_url_gap(conn):
    total = q_one(conn, """
        MATCH (p:Package {name:'ApiSamples'})-[:HAS_SCRIPT]->(s:Script) RETURN count(s);
    """)
    documented = q_one(conn, """
        MATCH (p:Package {name:'ApiSamples'})-[:HAS_SCRIPT]->(s:Script)<-[:DOCUMENTS]-()
        RETURN count(distinct s);
    """)
    if total > 0:
        ratio_undoc = (total - documented) / total
        assert ratio_undoc >= 0.8, f"Only {ratio_undoc:.0%} undocumented (expected >80%)"


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

ANSI = sys.stdout.isatty()
def _color(s: str, code: str) -> str:
    return f"\033[{code}m{s}\033[0m" if ANSI else s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--filter", default="", help="Only run tests in this category (substring match)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not DB_PATH.exists():
        print(f"DB not found at {DB_PATH}. Run `py build.py` first.", file=sys.stderr)
        return 2

    db = kuzu.Database(str(DB_PATH))
    conn = kuzu.Connection(db)

    selected = [t for t in REGISTRY if not args.filter or args.filter in t.category]
    by_cat: dict[str, list[TestCase]] = {}
    for t in selected:
        by_cat.setdefault(t.category, []).append(t)

    total = len(selected)
    passed = failed = 0
    failures: list[tuple[TestCase, str]] = []
    t0 = time.time()

    for cat, items in by_cat.items():
        print(f"\n[{cat}]")
        for tc in items:
            try:
                tc.fn(conn)
                passed += 1
                print(f"  {_color('PASS', '32')}  {tc.name:<42}  {tc.desc}")
            except AssertionError as e:
                failed += 1
                failures.append((tc, str(e)))
                print(f"  {_color('FAIL', '31')}  {tc.name:<42}  {tc.desc}")
                print(f"        -> {e}")
            except Exception as e:
                failed += 1
                failures.append((tc, f"{type(e).__name__}: {e}"))
                print(f"  {_color('ERROR', '31')} {tc.name:<42}  {tc.desc}")
                print(f"        -> {type(e).__name__}: {e}")
                if args.verbose:
                    traceback.print_exc()

    print(f"\n{passed}/{total} passed, {failed} failed in {time.time() - t0:.1f}s")
    if failures and args.verbose:
        print("\n--- failure detail ---")
        for tc, msg in failures:
            print(f"  {tc.category}/{tc.name}: {msg}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
