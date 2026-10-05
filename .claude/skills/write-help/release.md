# Updating the help for a release

Part of the `write-help` skill. The job: find what shipped, decide what the help needs and what
only the community post needs, then write and verify it (SKILL.md sections 2–6).

## Scope: what shipped

Collect one list with a ticket key per item:

- Jira: `project = GROK AND fixVersion = "<version>"`. Take all types and statuses, keys and
  summaries first. Bugs matter too: a fix can make the help wrong.
- The release plan in the private core repository (`core/docs/release/<version>.yaml`), if you
  have access.
- Commits since the previous release: in core, from the previous release branch point; in this
  repository, `git log --since=<date of the previous release> origin/master`.
- Plugin `CHANGELOG.md` files: dated versions are released, `## v.next` is not.
- The community post drafts, if the team keeps them: they show what matters to users.

For each item, with evidence:

| Question | How | If "no" |
|---|---|---|
| In master? | `git log origin/master --grep=<KEY>`, the code itself | Leave out, report where it is |
| New in this version? | Core: absent from the previous release branch (`git grep <symbol> origin/release/<previous>`, private core repo). Plugins: first appears in a version dated after the previous release | Not new. Document only if the help lacks it |
| User-facing? | A UI label, menu item, or property | Leave out |
| Finished and released? | Ticket status, plugin version | Ask. Don't document unfinished work |
| In the running product? | The stand the team reviews on | Report as not deployed |

Split a large audit by area (core viewers, core app, server and deploy, plugins, Jira) and run it
in parallel with read-only agents. Each returns, per item: in master (evidence, exact labels,
defaults), help coverage (page#anchor: covered, partial, missing, or wrong), and the best place.
Agent claims are leads: re-check everything you will write.

## Triage: help or community post

Document in the help, in this order:

1. Behaviour that changed so the help is **wrong**: a menu path that no longer exists, renamed
   commands, dialog fields that changed, removed options.
2. Capabilities the help doesn't describe at all, even old ones.
3. Settings whose meaning isn't obvious from their name in the properties table.

Leave to the community post:

- A capability every viewer shares, when adding it to one page would duplicate it.
- Small UX improvements that need no instructions.
- Properties already in the generated properties table and self-explanatory.
- Changes that affect existing results (recalculated statistics, models that must be retrained,
  upgraded runtimes): these are announcements, and belong in the post first.

For every item left out of the help, check the post covers it. If it doesn't, propose a note in
the post's style and add it only after the reviewer agrees.

Admin-facing changes (users, permissions, deployment, configuration) go into a separate commit,
so an administrator can review them on their own.

## Report

A candidates table with, per item: area, change, evidence (sha or key, date), UI labels, help
now (page#anchor and state), recommendation. List separately: help that describes unreleased
work, items waiting for owner confirmation, community-post-only items, unfinished items, and
exclusions decided by the reviewer.
