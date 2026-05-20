# Scenario 5 — Governance in Motion

**Theme:** Data access democratization + enterprise governance
**Length target:** 3 min
**Audience lean:** IT/platform leads (primary), data engineers, compliance/QA
**Recording risk:** Medium — needs two accounts and audit-log surfacing

This is the **IT-buyer scenario.** It exists because every IT lead who watches S1 will ask the same question: *"What stops the agent from leaking data?"* This scenario answers it on screen.

## The setup that makes this scenario unique

Two demo user accounts, both signed in (or quickly switched between):

- **`alice@demo`** — broad access. Member of `clinical-research-all`.
- **`bob@demo`** — project-limited. Member of `project-alpha-only`. Cannot see other trial data.

The point of the scenario is: **same prompt, same agent, different result by user.**

## Prompt (run twice)

> Show me all clinical trial outcomes for Drug X with patient demographics. Then share the dashboard with the Clinical-Ops team.

## Expected sequence — as Alice

1. `Query` provider routes to the curated `TrialOutcomesByCompound` template.
2. Full result set returned (e.g., 4 trials, 1800 patients).
3. Grokky builds a small dashboard: outcomes table + demographics summary + a viewer.
4. **Share confirmation.** Before the irreversible share action, Grokky asks in chat: *"About to share this dashboard with `clinical-ops` (23 members). Confirm?"* Alice replies `yes`. (See T1.6 in POLISH.md — safety-net confirmation is a system-prompt addition today; inline buttons are roadmap.)
5. Project created and shared. A shareable URL appears in chat.
6. **A small "audit" badge appears** in the chat indicating the action was logged.

## Expected sequence — as Bob

1. Same prompt. `Query` provider routes to the same curated template (verified: `query-matching.ts:54-57` only invokes registered queries with `meta.searchPattern`; no ad-hoc SQL).
2. The query executes under Bob's identity (verified token pass-through: `runtime-client.ts:135` → `claude-runtime/src/server.ts:149` → `mcp-server/src/shared-api-client.ts:15` enforces `x-user-api-key`). Row-level / project-level ACLs filter results to Project Alpha only (e.g., 1 trial, 200 patients).
3. **Today (T3.1 not yet in place):** the query returns 1 row; the agent has no signal that more rows existed pre-filter (verified absent — no `filtered_by_acl` flag in query results). The presenter voices it on stage: *"Bob is only in `project-alpha-only`; the platform filtered the other trials automatically — same way they'd be filtered if Bob ran the SQL directly."* **With T3.1:** Grokky says explicitly: *"You have access to 1 of 4 trials matching this query. Trials in projects you don't belong to (Beta, Gamma, Delta) were excluded. Contact your project admin for broader access."*
4. The partial dashboard is built. Same share confirmation beat fires; Bob confirms.
5. Sharing succeeds, but only with the subset Bob can see.
6. Audit badge appears (with T3.2 — today, audit entries exist server-side per `help/govern/audit/audit.md`, just not surfaced in chat).

## Wow moment

Two split-screen browser windows side by side, both running the same prompt. **Different results. No agent jailbreak, no over-share, no error.** The audience sees governance working without anyone having to claim it.

## Talking points

- **IT lead:** "The agent runs as the user. Permissions are enforced by the platform, not by the prompt. There's no prompt-injection class that lets a user see data they couldn't already reach via SQL — the agent has no credentials of its own."
- **IT lead (clarification):** "The share API won't fire until you confirm — the server refuses the call without your explicit yes. This isn't a rule the agent follows; it's a check at the server, and there's no path the agent can use to skip it. Same model as a PR that won't merge without approval."
- **Compliance / QA:** "Every agent action logs the user, timestamp, query template, parameters, and recipients. Same audit trail as anything else on the platform."
- **Data engineer:** "The agent uses your existing curated query templates. Whatever governance those carry, the agent inherits."

## Works today

- Datagrok ACL/permission system enforces at the platform level. The agent has no way around it.
- Curated query templates exist and are matched by the `Query` provider.
- Datagrok's audit log captures action provenance.
- Group sharing exists.

## Needs polish (see POLISH.md)

- **Graceful-degradation messaging** (T3.1) — the LLM needs to detect "result set was filtered by ACL" and surface that clearly. Today it might silently return the partial result without explaining. Add a structured response pattern.
- **Audit-log surfacing UI** (T3.2) — there should be a one-click "show audit entry for this action" affordance in the chat. May need a small UI addition to the chat panel.
- **Mechanical confirmation gate** (T1.6) — the load-bearing dependency for this entire scenario. The IT-lead talking point is only honest if L1 (gated MCP tools refusing without a confirmation token) and L2 (`datagrok-exec` sandbox closing the side door) are both shipped. **If T1.6 isn't shipped, cut S5 from the demo set** — the soft-form-only version was rejected as "ask and pray."
- **Demo-account preset** (T4.1) — Alice and Bob must exist with the right group memberships in the demo server before the take. Document this in a setup script.

## Backup plan

- If running side-by-side browsers is fragile, record Alice and Bob takes separately and cut them together with a clean split-screen edit.
- If audit-log surfacing isn't shipped, open the admin audit-log page in a separate tab and switch to it manually during the live demo.
- If the graceful-degradation message is rough, the presenter can voice it: "Notice Bob got the project-restricted subset — the platform filtered Beta, Gamma, Delta automatically." Still lands.

## Why this scenario closes deals

S1 makes scientists excited. S5 makes IT leads sign. Without S5, every S1 viewer asks the same anxious question. With S5, that question is pre-answered.
