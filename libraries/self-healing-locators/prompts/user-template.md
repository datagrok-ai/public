# User Message Template — Self-Healing Locators

> The runtime renders this template into the `messages[0].content` of
> the Anthropic API call. Sections wrapped in `{{...}}` are
> substituted; everything else is literal. The template is intended to
> be deterministic — same inputs yield byte-identical prompts so
> caching works.

---

```
<task>
A UI test selector failed. The deterministic resolver tried the
strategies it can attempt automatically and either could not find a
unique element or found one with low confidence. Propose alternative
selectors using the `propose_locators` tool.
</task>

<test_action>
{{action}}
</test_action>

<failed_selector>
{{originalSelector}}
</failed_selector>

<already_tried>
{{#each alreadyTriedSelectors}}
- {{this}}
{{/each}}
</already_tried>

<fingerprint>
{{fingerprintJson}}
</fingerprint>

<accessibility_tree>
{{accessibilityTreeText}}
</accessibility_tree>

{{#if includeHtmlSnippet}}
<html_snippet>
{{htmlSnippet}}
</html_snippet>
{{/if}}

{{#if includeScreenshot}}
<note>
A screenshot of the page region is attached as an image content block.
Use it only as a sanity check; the accessibility tree is the primary
source of truth.
</note>
{{/if}}

<negative_constraints>
- Selectors listed in <already_tried> have been verified to NOT work.
  Do not propose them again or trivial variations of them.
- The fingerprint reflects the element when the test last passed. The
  accessibility tree above reflects the page now. Differences between
  them are exactly what you are reasoning about.
- The platform is Datagrok. Common UI primitives (`ui.button`,
  `ui.input.*`, `ui.dialog`) typically expose ARIA roles and accessible
  names; viewers expose `dg-widget-type` markers. Prefer these over
  raw DOM selectors.
- The test action is "{{action}}". If a candidate element cannot
  receive this action (e.g. a static `<span>` for a click action),
  do not propose it.
</negative_constraints>

<output_instructions>
Call the `propose_locators` tool with up to 5 candidates ranked by
your estimate of stability. If you cannot propose any plausible
candidate, set `unable_to_match: true` with the appropriate
`unable_reason`. Do not write any text outside the tool call.
</output_instructions>
```

---

## Substitution variables

| Variable | Source | Notes |
|---|---|---|
| `action` | `healing.resolve()` options | One of `click`, `type`, `assert`, `read`, `hover`, etc. |
| `originalSelector` | failed call site | The exact selector string that the test passed in. |
| `alreadyTriedSelectors` | runtime resolver | Each candidate the deterministic tiers attempted, in order. Truncated to 20 entries — beyond that is noise. |
| `fingerprintJson` | registry | The fingerprint object, with `contextSnippet` stripped to avoid duplication. JSON.stringify with 2-space indent. |
| `accessibilityTreeText` | snapshot | A YAML-like flattening of the a11y tree, pruned per § "Pruning". |
| `includeHtmlSnippet` | feature flag | True when the a11y tree alone is judged insufficient (e.g. dominant role is `generic`). |
| `htmlSnippet` | snapshot | Sanitized outerHTML of the candidate region, truncated to 4 KB. Inline event handlers and data-* attributes preserved; class hashes left in (the LLM is told to ignore them). |
| `includeScreenshot` | feature flag + classification | True only when screenshots are enabled AND case is classified structural-or-worse. |

## Pruning rules for the accessibility tree

We do not send the entire a11y tree. The pruner walks the tree once
and keeps:

1. The current view's accessibility root.
2. All nodes whose `role` matches the fingerprint's `ariaRole`.
3. All nodes whose `name` is within Levenshtein distance 0.7 of the
   fingerprint's `name`, `ariaLabel`, or `visibleText`.
4. The 2-level neighborhood around each kept node.
5. Nodes whose `dg-widget-type` matches the fingerprint's
   `dgWidgetType`.

Output format is a stable, indented text representation:

```
[role="dialog" name="Confirm" dg-widget-type="Dialog"]
  [role="heading" level=1 name="Confirm"]
  [role="button" name="OK" testid="dlg-ok"]
  [role="button" name="Cancel" testid="dlg-cancel"]
```

The format is deterministic so prompts are reproducible and the
response cache works.

## Token-budget guards

Before sending, the prompt builder verifies:

- Total estimated tokens ≤ the model's budget for this case.
- If over budget, prune in this order: drop screenshot → drop HTML
  snippet → reduce a11y tree depth → reduce a11y tree breadth.
- Never drop the fingerprint or the negative constraints. Never
  drop the `already_tried` list (the LLM will repeat itself otherwise).

## Why this template shape

- **XML-style tags.** Claude follows these reliably for structure;
  reduces "the LLM ignored part of the message" failures.
- **No prose preamble.** "Hello, please consider this carefully" wastes
  tokens and changes nothing about output quality.
- **Negative constraints inline rather than in system.** Some
  constraints are case-specific (the already-tried list). System stays
  stable and cacheable; per-case detail goes here.
- **Output instructions repeated.** The system prompt sets the contract;
  the user message reinforces it. Helps when the prompt is long and
  attention drifts.

## What changes per call vs per case

Per call (system + user template structure): **stable**. Cacheable on
the Anthropic side via prompt caching.

Per case (variables): **changes**. Constitutes the cache miss
boundary. The case-id hash is also a cache key for our own
deduplication layer in the CLI.
