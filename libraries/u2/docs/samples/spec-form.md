---
title: A form as a dg-ui/1 spec
group: Specs
order: 10
summary: The UI as JSON — tags, props, binds into the app's signals, a command on click — rendered in one call.
live: yes
---

A spec is JSON: `{"$schema": "dg-ui/1", "root": {...}}`, each node a registered `tag` with `props`, `bind`
entries that point into the render context (`$.name`), and `on` entries that name a command (`cmd:save`) —
never code. `renderSpec(spec, context)` builds it; every bound prop follows its signal in place, two-way for an
input's value, so the app keeps its model in signals and the spec is the view. `manifest.json` lists every tag,
its props and events — what an assistant or a designer writes specs from. Call `registerAll()` once first. The
form's buttons go in a `u2-form-buttons` row (at the form's end unless `align` says otherwise), and a button's
`disabled` binds to why it is off.

## What it takes

{lines} lines, half of them the spec. u2 handles the registry and prop validation (a bad node renders as a
placeholder naming the mistake), the two-way binds, the button row, the disabled Save's reason (its tooltip, also
on keyboard focus), the command routing and disposing the whole tree.
