---
title: Your first control
group: Basics
order: 10
summary: A function component — everything it builds, binds and listens to goes away with it.
live: yes
---

Everything visual in u2 is a control: a `root` element and the scope that owns what it builds. Write one as a
function component: `component((self) => …, {name})` runs the builder with the new control as the owner, so the
heading's live text, the input and its listeners are all released when the control is disposed. State lives in
signals: `person.value` read inside a `computed` is followed, so the heading changes with every keystroke and no
event code.

Hand back `{root, …api}` and the control carries the api — here the input, as `greeting().person`. Name a part for
what it holds: a control's own members (`name`, `root`, `visible`, `dispose`, …) are taken. A piece that is only
part of its parent needs no component of its own: a plain function that builds it, called inside the parent's
builder, is owned by the parent.

## What it takes

{lines} lines. u2 handles the binding of the text to the input, the input's label, focus and accessible name,
the automation name (`data-u2-name="greeting"`) and the cleanup.
