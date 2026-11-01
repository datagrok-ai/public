---
title: A spec that reacts with no code
group: Specs
order: 20
summary: Nodes bound to each other — a heading follows a name, units follow a choice, a switch shows a field.
live: yes
---

A named node is a bind source: `$.name` is that input's value, so nodes can follow each other with no code at
all. Every declared prop is live — text, a postfix, `visible` — and follows its source in place: the node is
never rebuilt, so focus, caret and typed values survive. (A prop marked structural, such as an input's `inline`,
is read once; the console says so.) This is the form an assistant or the visual designer produces, and
`dump()` writes it back with the current values folded in.

## What it takes

{lines} lines, nearly all of them the spec itself. u2 handles resolving the names (a node may bind to one declared
later), the live props, the two-way value binds and disposing the tree.
