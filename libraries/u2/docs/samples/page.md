---
title: A view's page
group: Panels & layout
order: 35
summary: Two titled sections over a pinned action row: the body alone scrolls, Create stays in view and says why it is off.
live: yes
height: 360
---

`Page` lays a view's content out: an optional `header` that stays put, a body that alone scrolls, and a `footer`
pinned under it — the action row, its buttons on the right and a `statusLine` on their left saying what Create
makes, or why it cannot yet. Each section is a heading with its `actions` (here, a Clear icon) over its content,
and the forms on the page share one label column, so Name and Description start their editors at one place.
`maxWidth` caps the body; its forms go as wide. Spacing comes from the token scale; a package tunes it with the
`--u2-page-*` properties on an element of its own (README › Styling), never with a rule over u2's classes.

## What it takes

{lines} lines and no CSS: no flex column with a scrolling middle, no pinned footer, no `max-width` rule by hand.
