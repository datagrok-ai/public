---
title: Elements in a card
group: Basics
order: 20
summary: Badges, tags, a timestamp and a link, composed in a clickable Card with an action.
live: yes
---

The element factories (`badge`, `tag`, `timestamp`, `link`, `span`, `divH`…) return plain elements styled by
u2; `Card` lays them out as a header (icon, title, subtitle, actions), a body and a footer. A click anywhere on
the card opens it, while a click on its own controls — the star, the link — stays theirs. `timestamp` shows the
short date in the browser's zone and the full date and time on hover. A function like `compoundCard` that builds
one piece for whoever calls it needs no component of its own: the card is owned by the builder it is called in.

## What it takes

{lines} lines, no CSS: `divH(children, {gap: 's'})` spaces the tags on the token scale. u2 handles the card's layout,
its keyboard activation (Enter, Space) and focus ring, the star's pressed state (`aria-pressed`) and its tooltip,
and the link opening in a new tab.
