---
title: An app's home header
group: Panels & layout
order: 40
summary: The icon, what the app does as bullets (text or links), and Learn more — the same in every app.
live: yes
---

An app's home view opens with the same header across Datagrok: the app's big icon beside what it does, one
bullet each — a bullet is text, an element, or a link that runs your code — and Learn more under them. With a
`title`, the name and a subtitle stand beside a smaller icon; `bottomLine` rules it off from what follows. The
icon is an image URL, a Font Awesome name, or an element. What the user starts from (a catalog, a form) goes
right under it. The header follows its own width, so build it under an owner: the view's host in an app (see "An
app home with a catalog"), a component here.

## What it takes

{lines} lines. u2 handles the layout and its narrow-pane response (under a `u2-adaptive` host), the links and
their keyboard, the external-link marker and the header's state for tests.
