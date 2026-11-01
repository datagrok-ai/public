---
title: Fold what matters less
group: Panels & layout
order: 20
summary: Accordion panes built on first open, each with a hint of what it holds; the user's folds remembered.
live: yes
---

An `Accordion` shows what matters and folds the rest. A pane's content can be a function, built the first time it
opens; its `hint` says what it holds while folded (and follows a signal); a `validity` signal marks it — a problem
is never hidden by a fold; `persistKey` remembers whether the user opened or closed it. `flush: true` lines the
panes up as rows of a form.

## What it takes

{lines} lines. u2 handles the disclosure semantics (`aria-expanded`, Enter/Space), lazy building under the pane's
own scope, the hint and the mark, and the remembered state.
