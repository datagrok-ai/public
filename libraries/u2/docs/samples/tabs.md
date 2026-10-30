---
title: Tabs that mark their problems
group: Panels & layout
order: 10
summary: Each tab a form; a tab with a failing field is marked once the user tried to save.
live: yes
---

`TabStrip` takes its tabs as data: an id, a `title`, the content (a function is built on first open). A tab whose
content is a `Form` marks itself from the form's shown error — nothing on a fresh form, the first failing field
once Save validated it — and any other tab takes a `validity` signal. Save here validates every form, so the
marks appear on the tabs you have not opened yet.

## What it takes

{lines} lines. u2 handles the tab semantics (`role=tablist`, arrows, Home/End), lazy building, the marks with the
message as their tooltip, and keeping each tab's state while another is shown.
