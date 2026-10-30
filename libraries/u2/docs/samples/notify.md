---
title: Notifications
group: Feedback
order: 10
summary: Info, warning and error balloons; one kept to a single copy; one with details to copy.
live: yes
height: 300
---

`notify.info / warning / error` show a balloon in the top-right stack — inside Datagrok, the platform's own
stack. An info or a warning hides by itself after 10 s (hovering pauses it); an error stays until the user closes
it. Each closes on a click and copies its text — plus `copyText`, the details a bug report needs — from its copy
button. `singleKey` keeps one balloon per key on screen (press Warning twice); `oneTimeKey` shows a tip once ever. A
balloon is how a failure reaches the user: an action's `run` that throws or rejects is only logged to the console,
so say what failed with `notify.error`.

## What it takes

{lines} lines. u2 handles the stack, the timers and the pause on hover, the copy and close buttons, the
de-duplication, `role=status` / `role=alert` for screen readers, and routing into the platform's balloons.
