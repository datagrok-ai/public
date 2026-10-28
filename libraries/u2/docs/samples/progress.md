---
title: A long job with progress
group: Feedback
order: 20
summary: A progress bar with a live description, Run and Cancel that say why they are off, a timer owned by the job.
live: yes
---

`ProgressBar` shows a 0..1 `value` (or runs indeterminate) with a live description; `progress.value = x` moves
it. The timer is started on the component's own scope (`self.setInterval`): it stops by itself when the job is
disposed mid-run — there is no `clearInterval` to remember. Run and Cancel are disabled with their reason, so the
user always knows why.

## What it takes

{lines} lines. u2 handles `role=progressbar` with its values, the percent text, the live subtitle, the disabled
reasons and stopping the timer with its owner.
