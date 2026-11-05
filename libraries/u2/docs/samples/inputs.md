---
title: Every input, one contract
group: Inputs & validation
order: 10
summary: Text, number, choice, date, slider, switch and text area — each a value with the same options.
live: yes
---

Every u2 input is constructed from one options object and exposes the same surface: `value` (read inside a
`computed` or an effect it is followed; assigning it writes it), `peek()` to read it untracked, `signal` to wire it,
`validity`, `error`, `label`, and the live options `disabled` (the reason the user reads while it is off),
`readOnly`, `required` and `visible`, each a value or a signal it follows. A `Form` lines their captions up in one
right-aligned column and moves to the next field on Enter. The line under the form reads the values — or, when one
is out of range, what is wrong (type 2000 into Dose). Switch Active off: Notes turns off and says why. A date reads
in the browser's zone (`Format.value` writes it in the user's platform format).

## What it takes

{lines} lines. u2 handles the labels and their `for`, the number's range check and units, the date picker and
its keyboard, the slider's value text, the disabled reason, the caption column (it re-aligns on resize) and Enter
navigation.
