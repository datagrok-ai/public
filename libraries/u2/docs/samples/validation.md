---
title: Validation and Save
group: Inputs & validation
order: 20
summary: Required fields, a validator, a cross-field rule — Save disabled with the first problem as its reason.
live: yes
---

An input has two verdicts. `validity` judges the value from the start — what Save goes by; `error` is what the
input shows, only once it was edited, left or validated — so a fresh form is not covered in red. A `Form`
aggregates both, and `form.problem` names the field: "Email: Not an email address". `required` is checked first:
an empty value says "Value can't be empty" and never reaches the validators, so the email check need not pass on
nothing (blank text is empty too). A validator that reads another input — From before To — re-judges whenever that
input changes.

## What it takes

{lines} lines. u2 handles the touched state, the red underline and message under each field, `aria-invalid` and
`aria-describedby`, the hanging required star, the aggregate verdict and the reason on the disabled button.
