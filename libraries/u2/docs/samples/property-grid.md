---
title: A settings grid from descriptors
group: Forms
order: 30
summary: Typed properties in categories, edited in place; the values one record you can bind to.
live: yes
---

`PropertyGrid` is for many small typed settings — a viewer's look, a tool's options. Describe each property
(`string`, `int`, `double`, `bool`, `choice`, `string_list`, with `min`/`max`, `choices` and a `category`) and
hand it the values; it picks the editors, groups them under collapsible categories, and replaces its `values`
record on every edit, so anything can follow it — here, the line under the grid. An object that carries platform
metadata (`DG.Property`) gets the same through `propertyForm` (see Inside Datagrok).

## What it takes

{lines} lines and {cssLines} of CSS (the line under the grid). u2 handles the editors per type, the ranges, the
categories, the descriptions as tooltips, and the change notification (`onChanged` gets the property and its
value).
