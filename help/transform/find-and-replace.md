---
title: "Find and replace"
format: mdx
---

```mdx-code-block
import FindAndReplace from './find-and-replace-dialog.png';
```



To find and replace values in tables, press <kbd>Ctrl + H</kbd>. This opens the
**Find and replace** dialog. 

<img src={FindAndReplace} width="235"/>

This feature is similar to other text editors but with additional
options:
* Choose specific columns (current, selected, all) and rows (selected, filtered)
  to target for the operation.
* Automatically filter all matching rows.
* Use [search patterns](../explore/search-filter-select/data-search-patterns.md) for matching non-textual columns. For
  example, you can search for `this week` in a datetime column and replace it
  with a specific date like `Nov 7, 2000` or `7/11/2000`.
* In string columns, you can also match case, entire word, or use regular
  expressions for more complex searches.

All replace operations are logged in the [Console](../datagrok/navigation/panels/panels.md#console).