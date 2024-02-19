---
title: FAQs
sidebar_position: 11
---

*Question:*
How do I set up VS Code or WebStorm for debugging on Datagrok? Debugging doesn't work for me now.

*Answer:*
Check out our guides on [VS Code] and [WebStorm]. If you configured VS Code or WebStorm following the guides, but
debugging still doesn't work, ask our technical gurus on the [Community Forum].

---

*Question:*
I've installed `npm` and `datagrok-tools` packages, but `grok` tools aren't available on my paths.

*Answer:*
It's likely that the `datagrok-tools` library wasn't installed globally. Install `datagrok-tools`
globally:

```sh
# First, remove the locally installed package...
npm uninstall datagrok-tools

# And then install datagrok-tools globally using the -g flag
npm install -g datagrok-tools
```

---

*Question:*
Others don't see my published package on the selected Datagrok instance. How can I fix this?

*Answer:*
If you publish a package in debug mode using the standard `grok publish` command, only you will see the published
package. Run `grok publish --release` so others could see your package, too.

---

*Question:*
What is the best approach to synchronize custom written filters between multiple
filter viewers?

*Answer:*
We currently use a combination of events for synchronization and saving state
while filtering to dataFrame to initiate new viewers from this state:

* `d4-filter-criteria-changed` to notify other filters.
* `dataFrame.rows.filterStates` to store the filter's state - this array is
emptied before each call of `onRowsFiltering`, and filters add their states
while filtering.

Depending on where you use the filtering, you have two options for saving the
filter state:

* For filtering within native grok places (the **Filter Panel** and **Context
Pane**), save the filter state via the `RowList_AddFilterState` method. In
this case Datagrok automatically finds and applies the state.
* For filtering in other places, you can store the filter state in dataFrame
tag. In this case you need to add the code for finding the saved state.

---
