---
title: "Custom package settings editors"
---

Datagrok provides a convenient way to define [package settings](../develop.md#package-settings). You only need to define
the properties in the package.json file, and the platform takes care of the rest: the UI gets automatically generated,
and settings could be edited either individually or in a centralized manner.

However, sometimes you need to provide a custom UI for the settings editor. To do it, define a function that returns
a `widget` and is tagged as `packageSettingsEditor`. This is all!
Now, when a user clicks on the package and expands the "Settings" pane on the right, our custom UI gets shown.

```javascript
//output: widget kpi
//tags: packageSettingsEditor
export function powerPackSettingsEditor(): DG.Widget {
  return new PowerPackSettingsEditor();
}

export class PowerPackSettingsEditor extends DG.Widget {
  constructor() {
    super(ui.divText('I am a custom package settings editor'));
  }
}
```

* [Package settings](../develop.md#package-settings)
