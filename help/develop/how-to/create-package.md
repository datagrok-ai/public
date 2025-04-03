---
title: "Creating a package"
---

This documents explains how to create a simple package that contains an
[info panel](add-info-panel.md) displaying some simple statistics for a text in a selected grid cell.

## 1. Create a package

Before we start, [set up development environment](../dev-process/set-up-environment.md) if you haven't done it already.

Now, let's create a package:

```shell
grok create TextStats --test
```

This creates a folder `TextStats` with the default [package structure](../develop.md#package-structure). For the newly
created package, you have to install the dependencies. Run this from the `TextStats` folder:

```shell
npm install
```

:::note

The internal name of the package should follow **PascalCase** (e.g., `TextStats`), while the friendly name should use **Title Case** (e.g., `Text Stats`). Additionally, all publicly available plugins must have custom icons without exception.

:::

## 2. Implement the panel function

Add the actual panel's code at `TextStats/src/package.ts`:

```typescript
//name: TextStats
//tags: panel, widgets
//input: string str {semType: text}
//output: widget result
export function textStats(str: string) {
  // for 'gattaca', produces {"g": 1, "a": 3, "t": 2, "c": 1}
  const symbolCounts = Array.from(str).reduce((counts, ch) => {
    counts[ch] = (counts[ch] || 0) + 1;
    return counts;
  }, Object.create(null));
  return new DG.Widget(ui.divV([
    ui.divText("Counting characters:"),
    ui.divText(`${JSON.stringify(symbolCounts)}`)
  ]));
}
```

That's it! What creates a panel is this function plus the annotation contained in the comments preceding it. Datagrok
recognizes these comments and makes the function `textStats` become a `panel`
producing a `widget`, taking a `string` as an input.

## 3. Build

Run webpack from the TextStats folder:

```shell
npm run build
```

Note. The `build` script in `package.json` includes commands to process your package, e.g., `webpack`. You are not
supposed to modify this script.

## 4. Publish

Now, let's [publish](../develop.md#publishing) our package to the `dev` server:

```shell
grok publish dev
```

The return code should be `0` to indicate a successful deployment.

## 5. Test

Now go to Datagrok and open any data file with string columns. This could be a `demog` dataset from our demo datasets.
Find any text column, right-click on it, pick `Column Properties`, then add a tag: `quality : text`. This tag tells
Datagrok that there is a plain text stored inside this column. Navigate to a text cell and find your freshly added panel
with text stats on the right side of Datagrok UI. Note, that semantic type `text` of a column matches `{semType: text}`
in function parameter declaration. Datagrok also can detect semantic types automatically, or with help
of [Detectors](./define-semantic-type-detectors.md) .

We have completed the deployment. Let's navigate to `Manage | Packages` and find your package in the list. Note it has
your name prefixed with a `v.`, which means it's only published for you. This is called a Debug mode. To make it
available to the user or a group of interest, you can `Share` it to the group via right-click menu on the package. Don't
forget to publish the package as `grok publish public --release`: this now makes this package _released_
for these groups of interest.

See also:

* [Datagrok JavaScript development](../develop.md)
* [Test packages](test-packages.md)
