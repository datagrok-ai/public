<!-- TITLE: Define semantic type detectors -->

# Defining semantic type detectors

Datagrok helps you get the most out of your data by encoding its meaning into
[semantic types](../../discover/semantic-types.md). In addition to the types that
are [automatically detected](../../discover/semantic-types.md#automatic-semantic-type-detection)
by the platform, you can define your own semantic types. This can serve multiple purposes: from custom rendering and
viewers to info panels and predictive models.

## Function annotation

Semantic types detectors have a tag `semTypeDetector`, take a single `column`
as input, and return either a `string` containing the semantic type or `null`. Their names typically start with
the `detect` prefix, e.g., `detectNucleotides`
or `detectRDSmiles`. These functions should live in `detectors.js`, a special file inside
your [package](../develop.md#packages). There you simply define a class named `<package_name>PackageDetectors` that
subclasses `DG.Package`
and add one or more detectors:

```javascript
class SequencePackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.name.startsWith('nuc')) {
      col.semType = 'nucleotides';
      return col.semType;
    }
    return null;
  }
}
```

Notice how the column properties are used. Such properties as `col.name` help determine a semantic type,
but `col.semType` actually assigns the detected type to the column.

To get a template for a detector function, run this command from your package directory:

```shell
grok add detector <semantic-type-name>
```

A column can have only one semantic type. In cases when several detectors match a column, its semantic type depends on
the order in which the detectors were triggered. Standard platform detectors do not necessarily get a preference.

## Advanced checks

Semantic type detectors are uploaded separately from the rest of a package
(therefore, you might not see webpack warnings about syntax errors in
`detectors.js`). Datagrok calls detectors each time the user opens a table. To quickly inspect the data, these functions
have to be lightweight, simple, and efficient. Often checking the data type (`col.type`) and applying regular
expressions to the column name (`col.name`) will be enough for such tests. However, there are ways to carry out a more
advanced check. For instance, the following example makes use of column statistics:

```javascript
class ViewersPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMagnitude(col) {
    if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) &&
      (0 < col.min && col.max < 10) && col.name.toLowerCase() === 'magnitude') {
      col.semType = 'Magnitude';
      return col.semType;
    }
    return null;
  }
}
```

Apart from `col.min` and `col.max`, there are other
[descriptive statistics](https://public.datagrok.ai/js/samples/data-frame/stats)
calculated for a column.

If the columns you are working with contain a lot of unique categories, use a special sampling
method `DG.Detector.sampleCategories()` to run your checks on a random subset of column values (see an example in the
public package
[NglViewer](https://github.com/datagrok-ai/public/blob/master/packages/NglViewer/detectors.js)).

Empty values deserve special attention: make sure you don't match them with your semantic type, for example, empty
strings should not match your regular expressions. Otherwise, this type might be assigned to a column consisting solely
of nulls, which won't add any helpful insights to the data profile.

![Detected Types: Latitude, Longitude, Magnitude](semantic-type-detectors.gif "Detected Types: Latitude, Longitude, Magnitude")

See also:

* [JavaScript Development](../develop.md)
* [Semantic types](../../discover/semantic-types.md)
* [JavaScript API Samples: Override standard semantic types](https://public.datagrok.ai/js/samples/data-frame/advanced/semantic-type-detection)
* [JavaScript API Samples: Column statistics](https://public.datagrok.ai/js/samples/data-frame/stats)
* [How to add an info panel](add-info-panel.md)
