# Contributor's guide

## Code style

In this public repo, we follow the [Google JavaScript Style Guide](https://google.github.io/styleguide/jsguide.html).

In particular:

1. **Add** intermediate spaces in `a, b`, `1, 2, 3`, `if (...`, `...) {`, `a == b`, and similar

2. **Do not** add empty lines in the beginning and/or the ending of the code block.

3. The default is **2 spaces per tab**.

4. We are *not strict* about [braces for code blocks of a single statement](https://google.github.io/styleguide/jsguide.html#formatting-braces-all) and [trailing commas](https://google.github.io/styleguide/jsguide.html#features-arrays-trailing-comma).

5. When creating a package, use the `--eslint` flag to get an up-to-date [configuration file](https://github.com/datagrok-ai/public/blob/master/tools/package-template/.eslintrc.json).

If you are using WebStorm IDE, we recommend you to stick to its defaults for JS/TS formatting,
except for the spaces settings: change its default value of 4 to 2.

Thank you for following the style!