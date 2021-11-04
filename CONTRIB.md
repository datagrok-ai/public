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

## Using a linter

If you have created a package with `grok create ... --eslint`, the `package.json` file
will already have the required `eslint` dependencies which would be installed together
with others once you call `npm install`. Also the file `.eslintrc.json` with all necessary
settings will be precreated.

However, if you work with packages either not created with `grok create` or not having
all the conditions above met, you should still set up `eslint`. It is straightforward:

* Install `eslint`: call `npm install eslint --save-dev -g`
* Install `eslint` `google` settings: `npm install eslint-config-google -g`.

Make sure that your `.eslintrc.json` is actualized to using TypeScript (look for `"parser": "@typescript-eslint/parser"`). If that's not the case, populate your `.eslintrc.json` file
with the settings
[matching this file from `datagrok-tools`]().