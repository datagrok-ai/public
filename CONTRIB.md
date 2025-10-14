# Contributor's guide

## Requirements

1. Node version [18.12.x](https://nodejs.org/dist/v18.12.0/)
2. Npm version 9.x.x: `npm install -g npm@9.x.x`
3. Latest `typescript`.

We are only using pure JavaScript in the packages not yet converted to TypeScript, such as
`public/packages/Charts`.

We use **only** TypeScript in all actual and new packages. Avoid using pure JavaScript.

## Code style

In this public repo, we follow the [Google JavaScript Style Guide](https://google.github.io/styleguide/jsguide.html).

In particular:

1. **Add** intermediate spaces in `a, b`, `1, 2, 3`, `if (...`, `...) {`, `a == b`, and similar

2. **Do not** add empty lines in the beginning and/or the ending of the code block

3. The default is **2 spaces per tab**

4. We are *not strict*
   about [braces for code blocks of a single statement](https://google.github.io/styleguide/jsguide.html#formatting-braces-all)
   and [trailing commas](https://google.github.io/styleguide/jsguide.html#features-arrays-trailing-comma)
   .

5. When creating a package, use the `--eslint` flag to get an
   up-to-date [configuration file](https://github.com/datagrok-ai/public/blob/master/tools/package-template/.eslintrc.json)
   .

6. **Do not** delete `package-lock.json` from the repository. Update it when needed

7. Document your code when there is a need for it, but do not overdo it. For instance, there is no reason to include
   information that is already in the function/class signature, such as types of parameters. Often, a one-liner is
   enough.

If you are using WebStorm IDE, we recommend you to stick to its defaults for JS/TS formatting, except for the spaces
settings: change its default value of 4 to 2.

Thank you for following the style!

## Using a linter

If you have created a package with `grok create ... --eslint`, the `package.json` file will already have the
required `eslint` dependencies, which would be installed together with others once you call `npm install`. Also, the
file `.eslintrc.json` with all necessary settings will be pre-created.

However, if you work with packages either not created with `grok create` or not having all the conditions above met, you
should still set up `eslint`. It is straightforward:

* Install `eslint`: call `npm install eslint --save-dev -g`
* Install `eslint` `google` settings: `npm install eslint-config-google -g`.

Make sure that your `.eslintrc.json` is actualized to using TypeScript (look
for `"parser": "@typescript-eslint/parser"`). If that's not the case, populate your `.eslintrc.json`
file with the settings
[matching this file from `datagrok-tools`](https://github.com/datagrok-ai/public/blob/master/tools/package-template/.eslintrc.json).

## Git

Check our [git recommendations](https://datagrok.ai/help/develop/dev-process/git-policy) to work with the repository