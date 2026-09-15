# Contributor's guide

## Requirements

1. Node 22 or later.
2. pnpm, enabled with `corepack enable` (the version is pinned in the root `package.json`). The
   repository is one pnpm workspace: run `grok setup` once (it enables pnpm, installs, and removes npm-era leftovers); never `npm install` inside a
   package.
3. `datagrok-tools`: `npm install -g datagrok-tools` (inside the repository, package scripts use the
   workspace copy automatically).

TypeScript, the bundler and eslint come with the workspace (`@datagrok/build-config`); nothing else is
installed globally. See [packages/BUILD.MD](packages/BUILD.MD) for the build commands.

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

5. The eslint configuration is one file at the repository root (`.eslintrc.json`); packages carry none.
   `pnpm run lint` in a package (or `pnpm turbo run lint`) applies it.

6. `pnpm-lock.yaml` at the repository root is the only lockfile. Commit it whenever `pnpm install`
   changes it; never add a `package-lock.json`.

7. Document your code when there is a need for it, but do not overdo it. For instance, there is no reason to include
   information that is already in the function/class signature, such as types of parameters. Often, a one-liner is
   enough.

If you are using WebStorm IDE, we recommend you to stick to its defaults for JS/TS formatting, except for the spaces
settings: change its default value of 4 to 2.

Thank you for following the style!

## Using a linter

Inside the repository, eslint and the TypeScript parser are workspace devDependencies and the one
configuration lives in `.eslintrc.json` at the root: run `pnpm run lint` in a package, or
`pnpm turbo run lint` for everything. A standalone package created with `grok create` gets the same
rule set through `@datagrok/build-config`; run `npm run lint` there.

## Git

Check our [git recommendations](https://datagrok.ai/help/develop/dev-process/git-policy) to work with the repository