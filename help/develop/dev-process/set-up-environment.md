---
title: "Environment setup"
sidebar_position: 0
description: How to install Node.js, npm, pnpm, and datagrok-tools, and log in with `grok login` for package development.
keywords:
  - grok login
  - developer key
  - grok config
  - datagrok-tools
  - Node.js
  - npm
  - pnpm
  - local development setup
---

This article explains how to set up a development environment for developing Datagrok [packages](../develop.md#packages).

## Tools

_NOTE_: To avoid permission issues when installing packages globally via `-g`, use a version manager to install
both `Node.js` and `npm` following
the [instructions](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm#using-a-node-version-manager-to-install-nodejs-and-npm)
.

_NOTE_: On macOS and Unix systems, you may also need to use `sudo` at the beginning of the installation command and
enter the root password if prompted.

1. Install [Node.js](https://nodejs.org/en/) 22 or later (npm comes with it)
2. Install [datagrok-tools](https://www.npmjs.com/package/datagrok-tools): `npm install -g datagrok-tools`
3. If you work inside the [public repository](https://github.com/datagrok-ai/public), run `grok setup`
   once at the repository root. It enables [pnpm](https://pnpm.io) and installs dependencies for every
   package. See [Build system](build-system.md). A standalone package created with `grok create` uses
   `npm install` and needs nothing else. The bundler and TypeScript come with the package's
   `@datagrok/build-config` dependency.

_NOTE_: The `Node.js` version from [Snap](https://snapcraft.io/)
can produce issues with the datagrok tools installation.
We recommend avoiding Snap and following the installation instructions provided above.

## Configuration

Log in to your server with a keypair:

```bash
grok login <server>
```

The CLI opens your browser to approve the login and saves the server to `config.yaml`. This file is used for
[publishing](../develop.md#publishing) all your packages. See
[Keypair authentication](../../govern/access-control/keypair-authentication.md).

The developer key is deprecated. As a legacy fallback, copy it from **Developer key...** on your profile, then run
`grok config` and enter it. Datagrok 1.28 and later accept the developer key only from datagrok-tools 6.6.0 or later.

## Next steps

Now you are ready to [create your first package](../how-to/packages/create-package.md).

See also:

* [Datagrok JavaScript development](../develop.md)
