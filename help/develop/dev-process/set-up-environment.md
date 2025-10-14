---
title: "Environment setup"
sidebar_position: 0
---

This article explains how to set up development environment for developing Datagrok [packages](../develop.md#packages).

## Tools

_NOTE_: To avoid permission issues when installing packages globally via `-g`, use a version manager to install
both `Node.js` and `npm` following
the [instructions](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm#using-a-node-version-manager-to-install-nodejs-and-npm)
.

_NOTE_: On macOS and Unix systems, you may also need to use `sudo` at the beginning of the installation command and
enter the root password if prompted.

1. Install [Node.js](https://nodejs.org/en/)
2. Install [npm](https://www.npmjs.com/get-npm)
3. Install [webpack](https://webpack.js.org/guides/installation/)
4. Install [datagrok-tools](https://www.npmjs.com/package/datagrok-tools)

_NOTE_: The `Node.js` version from [Snap](https://snapcraft.io/)
can produce issues with the datagrok tools installation.
We recommend to avoid using Snap and follow installation instructions provided above.

## Configuration

1. Retrieve your developer key by opening your user profile, clicking on `Developer key` and copying it to the clipboard
2. Configure your environment with the following command:

   ```bash
   grok config
   ```

   Enter developer keys and set the default server. Your credentials will be stored locally in `config.yaml`. Once
   created, this file will be used for [publishing](../develop.md#publishing)
   all your packages. Administrators can manage existing keys and grant or revoke privileges.

## Next steps

Now you are ready to [create your first package](../how-to/packages/create-package.md).

See also:

* [Datagrok JavaScript development](../develop.md)
