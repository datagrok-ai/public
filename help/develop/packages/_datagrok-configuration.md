<!-- TITLE: Debugging -->
<!-- ORDER: 5 -->

# Datagrok configuration

After you install `datagrok-tools` and run the first command (probably, to [generate a package]), you need to configure
your environment.

## Overview

The Datagrok Tools CLI creates a configuration file under `%USERNAME%/.grok/config.yaml` to store the configurations for
your server, developer keys, and other metadata that's used by Datagrok.

For Windows, the actual path is `C:\Users\%Username%\.grok\config.yaml`, and for Linux &mdash; `~/.grok/config.yaml`.
For more details about the path to `config.yaml`, refer to the [Node.js > OS Homedir section].

You can configure this file from the terminal: just run `grok config` and then enter the needed configurations.

## Developer keys

To develop and publish your packages or applications on Datagrok, you need to obtain the developer keys. A developer key
is an alpha-numeric sequence that you obtain from your user profile on [public.datagrok.ai] or [dev.datagrok.ai].

You can also configure Datagrok locally to run in a Docker container. For that, obtain your unique key from Docker and
add this key to the `./grok/config.yaml` file.

The server specified as _default_ in `config.yaml` must have an up-to-date key.

It's possible to set the key and default server when you run `grok publish`;

```shell
grok publish
```

## See more

* [Datagrok Tools](https://github.com/datagrok-ai/public/tree/master/tools#commands)

## What's next?

* [Packages](#)

[generate a package]: ./_packages.md#creating-a-package

[node.js > OS Homedir section]: https://nodejs.org/api/os.html#os_os_homedir

[public.datagrok.ai]: https://public.datagrok.ai

[dev.datagrok.ai]: https://dev.datagrok.ai

[your user profile on the public server]: https://public.datagrok.ai/u

[your user profile on Datagrok Development server]: https://dev.datagrok.ai/u