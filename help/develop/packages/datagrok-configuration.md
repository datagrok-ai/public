# Datagrok configuration

After you install `datagrok-tools` and run the first command (probably, to [generate a package]), you need to configure
your environment.

## Table of Contents

* [Overview](#overview)
* [Developer Keys](#developer-keys)

## Overview

The Grok CLI creates a configuration file under `%USERNAME%/.grok/config.yaml` to store the configurations for your
server, developer keys, and other metadata that's used by Datagrok.

You can configure this file from the terminal by running `grok config` and then entering the needed configurations. 
Alternatively, you can change this file manually.

## Developer Keys

To develop your packages or applications on Datagrok, you need to obtain the developer keys.

## How to obtain the public developer key

To obtain the public developer key, open [your user profile on the public server] and then click **Developer key**.
Datagrok administrators can manage existing keys and grant or revoke privileges.

## How to obtain the development developer key

To obtain the developer key for development server, open [your user profile on the development server] and then click
**Developer key**.

## How to obtain the local developer key

You can configure Datagrok locally to run in a Docker container. For that, obtain your unique key from Docker and then
add this key to the `./grok/config.yaml` file.

[your user profile on the public server]: https://public.datagrok.ai/u
[your user profile on Datagrok Development server]: https://dev.datagrok.ai/u