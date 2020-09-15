# datagrok-tools

Utility to upload and publish packages to Datagrok.

## Installation

```
npm install datagrok-tools -g
```

## Usage

1. Create a configuration file common to all packages:

```
grok config
```

Enter developer keys and set the default server.

2. Get a package template for your current folder:

```
grok create
```

Alternatively, pass the package name:

```
grok create MyPackage
```

A new folder `MyPackage` will be created automatically as well as its contents.

3. To upload a package, use the following command:

```
grok publish
```

The default mode is set to `--debug --build`. You can easily change it with options `--release` and `--rebuild` respectively.

In addition, you can pass another server either as URL or server alias from the `config.yaml` file:

```
grok publish dev
grok publish https://dev.datagrok.ai/api -k key-for-dev
```

The developer key, if specified, gets updated in the `config.yaml` file. All new servers along with their keys will also be added to this file.

4. When having a package under development, please run

```
grok migrate
```

This command will convert your scripts in `package.json` and copy keys from `upload.keys.json` to `config.yaml`.

Run `grok` for instructions and `grok <command> --help` to get help on a particular command.
