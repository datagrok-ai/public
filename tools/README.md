# datagrok-tools

Utility to upload and publish packages to Datagrok.

## Installation

```
npm install datagrok-tools -g
```

## Usage

### datagrok

1. Create a folder:

```
mkdir MyPackage
cd MyPackage
```

2. Get a package template:

```
datagrok-init
```

Enter your package name (by default, it is the name of the current folder), the remote server's URI, and your developer key.

3. Install dependencies:

```
npm install
```

The only package listed in `package.json` is `datagrok-api`. You can add new dependencies later.

### grok

Discover our new package management tool, all while it's under development.

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

Run `grok` for instructions and `grok <command> --help` to get help on a particular command.
