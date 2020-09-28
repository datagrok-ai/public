# datagrok-tools

Utility to upload and publish [packages](https://datagrok.ai/help/develop/develop#packages) to Datagrok.

## Installation

```
npm install datagrok-tools -g
```

## Usage

1. Configure your environment with the following command:  
    ```
    grok config
    ```
    Enter developer keys and set the default server. The developer key can be retrieved from your user profile (for example, see https://public.datagrok.ai/u).
2. Create a new package by running this command:
    ```
    grok create MyPackage
    ```
    A new folder `MyPackage` will be created automatically as well as its contents.
3. Once you have completed the work on your package, upload it by running:
    ```
    grok publish
    ```

If you are developing a package using the old template, please run `grok migrate`. This command will convert your scripts in `package.json` and copy keys from `upload.keys.json` to `config.yaml`.

Run `grok` for instructions and `grok <command> --help` to get help on a particular command.

Read more about package development in [Datagrok's documentation](https://datagrok.ai/help/develop/develop).
