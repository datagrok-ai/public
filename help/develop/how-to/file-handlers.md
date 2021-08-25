<!-- TITLE: Develop Custom File Handlers -->

## File Handlers

To handle custom file formats, simply register a function with
the "file-handler-<extension>" tag (you can specify more than one).
Function's input is either a string or a list of bytes, the output is list of
[tables](../overview/table.md).

For example, the following function will get executed whenever a user
opens a file with the "fasta" extension:

```js
//input: string content
//output: list tables
//tags: file-handler
//meta.ext: fasta
function fastaFileHandler(content) {
    // ... processing files ...
    return tables;
}
```