<!-- TITLE: Develop custom file handlers -->

# File handlers

To handle custom file formats, register a function with the `file-handler` tag, and specify the comma-separated
extensions in the `meta.ext` parameter. Function's input is either a string or a list of bytes, the output is list of
[tables](../../overview/table.md).

For example, the following function will get executed whenever a user opens a file with the "fasta"
extension:

```javascript
//input: string content
//output: list tables
//tags: file-handler
//meta.ext: fasta
function fastaFileHandler(content) {
    // ... processing files ...
    return tables;
}
```
