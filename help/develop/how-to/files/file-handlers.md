---
title: "Develop custom file handlers"
---

To handle custom file formats, register a function with the `fileHandler` role, and specify the comma-separated
extensions in the `meta.ext` parameter. Function's input is either a string or a list of bytes, the output is list of
[tables](../../../datagrok/concepts/table.md).

For example, the following function will get executed whenever a user opens a file with the "fasta"
extension:

```javascript
//input: string content
//output: list tables
//meta.role: fileHandler
//meta.ext: fasta
function fastaFileHandler(content) {
    // ... processing files ...
    return tables;
}
```
