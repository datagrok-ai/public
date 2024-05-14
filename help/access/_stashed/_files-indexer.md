---
title: "Indexing File Shares"
---

For data sources that support file access, the Datagrok platform provides a file indexing feature. File indexer is a
data job with hourly recurrence that indexes files with
[supported](../files/supported-formats.md) extensions in the storage specified by data connection.

# Indexing file shares

>Key concept: _File share_
>
>We use the term _file share_ to refer to a _file storage connection_ or any folder or file shared with others. To learn more about file shares, see [File storage](../files/files.md).

For _file storage_<!--add link when ready--> Datagrok supports indexing of folders and _supported file formats_<!--add link when ready--> (including archives such as .tar or .zip).

Connections and folders are default indexed when you create a connection<!--add link when ready-->. File indexing is optional. To index files, select the **Index Files** option when creating a _file share_<!--add link when ready-->.

>Tip:  You can turn on file indexing at any time. Right-click the _file share_ and select **Edit...** Then, check the **Index file** checkbox in the dialog that appears. Click **OK** to save.

File indexing is a recurring _data job_<!--add link when ready--> that runs every hour. Datagrok extracts the following information from the indexed file:

* Filename
* File size, in bytes
* Number of rows and columns
* Column-level information such as name, data type, and [semantic type](../govern/catalog/semantic-types.md).

For example, thanks to indexing, you can browse columns within a CSV file residing within a ZIP file.

![File Explorer](./connectors/files-browser.gif "File Explorer")

Indexing helps you find datasets quicker as indexed files appear in the search results based on _metadata_ extracted. For example, you can search for dataframes matching the following criteria across specified or all data providers at once:

* Created in the last month
* Has a column that contains molecules, and
* Has a column named "activity."

>Tip: To learn how searching works in Datagrok, see [Smart search](../datagrok/smart-search.md).

See also:

* [File shares](../files/files.md)
* [Data connection](../access.md#data-connection)
* [Data job](_data-job.md)
