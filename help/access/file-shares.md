# File shares

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

Datagrok offers a convenient interface for working with files. You can connect to [popular file systems](supported-connectors.md), including the [Amazon S3 bucket](connectors/s3.md), [Dropbox](connectors/dropbox.md), [Google Drive](connectors/googlecloud.md), and [Git](connectors/git.md), as well as [Windows and Linux network shares](connectors/files.md). In addition, each Datagrok user has a personal directory called **Home** for uploading files using the drag-and-drop feature. Datagrok automatically creates this directory upon signup.

:::note

Connecting to an SMB file storage is only available for on-premise deployment and is not available on the public Datagrok instance (public.datagrok.ai).

:::

:::note developers

You can [create custom connectors](../develop/how-to/access-data.md/#connections).

:::

## Connecting to file storage

### Adding connection

![File share connection parameters](add-a-file-share.gif)

To connect to your file storage, follow these steps:

1. Go to **Data** > **Files**.
1. Open the **New file share** dialog (**Toolbox** > **Actions** > **New file share**). Alternatively, click the **New file share** icon on the **Menu Riboon**.
1. In the the **New file share** dialog, choose the data source from the **Data Source** dropdown. The dialog updates with connection-specific parameters.
1. Set the parameters.
1. Click **TEST** to test the connection, then click **OK** to save it. A pop-up notification displays the connection status.

Some connection parameters have unique characteristics, and it's important to specify them correctly:

* _Directory path_. When connecting to the root directory, leave the **Dir** field empty. Otherwise, enter a directory path.
* _Credentials_. You can specifying credentials manually or using the [Secrets Manager](data-connection-credentials.md), such as the AWS Secrets Manager. When entered manually, Datagrok stores secrets in a [secure privilege management system](../govern/security.md/#credentials). To specify who can change the connection credentials, click the **Gear** icon and select from the **Credential owner** dropdown.

  :::caution

  When connecting to _public buckets_ in AWS S3, always check the **Anonymous** checkbox.

  :::

Once you have established a connection to a folder in your file system, the folder appears in the **File Manager** under the corresponding data source. This connection is referred to as a _file share_. You can view the files and subfolders within the _file share_ by expanding it.

:::note

Like other objects in Datagrok, newly created connections are only visible to the user who created them. To let others access the file share, you must share it (right-click the connection and select **Share...** from the list of options).

:::

### Modifying connection

To modify a connection, right-click it and select **Edit...** from the list of options. To quickly create a connection similar to an existing one, right-click it and select **Clone...**

### File indexing

For file shares, Datagrok supports indexing of folders and [supported file formats](supported-formats.md), including archives such as .tar or .zip.

Connections and folders are indexed by default when you create a connection. File indexing is optional. To index files, select the **Index Files** option when creating a file share.

:::tip

You can enable file indexing at any time. Right-click the file share and select **Edit...** Then, check the **Index file** checkbox in the dialog that appears. Click **OK** to save.

:::

File indexing is a recurring [data job](data-job.md) that runs every hour. Datagrok extracts the following information from the indexed file:

* Filename
* File size, in bytes
* Number of rows and columns
* Column-level information such as name, data type, and [semantic type](../discover/semantic-types.md).

For instance, with indexing, you can browse columns within a CSV file inside a ZIP file:

![File Explorer](./connectors/files-browser.gif "File Explorer")

Indexing helps you find datasets quicker as indexed files appear in the search results based on metadata extracted. For example, you can search for dataframes matching the following criteria across specified or all data providers at once:

* Created in the last month
* Has a column that contains molecules, and
* Has a column named "activity."

:::note

To learn how searching works in Datagrok, see [Smart search](../datagrok/smart-search.md).

:::

## File Manager

The **File Manager** is an interface that allows you to manage connections, browse and preview file content, and perform standard file and folder actions such as opening, downloading, deleting, and renaming. To access an object's context actions, you can right-click it or left-click and expand the **Actions** pane in the **Context Panel** on the left. By clicking a file or folder in the **File Manager**, you can open its preview. Double-clicking a file opens it in Datagrok, and double-clicking a folder expands its content.

:::note

If you don't see a certain action, it may be due to insufficient permissions. For files and folders shared with you, contact the credentials owner. If you are a credentials owner, contact the data source owner.

:::

In addition to the hierarchical browsing, the **File Manager** offers an array of advanced preview and data augmentation capabilities using **Directory**, **Preview**, and **Context Panel**.

The **Directory** section shows the contents of your current folder with three viewing modes: icons, cards, and grid. You can click an object to view its content in the **Preview**, or right-click it to access available context actions. You can also use the search bar located above the **Directory** to search for files and folders within your current directory. The search bar allows you to search for items by name, file extension, or metadata.

The **Preview** is a context-sensitive view that adapts to the the selected object. For folders, the **Preview** generates a [treemap](../visualize/viewers/tree-map.md) that highlights the largest items. For files, the functionality varies based on the file's format and data properties. It includes custom viewers for [supported formats](supported-formats.md), such as interactive spreadsheets for displaying tabular data, cell and image renderers, and chemical and biological structure viewers. You can also view the content of ZIP files and edit Markdown, TXT, and HTML files.

![File browsing and preview](file-manager-file-browsing.gif)

:::note

File preview is limited to files under 10MB. The platform won't display larger files. Unsupported file formats cannot be previewed, but you can download them.

:::

:::note developers
  
You can [add custom formats using package extensions](../develop/how-to/create-package.md). In addition, you can create organization-specific previews:

<details>
<summary> Example: Create custom file viewers </summary>

In this example, a script is executed against the folder content. If the folder contains files that match the file extension parameter PDB, the **Preview** displays a custom NGL viewer to visualize the molecule.

![preview using custom viewer](preview-with-custom-viewer.gif)

To add a custom viewers, you have two options:

* Develop in JavaScript using the [Datagrok JavaScript API](../develop/js-api.md).
* Use the visualizations available for popular programming languages like Python, R, or Julia.

To learn more about each option, see [Develop custom viewer](../develop/how-to/develop-custom-viewer.md).

</details>

<details>
<summary> Example: Create custom folder viewers </summary>

In this example, a [script](../develop/how-to/folder-content-preview.md) is executed against the folder content. If the folder contains files matching the file extension parameter, the **Preview** shows a custom [widget](../visualize/widgets.md) (in this case - the application launch link) every time the folder is opened.

![Suggest an application based on file types](clinical-case-file-manager.gif)

</details>
<details>
<summary> Example: Create custom cell renderers </summary>

In this example, a [script](/develop/how-to/custom-cell-renderers.md) is executed against the [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) strings within the CSV file. The script computes the structure graph and 2D positional data, and renders the structure graphically.

![Smiles renderer](Smiles-renderer.png)

</details>

:::

Lastly, the [Context Panel](../datagrok/navigation.md#properties) provides specific information about a selected object and available actions. It works with folders, files, and data contained within these files. For example, when you click a CSV file, the **Context Panel** updates to show the file's metadata, available context actions, and other relevant information. If you subsequently click any of the dataframe's columns in the **Preview**, the **Context Panel** will update to display information and actions specific to that column, such as summary statistics for the column under **Stats**, or its data and semantic types under **Details**.

![Details on demand](file-manager-details.gif)

:::note developers

**Context Panel** can be extended. You can add custom [info panes](../develop/how-to/add-info-panel.md) and [context actions](../develop/how-to/context-actions.md).

<details>
<summary> Example: Image augmentation </summary>

In this example, a [Python script](/develop/how-to/create-custom-file-viewers.md) creates a custom _info pane_ called **Cell Imaging Segmentation**. This script executes against JPEG and JPG files during the indexing process and extracts custom metadata (such as the number of cells) and performs predefined transformations (such as cell segmentation). When a user selects the corresponding image, the **Context Panel** shows a custom info panel that displays the augmented file preview and the number of detected cell segments.

![Cell image segmentation](Cell-image-segmentation.gif)

</details>

:::

## File sharing and access control

In Datagrok, you can share files in two ways: by sharing the actual file (or folder), or by sharing an URL that points to it. To share an URL, open the file in Datagrok and copy the URL from the address bar. To access the file from the link provided, users must have permissions to open it. Once the file is open, users can download the file and then upload it to their **Home** directory, or [save the file as a project](/datagrok/create-project.md). The URL links never expire and can't be revoked.

:::tip

For tabular formats, you can create dynamic dashboards and share them with others via URL or reference them on external websites. To learn more about dynamic data updates, see [Dynamic data](../datagrok/project.md/#dynamic-data).

:::

Another option is to share directly in Datagrok by creating a _file share_ and specifying access privileges for each shared item, such as separate files and subfolders. Once the item is shared, it appears in the recipient's **File Manager**. When using this method, you have several options:

* Share a _connection_ (root folder) to give access to the entire directory.
* Share a _folder_ to give access to the content of individual folders in your directory.
* Share a _file_ to limit access to individual files within a folder.<!--future-looking-->

To share, follow these steps:

1. Right-click the item you want to share and select **Share...** from its context menu. This action opens the **Share...** dialog.
2. In the identity/email field, start typing a person's name, username, email, or group name, and pick from the list of matching identities.
3. From the respective dropdowns, select access privileges for either or both: (1) the connection and (2) individual
   files/folders. You can select any or all of the following options<!--TBU-->:

    * _Can view_: Users can view, open, and download
    * _Can edit_: Users can rename and edit
    * _Can delete_: Users can delete
    * _Can share_: Users can reshare with any other user or group.<!--how does it work with URL links? -->

    :::caution

    A file's _name_ and _namespace_ are encoded within the URL. When you rename a file (or its location), the link changes accordingly, which may lead to may cause broken URL links, script errors, and similar issues.

    :::

4. Optional. Enter a description in the text field provided. You may also notify the users you share with. If you don’t want to send a notification, uncheck the **Send notification** checkbox.

  :::note
  
  To send an email notification, enter the user's email in the identity/email field. The email notification contains a link to the shared item and entered description. If you enter a user or group name, they will be notified via the Datagrok interface.

  :::

5. Click **OK** to share. Once shared, the shared item appears in the recipient's **File Manager**.

   ![Share a folder](share-the-folder.gif)


You can use the **Sharing** info panel in the **Context Pane** to inspect and quickly adjust access permissions to your _file shares_, send comments to those you're sharing with, and more. The same actions are available from the context menu.

<!--TBD: GIF pending changes in the UI-->

## Data reproducibility

Datagrok provides a visual interface to automate manual, repetitive data ingestion, and data transformation tasks. For more information on workflow automation, see [Data preparation pipeline](data-pipeline.md).

## See also

* [Data access](access.md)
<!--* [Databases](link)
*[Web services](link)
*[Context Pane](link)
*[Indexing](link)-->

## Resources

[![Data Access - File Shares](file-manager-youtube.png)](https://www.youtube.com/watch?v=dKrCk38A1m8&t=417s)
