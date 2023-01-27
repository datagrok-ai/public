---
title: File shares. File Manager
---

**File Manager** lets you work with files and directories on your system from the convenience of your web browser. You
can browse, preview, open, create, delete, rename, download, clone, and share files and directories.

To open **File Manager**, from the main menu, select **Data**>**Files**.

In this article:

* [Getting data](#getting-data)
* [Sharing files](#sharing-files)
* [Manage file shares](#manage-file-shares)
* [Manage files](#manage-file-shares)
* [Preview file and folders](#browse-and-preview)
* [Customizations](#customizations)
* [Automation](#automation)

> Key concept: _file share_
>
> We use the term _file share_ when referring to the following:
>
>* Connections to folders and files located on a Datagrok server that you access with a Datagrok client.
>* Any folder or file in the **File Manager** that has been shared with others.
>
>_File shares_ are entities<!--inseert a link-->, meaning you can perform a standard set of operations against them (for
> example,
> annotate, set access privileges, [use in automation workflows](data-pipeline.md), or enable discovery by other
> Datagrok
> users).
>
> Note: For enterprise users, the Datagrok administrator defines which file shares you can access and which privileges
> you have in them.

## Getting data

Each Datagrok user has a personal directory for files and folders they want to upload directly to Datagrok. This
personal directory is accessible in the **Folder Tree** under the name _Home_. The **Home** directory is created
automatically upon signup. By default, files and folders in your **Home** directory are visible to you only.

In addition, you can connect network folders hosted on supported data sources, and access files and directories shared
with you.

### Supported data sources

You can connect to files and folders located on a remote server or a mapped drive. Out-of-the-box, Datagrok provides
connectors to the following data sources:

* [An Amazon Simple Storage Service (S3) bucket](connectors/s3.md)

* [Windows and Linux network shares](connectors/files.md)

  > Note: Connecting to an SMB file share (that is, a file share mounted on a server, such as Linux or Windows network
  > shares) is only available for on-premise deployment and is not available on the public Datagrok instance
  > (public.datagrok.ai).

* [Dropbox](connectors/dropbox.md)

* [Git](connectors/git.md)

* [Google Drive](connectors/googlecloud.md)

> Developers: You can [create custom connectors](../develop/how-to/access-data.md/#connections).

### Add new connection

To add a connection, follow these steps:

1. From the main menu on the left, click **Data** > **Files**.
1. Open the **New file share** dialog by either: (1) expanding the **Actions** panel and clicking **New file share**, or
   (2) Clicking the **New file share** icon on the **Menu Ribbon**.
1. From the **Data Source** dropdown list, select the desired data source. This action updates the dialog with
   connection-specific parameter fields.

1. Fill in all dialog fields displayed.

   > Notes:
   >
   >When connecting to the root directory, leave the **Dir** field empty. Otherwise, enter a directory path within the
   > file share.
   >
   >You can enter _credentials_ (typically, login/password) manually. When entered manually, Datagrok stores secrets in
   > a secure [privilege management system](../govern/security.md#credentials). You can also connect using Datagrok's
   > integration with the AWS Secrets Manager (see [Secrets Managers](data-connection-credentials.md/#secrets-managers))
   .
   >
   >You can define connection credentials for each user or group. To do so, select a user or group from the **Credential
   > owner** dropdown, then enter the appropriate credentials in the fields provided. IMPORTANT: When connecting to
   > _public buckets_ in AWS S3, always check the **Anonymous** checkbox.

1. Click **TEST** to the connection, then click **OK** to save it.

   ![File share connection parameters](add-a-file-share.gif)

Once you connect the network folder, it becomes the root folder for this connection. You can expand and view its content
by double-clicking it.

> Note: When you have a connection set up precisely the way you want it, you can clone it and make additional changes as
> needed:
>
> 1. Right-click the connection and select **Clone...**. This action opens the **Edit Connection** dialog.<!--confusing
>    UI, discuss-->
> 1. Type in a **Name** for the new connection and make other changes as needed.
> 1. Re-enter password or access keys.
> 1. Click **OK** to save the new connection.

### Modify a connection

You can modify a connection at any time. To modify a connection:

1. Right-click the connection to open its context menu, then click **Edit...**. This action opens the **Edit
   connection** dialog.
1. In the **Edit Connection** dialog, change the connection name, parameters, and credentials as needed.
1. Click **TEST** to the connection, then click **OK** to save the changes.

## Sharing files

### Sharing methods

You have two sharing methods:

1. **Sharing files** with individual users and groups directly from the **File Manager**. When you use this method, you
   can specify access privileges for each shared item (such as separate files and subfolders). Once the item is shared,
   it appears in the recipient's **File Manager**.
1. **Distributing a file link**. Datagrok uses unique URLs for each file and folder. You can use these URLs as a quick
   way to point users to a file or reference the file on external websites.

   Unlike sharing _files_, when you share a _file link_, the shared file doesn't appear in the recipient's **File
   Manager**. Instead, clicking the link _opens_ the file in Datagrok. From there, users can either download the file
   and then upload it to Datagrok or save the file as a [project](../datagrok/create-project.md).

   To access the file or folder from the link provided, users must have access privileges for this file or folder.

### Share files

When your access privileges allow it, you can share folders and files available to you. You have several options:

* Share a _connection_ (root folder) to give access to the entire directory.
* Share a _folder_ to give access to the content of individual folders in your directory.
* Share a _file_ to limit access to individual files.<!--future-looking-->

To share an item, do the following:

1. Right-click the item and select _share_ from the context menu. This action opens the **Share...** dialog.
2. Enter the user or user group you want to share it with.

   In the identity/email field, start typing a person’s name, username, email, or group name. Pick from the list of the
   matching identities.

3. From the respective dropdowns, select access privileges for either or both: (1) the connection and (2) individual
   files/folders. You can select any or all of the following options<!--TBU-->:

   * _Can view_: Users can view, open, and download
   * _Can edit_: Users can rename and edit
   * _Can delete_: Users can delete
   * _Can share_: Users can reshare with any other user or group.<!--how does it work with URL links? -->

   > Note: For each _file share_, Datagrok automatically assigns a _friendly name_ (a name displayed in the UI) that
   > specifies all the directory names starting from the root folder. While you _can_ change the _friendly name_, we
   > encourage you to exercise caution. Changing the item's _friendly name_ also changes its unique
   > _namespace-qualified name_ used in URL links, automation workflows, and more.

4. Optionally, enter a description in the text field provided. You may also notify the users you share with. If you
   don’t want to send a notification, uncheck the **Send notification** checkbox.

   > Note: To send an email notification, enter the user's email in the identity/email field. The email notification
   > contains a link to the shared item and entered description. If you enter a user or group name, they will be
   > notified via the Datagrok interface.

5. Click **OK** to share. Once shared, the shared item appears in the recipient's **Folder Tree**.

   ![Share a folder](share-the-folder.gif)

   > Note: When you share a file or a folder, Datagrok automatically indexes folders <!--link out to indexing when ready-->
   > and working--> and extracts basic metadata (like the date created or its size). File indexing is optional (to index
   > files, toggle **Index Files** in the **Share...** dialog).

## Manage file shares

Subject to your privileges, you can use the **Data Explorer** panel to inspect and quickly adjust access permissions to
file shares, send comments to those you're sharing with, and more.

1. First, select the connection, file, or folder.
2. Then navigate to **Data Explorer** on the left and expand the **Sharing** info panel to see the complete list of
   users with access and their privileges. From here, you can click any user or group to see their profile, your
   conversation history with them, send them a note, and more.
3. Use the buttons provided to share access with more users, revoke access, or edit permissions to the shared item.

   > Note: The same actions are available from the context menu.

<!--TBD: GIF pending changes in the UI-->

## Managing folders and files

The **File Manager** lets you perform standard file and folder actions such as open, download, delete, or rename. To see
a complete list of available actions, right-click the file or folder.

> Tip: The same list of actions is available from the **Data Explorer** on the right under the **Actions** info panel.

Depending on your privileges, certain actions may not be available to you. For files and folders shared with you,
contact the _credentials owner_. If you are a _credentials owner_, contact the _data source_ owner.

## Browse and preview

**File Manager** has multiple views to help you navigate and browse the content of your directory:

1. A **Folder Tree** for hierarchical browsing
2. A **Directory View** with three viewing modes (icons, cards, and grid) and built-in search
3. A dynamic file/folder **Preview** (here, you can also do lightweight editing for CSV and TXT file formats).

   You can reposition or close any pane<!--link to article when ready-->. You can also hide the **Preview** pane by
   clicking the **Toggle file preview** icon above the **Folder Tree**.

   <!--gif-->

   Clicking a file or a folder opens its preview, double-clicking a file opens it in the platform, and double-clicking a
   folder expands its content.

   > Tip: Use built-in [smart search](../datagrok/smart-search.md) to find files and folders of interest quickly. You can
   > type a word or part of a word to search for the following:
   >
   >* File and folder names
   >* File extension
   >* File and folder content metadata.<!--is this a complete list?-->

Datagrok has multiple tools to help you find information you need precisely when you need it.

* A [Treemap](../visualize/viewers/tree-map.md) helps you visually profile the content of your current folder and identify
  the files that take up
  the most space.
  > Developers: You can [create custom folder viewers](../develop/how-to/folder-content-preview.md).
  >
  ><!--gif clinical data-->
  >
  >In the example above, a script is invoked against the folder content. If the folder contains files matching the file
  > extension parameter, the **Preview** shows a custom widget (in this case - the application launch button) every time
  > the folder is opened.

* You can browse and preview the content of individual files in [30+ formats](supported-formats.md). What you see is
  dynamically adjusted based on the file type. For example:

  * You can view the content of ZIP files just like regular folders.
  * View and edit Markdown for HTML files.
  * View image data in cell output and as files preview.
  * Scroll through an interactive spreadsheet for CSV files.
  * Look at protein structure from multiple angles using a 3D protein browser.

   <!--gif-->

  > Developers: You can add custom formats via [package extensions](../develop/how-to/create-package.md).
  > Note: File preview is limited to files less than 10MB in size. The platform won't display bigger files.
  >
  >You can't preview unsupported file formats, but you can download these files.

* You can get on-demand information about folders and files, and data within these files:
  * Hover over objects to see context-driven tooltips.

  * Select an object and navigate to **Data Explorer** on the left to see all applicable [info
    panels](../discover/info-panels.md) (for example, **Details** shows basic metadata about files and folders, **Stats**
    shows summary column statistics for tabular datasets, and **Chem** lets you quickly find molecules of interest using
    sketcher).
  > Developers: You can augment file preview with custom metadata and [create organization-specific file
  > viewers](../develop/how-to/create-custom-file-viewers.md).
  >
  > ![Cell image segmentation](Cell-image-segmentation.gif)
  >
  >In the example above, a Python script is invoked against the JPEG and JPG files, cells are automatically segmented in
  > file preview, and the number of detected cell segments is displayed on the info panel.

## Customizations

Fundamentally, Datagrok is designed as an extensible environment. Datagrok extensions can customize or enhance any part
of Datagrok. You can provide custom UI, file and folder viewers and editors, or [cell
renderers](../develop/how-to/custom-cell-renderers.md) for rich preview. You can add new data types and extract
organization-specific metadata. Extensions can add items to the menus and [context
actions](../develop/how-to/context-actions.md), and so much more.

To learn more about extensions, see [Extending and customizing Datagrok](../develop/packages/extensions.md).

## Automation

Any action performed in **File Manager** (other that uploading local files) is reproducable and can be used in
automation workflows. For example, you can use data preparation pipeline to define jobs for data ingestion,
postprocessing, and transformations.

To learn more about automating workflows usung data preparation pipelines, see [Data preparation
pipeline](data-pipeline.md)

## Resources

<!--insert link to YouTube video-->

See also:

<!--* [Databases](link)-->

<!--* [Web services](link)-->

<!--*[Data explorer](link)-->

<!--*[Indexing](link)-->
