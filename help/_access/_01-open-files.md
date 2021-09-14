# Open files

One of the easiest ways to open a file is to simply drag-and-drop it from a file system. After
you start dragging, there is a message that you can drop it here, and just drop it. In a second
the file is opened and available in the platform.

Among many data ingestion steps available within the platform, opening a file is unique in a
sense that it is the only one which cannot be automated. This is because the browser does not
have access to a local computer's file system. This cannot be easily recorded and reproduced.
All the other steps are fully reproducible.

Additionally files could be opened in a different manner. There is a whole section for data
access on top of the left bar with a "Folder" 📁 icon. After clicking on it, different options may
be explored. It starts with files living in a file system that is accessible from a server. When
one of these folders is open, what is happening behind the scenes is that Datagrok which runs in
your browser is talking to the server which is installed on Amazon (or in your private cloud, if
you are working with your company's instance of Datagrok). This server has access to a number
of folders that are shared.

In order to add a share, click the "New File Share", then we can connect to multiple types of
file storages, such as Dropbox, Google Cloud, S3 and others. Then, depending on the storage type
you are asked to enter login credentials. Once everything is entered, the connection will be
created.