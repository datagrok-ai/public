# Open files

## Open a file from a local file system

One of the easiest ways to open a file is to simply drag-and-drop it from a file system. After you
start dragging, there is a message that you can drop it here, and just drop it. In a second the file
is opened and available in the platform.

Among many data ingestion steps available within the platform, opening a file is unique in a sense
that it is the only one which cannot be automated. This is because the browser does not have access
to a local computer's file system. This cannot be easily recorded and reproduced. All the other
steps are fully reproducible.

## Open a file from a file share

Additionally files could be opened in a different manner. There is a whole section for data access
on top of the left bar with a "Folder" 📁 icon. After clicking on it, different options may be
explored. It starts with files living in a file system that is accessible from a server. When one of
these folders is open, what is happening behind the scenes is that Datagrok which runs in your
browser is talking to the server which is installed on Amazon (or in your private cloud, if you are
working with your company's instance of Datagrok). This server has access to a number of folders
that are shared.

In order to add a share, click the "New File Share", then we can connect to multiple types of file
storages, such as Dropbox, Google Cloud, S3 and others. Then, depending on the storage type you are
asked to enter login credentials. Once everything is entered, the connection will be created.

The platform allows us to explore the file system just like you would do regularly with your file
system browser. There is a hierarchical tree of directories. The files could be dragged and dropped
from your desktop to the corresponding server locations. You can right-click and see available
actions. It's a file management platform as well which you can control your file system with.

It is all covered with a system of priviledges. A connection could be shared to multiple individuals
or multiple groups with appropriate set of priveledges.

The way the browser is structured is that on top it has a list of files and folders in it, and then
there is also a preview of the current folder. By default it shows us what takes the most space,
using tree map which is area-coded by the file size. As we click on the file, the content gets shown
interactively.

To bring the file to the platform, you would either double-click the file, or right-click and
select `Open`. This particular action of opening a file, as opposed to the action of opening a file
from a local file system with a drag-and-drop, is replayable. Behind the scenes, in [Console]() we
can find that the action was recorded with a name `OpenServerFile` and a path to it. A path is
qualified with a namespace. What is great about it is that it could be replayed at any time as part
of some other script.

## Import file from text

...