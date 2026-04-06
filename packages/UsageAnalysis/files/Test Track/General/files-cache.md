1. Go to Files > Demo. Check if there is cache enabled for Demo connection. (Open context menu on Demo connection, press 'Edit...', check 'Cache' checkbox)
2. Create new folder in Files > Demo root directory with name "Folder cache test"
3. Create cache mapping for newly created directory. Open context menu on Demo connection, press 'Cache...', choose new foolder, add cron invalidtion: `*/2 * * * *` (every 2 minutes)
4. Create file test.txt in "Folder cache test" folder
5. Write "Hello world!" in the file and save it.
6. Rename file to 'test1.txt'
7. Rename folder to "Folder cache test1"
8. Delete "Folder cache test1" folder. Make sure the folder is deleted.
9. Check cache mappings, the mapping for created folder should be deleted also (Context menu, Cache...)

---
{
  "order": 4
}
