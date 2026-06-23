File-share folder/file CRUD lifecycle on a connection:

1. Create a new folder "Folder cache test" in a connection's directory.
2. Create file test.txt in the folder and write "Hello world!" into it.
3. Rename the file to 'test1.txt'.
4. Rename the folder to "Folder cache test1".
5. Delete the "Folder cache test1" folder. Make sure the folder is deleted.

> Automated by `files-cache-spec.ts` (folder/file create-write-read-rename-delete
> lifecycle via `grok.dapi.files`).
>
> The connection-level **cache** configuration — enabling "Cache" on the
> connection, creating a "Cache..." mapping with cron invalidation
> `*/2 * * * *`, and verifying the mapping is removed together with its folder —
> is manual (no JS API): see `files-cache-ui.md`.

---
{
  "order": 4
}
