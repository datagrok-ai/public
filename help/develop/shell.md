<!-- TITLE: Datagrok Shell -->

# Datagrok Shell

Datarok visual shell is used to get access to top-level views, tables, methods,
and platform states.

## Error handling

Use `grok.shell.lastError` to get the last unhandled error state which the platform is in.

In the Datagrok UI this information is found in the left-bottom corner, where access
to help is usually provided via a question mark icon. If the platform intercepts an exception
not handled elsewhere (by the platform itself or by the applied code), or gets into some
other erroneous states, this question sign turns into a red exclamation sign. Hovering
over this sign an error message may be found.

To know that the platform is in such erroneous state, use `grok.shell.lastError`
read-only property. If there is no error, the property returns an empty string.
If there is an error, the property returns a non-empty error description.
