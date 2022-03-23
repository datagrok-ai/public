# Frequently asked questions

*Question:*
How do I set up VS Code or WebStorm for debugging on Datagrok? Debugging doesn't work for me now.

*Answer:*
Check out our guides on [VS Code] and [WebStorm]. If you configured VS Code or WebStorm following the guides, but
debugging still doesn't work, ask our technical gurus on the [Community Forum].

*Question:*
I've installed `npm` and `datagrok-tools` packages, but `grok` tools aren't available on my paths.

*Answer:*
It's likely that the `datagrok-tools` library wasn't installed globally. Install `datagrok-tools`
globally:

```sh
# First, remove the locally installed package...
npm uninstall datagrok-tools

# And then install datagrok-tools globally using the -g flag
npm install -g datagrok-tools
```

*Question:*
Others don't see my published package on the selected Datagrok instance. How can I fix this?

*Answer:*
If you publish a package in debug mode using the standard `grok publish` command, only you will see the published
package. Run `grok publish --release` so others could see your package, too.
