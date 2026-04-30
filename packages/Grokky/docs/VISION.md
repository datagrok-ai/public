# Agentic force for scientists

## Plugins

Plugins ship special instructions (collection of .md files?) that describe the plugin's
functionality and what could be done with it for the user.

Example 1: Chem plugin might ship a document with this

```
for chemical space, call the ChemicalSpace function
```

Example 2: Company-specific: A plugin (or file) with company-specific procedures or
instructions that get included to the session

```
When registering a new chemical compound, do the following
- make sure it's a valid structure
- ask user to associate it with a project
- register in Global Compound Registration System
- kick off ADME prediction
- email John Smith the link to ADME results
```

At the Grokky Claude Code session startup, these folders get loaded in the Claude Code
context. Perhaps as one CLAUDE.md with the index of available plugins and short
descriptions.

So this becomes both a mechanism that defines what is available to the user (via plugin
privileges), and a mechanism for efficient working.

---

For the implementation of the skills & knowledge sync system, see
[IMPLEMENTATION.md](IMPLEMENTATION.md).
