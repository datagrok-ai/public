## Plugin changelog rule

When committing a sizable feature or bug fix to a plugin under `packages/`,
add a one-line entry to the `## v.next` section at the top of that
plugin's `CHANGELOG.md`. Create the section if it doesn't exist yet.

Format:

```markdown
# PackageName changelog

## v.next

* Component: Added new feature description
* GROK-12345: Fixed something

## 1.17.2 (2026-03-23)
...
```

- One bullet per logical change, flat list (no subsections)
- Start with ticket ID if available (`GROK-NNNNN:` or `[#NNN](url):`)
- Use past-tense verbs: Added, Fixed, Improved, Introduced, Implemented
- Skip trivial changes (typos, formatting, dependency bumps)
- When the plugin is published, `/plugin-changelog` converts `v.next` to a versioned entry
