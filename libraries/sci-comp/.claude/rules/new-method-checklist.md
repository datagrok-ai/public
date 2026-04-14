# Adding New Methods Checklist

When adding any new numerical method to this library, ALL of the following are required. Do NOT report the task as complete until every item is done.

1. **Tests** — full test file in `__tests__/` following the reference test structure for that domain (sync + async variants where applicable)
2. **TSDoc** — all public types, interfaces, functions, and classes must have TSDoc comments
3. **Section README** — update the README.md in the same directory (or nearest parent directory) as the new method's source files
4. **Project README** — update the top-level `README.md`
5. **CLAUDE.md** — update the architecture tree with new files
6. **Verification** — `npm run lint-fix && npm run build && npm test` must all pass
