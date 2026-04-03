---
paths:
  - "**/src/tests/**"
  - "**/package-test.*"
---

## Datagrok Test Conventions

Tests run in a Puppeteer-controlled browser against a running Datagrok instance.

```bash
grok test                          # Run tests against default server
grok test --host localhost         # Test against local instance
grok test --gui                    # Visual browser (not headless)
grok test --gui --debug            # Debug breakpoints
grok test --test "TestName"        # Run specific test by prefix
grok test --category "Category"    # Run tests in a category
grok test --skip-build             # Skip building before test
grok test --verbose                # Detailed output
```

Test files use `@datagrok-libraries/test` utilities. The test entry point is `src/package-test.ts` which exports a tests array.
