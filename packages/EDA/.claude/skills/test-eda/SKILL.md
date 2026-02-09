---
name: test-eda
description: Run grok test for the EDA package against a server host and report results
argument-hint: <EDA> [--host host] [--category cat] [--test name] [--verbose] [--skip-build] [--no-retry]
---

# Test EDA

Run `grok test` for the EDA package and provide a summary of results.

## Instructions

1. Parse the user's input to extract:
   - `--host` value (default: `local`)
   - Any extra `grok test` flags (`--category`, `--test`, `--verbose`, `--skip-build`, `--gui`, `--no-retry`, `--benchmark`, `--tag`)

2. Run tests from the package EDA directory using the Bash tool:
   ```
   grok test --host <host> --no-retry [extra flags]
   ```
   Use a timeout of 600000ms (10 minutes). Print output to console.

3. After all packages finish, provide a **summary table** with:
   - Total passed / failed / skipped counts (parse csv output file test-results.csv)
   - List of any failed categories and test names
   - Overall verdict: PASSED or FAILED

4. Cleanup csv output file.