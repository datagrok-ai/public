#### Databases node visibility & Test packages check

- Verify that only the required demo databases are visible in the public environment, and that no test packages or test connections affect the user experience.

- **Preconditions:**
- Access to the public environment.
- Databases is visible in the Browse.
- Demo databases should be: **Chembl, Northwind, Starbucks, World,** and optionally **SureCHEMBL**.

#### Steps:

1. Public: Check installed packages: ensure that DbTests and any other test packages are not installed.
2. Navigate to the Databases in Browse. Perform the initial expand of the Databases node:
- Only providers with at least one active connection should be displayed.
- Providers without any connections must not be shown in the first expand.
- Following databases are visible: Chembl, Northwind, Starbucks, World, SureCHEMBL (if present, should be properly named).
- No test databases (e.g., Athena, other disconnected or experimental databases) are shown.

