# Package Guidelines

## Naming convention

Use kebab-case for package name (`power-pack`) and pascal-case for folder name (`PowerPack`).
Use kebab-case for file names.

## Metadata

`package.json` file must contain:
1. Proper package name with scope:

`@datagrok/power-pack` or `@datagrok-libraries/utils`
2. Description field
3. Repository section:
```
  "repository": {
    "type": "git",
    "url": "https://github.com/datagrok-ai/public.git",
    "directory": "packages/CustomML"
  }
  ```
4. Dependencies section with `datagrok-api` dependency. Refer to Semantic Versioning concept
5. Proper package version, according to Semantic Versioning concept

## Tests

## Code

### Linter