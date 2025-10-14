---
title: "Package API"
---

# Strongly Typed Function Calls in Grok

There is a way to call package functions using the `grok.functions.call` method.  
Unfortunately, this approach is not reliable because it is not strongly typed, making it easy to introduce type errors in the input.

To enforce strict typing, packages can share their API.

## Generating a Strongly Typed API

The easiest way to generate an API is by using the following command:

```bash
    grok api
```

This will create a `package-api.ts` file with all strongly typed functions.  
From this point forward, you can use the API within the current package.

## Making the API Available to Other Packages

To make the package API accessible from other packages, follow these steps:

- Update `.npmignore`

Ignore all `.d.ts` files except the generated `package-api.d.ts`:


```text
    **/*.d.ts
    !src/package-api.d.ts
```

- Update `.gitignore`

Ignore all `.d.ts` files to avoid committing unnecessary declaration files:

```text
    **/*.d.ts
```

- Update `tsconfig.json`

Enable TypeScript to generate declaration files by modifying `tsconfig.json`:

```json
    {
      "compilerOptions": {
        "declaration": true
      }
    }
```

## Using the Shared API

Once all the above steps are complete, you can:

- Add the package as a dependency in other projects  
- Use its **strongly typed API** safely and reliably

This setup ensures better type safety and developer experience across packages.