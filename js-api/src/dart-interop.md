# Dart Interop Guide

This document explains how the Datagrok JavaScript API communicates with the Dart backend.

## Overview

The Datagrok platform is built with Dart on the backend. The JavaScript API provides TypeScript/JavaScript bindings that wrap Dart objects and communicate through an interop layer.

## The `dart` Property Convention

Every JavaScript wrapper class holds a reference to its underlying Dart object in a property called `dart`:

```typescript
class DataFrame {
  public readonly dart: any;  // Dart handle

  constructor(dart: any) {
    this.dart = dart;
  }
}
```

When working with the API:
- **Do not** modify the `dart` property directly
- **Do** pass wrapper objects to API methods (they will extract the `dart` handle automatically)
- The `dart` property is primarily used internally for calling Dart methods

## Type Conversion: `toJs()` and `toDart()`

Two functions handle conversion between Dart and JavaScript objects:

### `toJs(dart, check?)` - Dart to JavaScript

Converts a Dart object to its JavaScript wrapper. The function handles these types:

| Dart Type | JavaScript Result |
|-----------|-------------------|
| `null` | `null` |
| `undefined` | `undefined` |
| `FLOAT_NULL` | `null` (special float null sentinel) |
| `Map` | Plain object with recursively converted values |
| `List` | Array with recursively converted values |
| `DateTime` | `dayjs` object |
| `ByteArray` | Passed through unchanged |
| `BigInt` | JavaScript `BigInt` |
| `Property` | `Property` wrapper |
| Objects with wrappers | Returns the existing JavaScript wrapper |

If `check` is `true` and the type cannot be converted, an exception is thrown.

### `toDart(x)` - JavaScript to Dart

Converts a JavaScript object to its Dart representation:

| JavaScript Type | Dart Result |
|-----------------|-------------|
| `null`/`undefined` | Passed through |
| `dayjs` object | Dart `DateTime` |
| Object with `toDart()` method | Result of calling `toDart()` |
| Object with `dart` property | The `dart` property value |
| Plain object `{}` | Dart `Map` |
| `BigInt` | Dart `BigInt` |
| Primitives | Passed through unchanged |

## The `wrappers.ts` vs `wrappers_impl.ts` Split

There are two wrapper files due to bootstrapping requirements:

- **`wrappers.ts`**: Delegates to `DG.toJs()` and `DG.toDart()`. Used during initial load when the `DG` global is being set up.
- **`wrappers_impl.ts`**: Contains the actual implementation. The `DG` global eventually calls these functions.

This split allows the API to work during the bootstrapping phase before all globals are initialized.

## API Function Naming Convention

Dart functions exposed to JavaScript follow a naming convention:

```
grok_<ClassName>_<MethodName>
grok_<ClassName>_Get_<PropertyName>
grok_<ClassName>_Set_<PropertyName>
```

Examples:
```typescript
api.grok_DataFrame_Columns(dart)        // Get columns from DataFrame
api.grok_Entity_Get_Id(dart)            // Get id property of Entity
api.grok_Entity_Set_Id(dart, value)     // Set id property of Entity
api.grok_BigInt_To_BigIntJs(dart)       // Convert BigInt to JS BigInt
```

These functions are defined in `src/api/grok_api.g.ts` (auto-generated from Dart code).

## The `api` Object

Each source file accesses Dart functions through an `api` object:

```typescript
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;
```

This is because Dart interop functions are attached to the global `window` object. The `IDartApi` interface (in `grok_api.g.ts`) provides type hints for these functions.

## Event System

Events are bridged using RxJS observables. The `__obs()` function creates an observable from a Dart event stream:

```typescript
function __obs<T>(eventId: string, object?: any): Observable<T>
```

Event IDs follow patterns like:
- `d4-*` - Platform events (e.g., `d4-context-menu`, `d4-table-added`)
- `grok-*` - Grok-specific events (e.g., `grok-view-added`, `grok-project-saved`)
- `ddt-*` - DataFrame events (e.g., `ddt-values-changed`, `ddt-selection-changed`)

## Common Patterns

### Getting a property
```typescript
get name(): string {
  return api.grok_Entity_Get_Name(this.dart);
}
```

### Setting a property
```typescript
set name(value: string) {
  api.grok_Entity_Set_Name(this.dart, value);
}
```

### Calling a method
```typescript
close(): void {
  api.grok_View_Close(this.dart);
}
```

### Converting results
```typescript
get columns(): ColumnList {
  return toJs(api.grok_DataFrame_Columns(this.dart));
}
```

### Passing arguments
```typescript
addColumn(column: Column): Column {
  return toJs(api.grok_DataFrame_Add_Column(this.dart, column.dart));
}
```
