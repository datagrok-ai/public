# webcomponents-vue changelog

## 0.3.4 (2026-08-07)

* `Viewer`: patch `dataFrame` before `options` so options referencing new columns don't hit the old frame
* Added `ResizeHandle` component (deferred mode, non-live by default)
* Added `wheelGuard` directive

## 0.3.1 (2026-06-03)

* Default float format for funccall `InputForm` doubles; exposed default float mask
* Extracted `useUnwrappedCallMeta` composable
* `DockManager`: keep preferred tab instead of switching on tab add; tab persistence for hidden outputs
