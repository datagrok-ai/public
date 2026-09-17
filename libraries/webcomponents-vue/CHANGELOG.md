# webcomponents-vue changelog

## v.next

* Introduced an injectable per-view service (`provideDgViewService`/`useViewService`/`useDgView`) that owns ribbon rendering and diffing
* RibbonPanel/RibbonMenu: rewritten as thin descriptor-based wrappers over the view service; `view` prop removed, added `priority` ordering and built-in disabled state (debounced style, reason tooltip/popup)
* wheelGuard: hint overlay colors use design tokens
* ifOverlapping: a show scheduled within the debounce window no longer lands after unmount; the loader is removed on unmount

## 0.3.4 (2026-08-07)

* `Viewer`: patch `dataFrame` before `options` so options referencing new columns don't hit the old frame
* Added `ResizeHandle` component (deferred mode, non-live by default)
* Added `wheelGuard` directive

## 0.3.1 (2026-06-03)

* Default float format for funccall `InputForm` doubles; exposed default float mask
* Extracted `useUnwrappedCallMeta` composable
* `DockManager`: keep preferred tab instead of switching on tab add; tab persistence for hidden outputs
