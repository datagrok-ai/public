# WebComponents changelog

## v.next

- DockManager: Esc no longer closes the active panel
- Viewer: An options change no longer re-sets the frame and fires `viewer-data-frame-changed` for frames created in JS
- Viewer: A frame set on the viewer itself is kept on later options changes, and its `viewer-data-frame-changed` fires outside the Dart event
- Viewer: Removed the unused `viewer-type-changed` event and `getViewerEventObservable`
- Viewer: A failed viewer creation, or a throwing frame or options apply, no longer stops the element from applying later changes
- InputForm: `destroy()` cancels a form still being built
- ToggleInput: A value set by the host is no longer sent back as `value-changed`
- ValidationIcon: A `disabled` or `restricted` input whose value is consistent shows a link icon
- InputForm: Added the `skipLogic` property
- DockManager: Fixed views opened, or re-rendered, in the background splitting their panels evenly instead of by the declared ratios
- ValidationIcon: Advice actions send their `additionalParams` with the action request
- InputForm: Fixed a form opened with invalid or missing required values reporting itself valid

## 1.4.6 (2026-09-23)

- Bundles webcomponents 0.3.7: InputForm funcCall swap and cleared-funcCall fixes, ValidationIcon retention, clearing and severity colors

## 1.4.5 (2026-08-14)

- Viewer: Fixed `viewer-data-frame-changed` firing inside the Dart event (grid mutations in handlers crashed)

## 1.4.4 (2026-08-07)

- Viewer: apply dataframe and options set in one tick as a single transition, fixing look/frame crashes (e.g. split column deselect) and removing the property-order dependency

## 1.4.3 (2026-07-02)

- Fixed validation popover scrolling, scrollbar visibility, and above/below positioning

## 1.3.0 (2025-03-19)

- 1.27 platform support

## 1.2.10 (2025-11-20)

- Tweak consistency icons in forms
- Tweak dock spawn webcomponent api

## 1.2.9 (2025-10-30)

- Update deps

## 1.2.8 (2025-09-09)

- Support restoring focused tab in dock

## 1.2.7 (2025-09-04)

- Small bug fixes

## 1.2.6 (2025-08-28)

- Compute2 related updates

## 1.2.5 (2025-07-28)

- API updated for 1.26.0

## 1.2.4 (2025-06-03)

- Allow all users to use

## 1.2.0 (2025-03-30)

- API updated for 1.25.0

## 1.1.0 (2025-02-19)

- Release for 1.24.0

## 1.0.0 (2024-12-25)

- Initial release
