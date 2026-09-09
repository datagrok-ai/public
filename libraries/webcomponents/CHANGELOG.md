# webcomponents changelog

## v.next

* InputForm: Fixed rapid funcCall swaps racing (a stale form could land in the DOM and win over the newer one)
* InputForm: Setting funcCall to undefined no longer permanently kills the input/validation event streams

## 0.3.6 (2026-08-11)

* Viewer: Fixed `viewer-data-frame-changed` firing inside the Dart event (grid mutations in handlers crashed)

## 0.3.5 (2026-08-07)

* Viewer: apply dataframe and options set in one tick as a single transition (new `ViewerHost` core), fixing look/frame crashes and removing the property-order dependency

## 0.3.4 (2026-07-02)

* Fixed validation popover scrolling, scrollbar visibility, and above/below positioning

## 0.3.2 (2026-06-03)

* Fixed `viewer-data-frame-changed` event detail
* Annotated `ValidationIcon` for click tracking
