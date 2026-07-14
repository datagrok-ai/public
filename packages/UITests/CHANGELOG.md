# UI Tests changelog

## v.next

* Tests: name the projects created in `tag.add` and `projects.api` — projects now require a non-empty name on save (datlas projects_service guard), which was throwing `Project name cannot be empty`.
* Tests: `pca` now skips gracefully when the ML menu (EDA package) isn't deployed, instead of crashing with `Cannot read properties of undefined (reading 'find')`.
* Tests: `UI: Sharing` uses a unique per-run entity name so a leftover from a prior aborted/failed run can't make the gallery search return >1 result ("more than one testing entity present").
* Tests: `View: Events` onViewLayoutApplying/onViewLayoutApplied now search the whole dock tree for the viewer — they only checked the top-level `children` (splitters), never the nested viewer leaves, so `viewerElement` was always null (chronic failure on dev too).

## 1.0.12 (2024-09-02)

* Test fixes update 

## 1.0.4 (WIP)
