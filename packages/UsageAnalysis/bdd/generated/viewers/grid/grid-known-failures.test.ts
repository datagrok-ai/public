/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1198, core#1195, core#1193, core#1194]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, pressKey, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, filterBetween, noneSelected, rowCount, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, pickFromAreaContextMenu, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid known failures", () => {
  const session = feature(test, "features/viewers/grid/grid-known-failures.feature", import.meta.url);
  test("A malformed number format marks the format field invalid", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:22", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(19, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(20, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(21, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(22, "Then Custom input in Format AGE dialog should be invalid", () => shouldBe(page, el("Custom input in Format AGE dialog"), "invalid"));
    }, 22);
  });
  test("A malformed number format is not taken and the old format stays", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:30", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(26, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(27, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(28, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(29, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(30, "Then \"AGE\" column should have tag \"format\" equal to \"\"", () => columnTag(page, "AGE", "format", ""));
    }, 30);
  });
  test("A column group named with spaces only is not made", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:42", "@realizes:core#1195"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(34, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(35, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(36, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(37, "And Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(38, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(39, "Then Group input in dialog should be visible", () => shouldBe(page, el("Group input in dialog"), "visible"));
      await session.step(40, "When user types \"   \" into Group input in dialog", () => typeInto(page, "   ", el("Group input in dialog")));
      await session.step(41, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(42, "Then \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
    }, 42);
  });
  test("A read-only grid keeps its rows on Shift+Delete", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:50", "@realizes:core#1193"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(46, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
      await session.step(47, "And user clicks on the \"row header 5\" area of grid holding Control", () => clickAreaHolding(page, "row header 5", el("grid"), "Control"));
      await session.step(48, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(49, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(50, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    }, 50);
  });
  test("Shift+Enter selects nothing while Allow Row Selection is off", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:61", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(54, "When user clicks on the \"cell 4 of SEX\" area of grid", () => clickArea(page, "cell 4 of SEX", el("grid")));
      await session.step(55, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(56, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(57, "When user clears the row selection", () => clearSelection(page));
      await session.step(58, "And user sets \"Allow Row Selection\" property of grid to \"false\"", () => setProperty(page, "Allow Row Selection", el("grid"), "false"));
      await session.step(59, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(60, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(61, "Then no rows should be selected", () => noneSelected(page));
    }, 61);
  });
  test("Shift+Enter selects only the same rows that pass the filter", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:70", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(65, "When user filters rows where \"AGE\" is between 20 and 30", () => filterBetween(page, "AGE", 20, 30));
      await session.step(66, "And user clicks on the \"cell 2 of SEX\" area of grid", () => clickArea(page, "cell 2 of SEX", el("grid")));
      await session.step(67, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(68, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(69, "Then some rows should be selected", () => someSelected(page));
      await session.step(70, "And 82 rows should be selected", () => selectedRowCount(page, 82));
    }, 70);
  });
});
