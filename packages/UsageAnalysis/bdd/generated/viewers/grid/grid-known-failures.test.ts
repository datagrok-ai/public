/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1200, core#1198, core#1195, core#1193, core#1194]
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
import {clickArea, clickAreaHolding, hoverArea, pickFromAreaContextMenu, propertyShouldBe, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid known failures", () => {
  const session = feature(test, "features/viewers/grid/grid-known-failures.feature", import.meta.url);
  test("Row Height typed in the context panel as 0 is not taken [typed=0]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:25", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(19, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(20, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(21, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(22, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(23, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(24, "When user enters \"0\" into \"Row Height\" property in context panel", () => enterInto(page, "0", el("\"Row Height\" property in context panel")));
      await session.step(25, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, 25);
  });
  test("Row Height typed in the context panel as 501 is not taken [typed=501]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:25", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(19, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(20, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(21, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(22, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(23, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(24, "When user enters \"501\" into \"Row Height\" property in context panel", () => enterInto(page, "501", el("\"Row Height\" property in context panel")));
      await session.step(25, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, 25);
  });
  test("A malformed number format marks the format field invalid", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:37", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(34, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(35, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(36, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(37, "Then Custom input in Format AGE dialog should be invalid", () => shouldBe(page, el("Custom input in Format AGE dialog"), "invalid"));
    }, 37);
  });
  test("A malformed number format is not taken and the old format stays", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:45", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(41, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(42, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(43, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(44, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(45, "Then \"AGE\" column should have tag \"format\" equal to \"\"", () => columnTag(page, "AGE", "format", ""));
    }, 45);
  });
  test("A column group named with spaces only is not made", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:57", "@realizes:core#1195"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(49, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(50, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(51, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(52, "And Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(53, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(54, "Then Group input in dialog should be visible", () => shouldBe(page, el("Group input in dialog"), "visible"));
      await session.step(55, "When user types \"   \" into Group input in dialog", () => typeInto(page, "   ", el("Group input in dialog")));
      await session.step(56, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(57, "Then \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
    }, 57);
  });
  test("A read-only grid keeps its rows on Shift+Delete", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:65", "@realizes:core#1193"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(61, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
      await session.step(62, "And user clicks on the \"row header 5\" area of grid holding Control", () => clickAreaHolding(page, "row header 5", el("grid"), "Control"));
      await session.step(63, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(64, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(65, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    }, 65);
  });
  test("Shift+Enter selects nothing while Allow Row Selection is off", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:76", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(69, "When user clicks on the \"cell 4 of SEX\" area of grid", () => clickArea(page, "cell 4 of SEX", el("grid")));
      await session.step(70, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(71, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(72, "When user clears the row selection", () => clearSelection(page));
      await session.step(73, "And user sets \"Allow Row Selection\" property of grid to \"false\"", () => setProperty(page, "Allow Row Selection", el("grid"), "false"));
      await session.step(74, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(75, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(76, "Then no rows should be selected", () => noneSelected(page));
    }, 76);
  });
  test("Shift+Enter selects only the same rows that pass the filter", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:85", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(80, "When user filters rows where \"AGE\" is between 20 and 30", () => filterBetween(page, "AGE", 20, 30));
      await session.step(81, "And user clicks on the \"cell 2 of SEX\" area of grid", () => clickArea(page, "cell 2 of SEX", el("grid")));
      await session.step(82, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(83, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(84, "Then some rows should be selected", () => someSelected(page));
      await session.step(85, "And 82 rows should be selected", () => selectedRowCount(page, 82));
    }, 85);
  });
});
