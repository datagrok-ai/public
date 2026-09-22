/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, noneSelected, selectAllRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, loadLayout, noBalloons, noErrors, propertiesShouldBe, propertyShouldBe, readingIs, readingReads, saveLayoutToServer, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot through a layout and a project", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-persistence.feature", import.meta.url);
  test("Two trellises keep their own settings through one layout", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(24, "Then the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(27, "Given user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
    await session.step(28, "When user sets properties of second trellis plot viewer:", () => setProperties(page, el("second trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]);
    await session.step(32, "Then the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
    await session.step(33, "And properties of first trellis plot viewer should be:", () => propertiesShouldBe(page, el("first trellis plot viewer"), [["X Column Names","SEX, CONTROL"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Viewer Type","Pie chart"]]);
    await session.step(36, "And properties of second trellis plot viewer should be:", () => propertiesShouldBe(page, el("second trellis plot viewer"), [["X Column Names","RACE"],["Viewer Type","Bar chart"]]), [["X Column Names","RACE"],["Viewer Type","Bar chart"]]);
    await session.step(39, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(40, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(41, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(42, "When user loads the saved layout", () => loadLayout(page));
    await session.step(43, "Then scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
    await session.step(44, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
    await session.step(45, "And properties of first trellis plot viewer should be:", () => propertiesShouldBe(page, el("first trellis plot viewer"), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(49, "And properties of second trellis plot viewer should be:", () => propertiesShouldBe(page, el("second trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]);
    await session.step(53, "When user sets \"Viewer Type\" property of first trellis plot viewer to \"Histogram\"", () => setProperty(page, "Viewer Type", el("first trellis plot viewer"), "Histogram"));
    await session.step(54, "Then the \"inner viewer type\" reading of first trellis plot viewer should be \"Histogram\"", () => readingReads(page, "inner viewer type", el("first trellis plot viewer"), "Histogram"));
    await session.step(55, "And the \"inner viewer type\" reading of second trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("second trellis plot viewer"), "Bar chart"));
    await session.step(56, "When user sets \"Viewer Type\" property of first trellis plot viewer to \"Pie chart\"", () => setProperty(page, "Viewer Type", el("first trellis plot viewer"), "Pie chart"));
    await session.step(57, "And user clicks on close icon of second trellis plot viewer", () => clickOn(page, el("close icon of second trellis plot viewer")));
    await session.step(58, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(59, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(63, "And the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("A selection behind Row Source Selected does not come back with the project, and the trellis still draws", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(24, "Then the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(67, "When user sets \"On Click\" property of trellis plot viewer to \"Select\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Select"));
    await session.step(68, "And user clicks on the \"cell F, false | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F, false | Caucasian", el("trellis plot viewer")));
    await session.step(69, "Then 475 rows should be selected", () => selectedRowCount(page, 475));
    await session.step(70, "When user sets \"Row Source\" property of trellis plot viewer to \"Selected\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Selected"));
    await session.step(71, "Then trellis plot viewer should show 475 rows", () => showsRows(page, el("trellis plot viewer"), 475));
    await session.step(72, "When user saves the current view as project \"bdd-trellis-selected\"", () => saveAsProject(page, "bdd-trellis-selected"));
    await session.step(73, "And user closes all views", () => closeAllViews(page));
    await session.step(74, "And user opens the \"bdd-trellis-selected\" project", () => openProject(page, "bdd-trellis-selected"));
    await session.step(75, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(76, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"],["On Click","Select"],["Row Source","Selected"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"],["On Click","Select"],["Row Source","Selected"]]);
    await session.step(82, "And no rows should be selected", () => noneSelected(page));
    await session.step(83, "And trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
    await session.step(84, "When user selects all rows", () => selectAllRows(page));
    await session.step(85, "Then trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
    await session.step(86, "And the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(87, "When user clears the row selection", () => clearSelection(page));
    await session.step(88, "Then trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
    await session.step(89, "When user sets \"Row Source\" property of trellis plot viewer to \"All\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "All"));
    await session.step(90, "Then trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
    await session.step(91, "And the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
  });
  test("With nothing selected a Selected trellis comes back from a project alive", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX, CONTROL"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(24, "Then the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(95, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["On Click","Select"],["Row Source","Selected"]]), [["On Click","Select"],["Row Source","Selected"]]);
    await session.step(98, "Then no rows should be selected", () => noneSelected(page));
    await session.step(99, "And trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
    await session.step(100, "When user saves the current view as project \"bdd-trellis-empty\"", () => saveAsProject(page, "bdd-trellis-empty"));
    await session.step(101, "And user closes all views", () => closeAllViews(page));
    await session.step(102, "And user opens the \"bdd-trellis-empty\" project", () => openProject(page, "bdd-trellis-empty"));
    await session.step(103, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
    await session.step(104, "And \"Row Source\" property of trellis plot viewer should be \"Selected\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Selected"));
    await session.step(105, "And no rows should be selected", () => noneSelected(page));
    await session.step(106, "And trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
    await session.step(107, "When user selects all rows", () => selectAllRows(page));
    await session.step(108, "Then trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
    await session.step(109, "And the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
    await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(111, "And no errors should have been logged", () => noErrors(page));
  });
});
