/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-row-source.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cellsWideTall} from '../../../bindings/trellis-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, filterPassesAll, noneSelected, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, propertyShouldBe, readingIs, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot row source", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-row-source.feature", import.meta.url);
  test("Trellis plot row source", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Pack Categories","false"]]));
    await session.step(20, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(21, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
    await run.scenario("On Click Filter moves Row Source to All, and Filtered moves On Click back to None", async () => {
      await session.step(24, "Then \"Row Source\" property of trellis plot viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(25, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(26, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(27, "Then \"Row Source\" property of trellis plot viewer should be \"All\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(28, "And \"On Click\" property of trellis plot viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(29, "When user sets \"Row Source\" property of trellis plot viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(30, "Then \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(31, "And \"Row Source\" property of trellis plot viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter and the selection every rung is read against", async () => {
      await session.step(35, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(36, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(37, "When user selects rows where \"RACE\" is \"Caucasian\"", () => selectWhereIs(page, "RACE", "Caucasian"));
      await session.step(38, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(39, "When user sets \"Row Source\" property of trellis plot viewer to \"All\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(40, "Then trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(41, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(42, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source All feeds the cells 1000 rows [source=All, shown=1000, blank=0, pictures=8]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"All\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"All\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(48, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 8", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Filtered feeds the cells 553 rows [source=Filtered, shown=553, blank=4, pictures=5]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(48, "And trellis plot viewer should show 553 rows", () => showsRows(page, el("trellis plot viewer"), 553));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 4", () => readingIs(page, "blank cells", el("trellis plot viewer"), 4));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 5", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 5));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected feeds the cells 896 rows [source=Selected, shown=896, blank=6, pictures=3]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"Selected\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Selected"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"Selected\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Selected"));
      await session.step(48, "And trellis plot viewer should show 896 rows", () => showsRows(page, el("trellis plot viewer"), 896));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 6", () => readingIs(page, "blank cells", el("trellis plot viewer"), 6));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 3", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 3));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source FilteredSelected feeds the cells 480 rows [source=FilteredSelected, shown=480, blank=7, pictures=2]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"FilteredSelected\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "FilteredSelected"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"FilteredSelected\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "FilteredSelected"));
      await session.step(48, "And trellis plot viewer should show 480 rows", () => showsRows(page, el("trellis plot viewer"), 480));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 7", () => readingIs(page, "blank cells", el("trellis plot viewer"), 7));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 2", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source SelectedOrCurrent feeds the cells 896 rows [source=SelectedOrCurrent, shown=896, blank=6, pictures=3]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"SelectedOrCurrent\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "SelectedOrCurrent"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"SelectedOrCurrent\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "SelectedOrCurrent"));
      await session.step(48, "And trellis plot viewer should show 896 rows", () => showsRows(page, el("trellis plot viewer"), 896));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 6", () => readingIs(page, "blank cells", el("trellis plot viewer"), 6));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 3", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 3));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source CurrentRow feeds the cells 1 rows [source=CurrentRow, shown=1, blank=7, pictures=2]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"CurrentRow\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "CurrentRow"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"CurrentRow\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "CurrentRow"));
      await session.step(48, "And trellis plot viewer should show 1 rows", () => showsRows(page, el("trellis plot viewer"), 1));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 7", () => readingIs(page, "blank cells", el("trellis plot viewer"), 7));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 2", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source MouseOverGroup feeds the cells 0 rows [source=MouseOverGroup, shown=0, blank=8, pictures=1]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"MouseOverGroup\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "MouseOverGroup"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"MouseOverGroup\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "MouseOverGroup"));
      await session.step(48, "And trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "blank cells", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 1", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 1));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source MouseOverRow feeds the cells 0 rows [source=MouseOverRow, shown=0, blank=8, pictures=1]", async () => {
      await session.step(46, "When user sets \"Row Source\" property of trellis plot viewer to \"MouseOverRow\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "MouseOverRow"));
      await session.step(47, "Then \"Row Source\" property of trellis plot viewer should be \"MouseOverRow\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "MouseOverRow"));
      await session.step(48, "And trellis plot viewer should show 0 rows", () => showsRows(page, el("trellis plot viewer"), 0));
      await session.step(49, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(51, "And the \"blank cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "blank cells", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 1", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 1));
      await session.step(53, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Putting the ladder back", async () => {
      await session.step(68, "When user sets \"Row Source\" property of trellis plot viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(69, "And user clears the row selection", () => clearSelection(page));
      await session.step(70, "Then no rows should be selected", () => noneSelected(page));
      await session.step(71, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(72, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(73, "Then \"SEX\" filter card should be absent", () => shouldBe(page, el("\"SEX\" filter card"), "absent"));
      await session.step(74, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(75, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(76, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
