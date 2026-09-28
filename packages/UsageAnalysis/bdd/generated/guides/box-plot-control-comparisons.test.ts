/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/box-plot-control-comparisons.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, selectIn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hoverArea, openToolbox, pickFromAreaContextMenu, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tell whether groups differ, and compare each group with a control", () => {
  const session = feature(test, "features/guides/box-plot-control-comparisons.feature", import.meta.url);
  test("Test the groups of a box plot, then compare each group with a control group", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "When user opens toolbox", () => openToolbox(page));
    await session.step(19, "And user clicks on box plot icon on toolbox", () => clickOn(page, el("box plot icon on toolbox")));
    await session.step(20, "Then box plot viewer should be added to the open tableview", () => viewerAdded(page, "box plot"));
    await session.step(21, "When user picks \"HEIGHT\" in the \"Value\" column selector of box plot viewer", () => pickInColumnSelector(page, "HEIGHT", "Value", el("box plot viewer")));
    await session.step(22, "And user picks \"RACE\" in the \"Category 1\" column selector of box plot viewer", () => pickInColumnSelector(page, "RACE", "Category 1", el("box plot viewer")));
    await session.step(23, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
    await session.step(24, "Then tooltip should contain text \"Alexander and Govern\"", () => shouldContainText(page, el("tooltip"), "Alexander and Govern"));
    await session.step(25, "When user clicks on show group stats icon in box plot viewer", () => clickOn(page, el("show group stats icon in box plot viewer")));
    await session.step(26, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
    await session.step(27, "And user selects \"Caucasian\" in control group choice input in box plot viewer", () => selectIn(page, "Caucasian", el("control group choice input in box plot viewer")));
    await session.step(28, "Then box plot viewer should have a \"p value of Asian\" area", () => hasArea(page, el("box plot viewer"), "p value of Asian"));
    await session.step(29, "When user picks \"Add Control Comparisons Table\" from the context menu of the \"group comparison\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Add Control Comparisons Table", "group comparison", el("box plot viewer")));
    await session.step(30, "Then table \"Control Comparisons: HEIGHT by RACE vs Caucasian\" should be open", () => tableOpen(page, "Control Comparisons: HEIGHT by RACE vs Caucasian"));
    await session.step(31, "And table \"Control Comparisons: HEIGHT by RACE vs Caucasian\" should have 3 rows", () => tableRows(page, "Control Comparisons: HEIGHT by RACE vs Caucasian", 3));
    await session.step(32, "And the value of \"Conclusion\" column in row 1 should be \"Significant\"", () => valueInRow(page, "Conclusion", 1, "Significant"));
    await session.step(33, "And the value of \"Group\" column in row 2 should be \"Black\"", () => valueInRow(page, "Group", 2, "Black"));
    await session.step(34, "And the value of \"Conclusion\" column in row 2 should be \"Not significant\"", () => valueInRow(page, "Conclusion", 2, "Not significant"));
    await session.step(35, "And the value of \"Conclusion\" column in row 3 should be \"Significant\"", () => valueInRow(page, "Conclusion", 3, "Significant"));
  });
});
