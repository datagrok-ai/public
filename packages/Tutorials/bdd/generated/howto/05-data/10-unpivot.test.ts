/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/10-unpivot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, displayedInRow, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/10-unpivot.feature", import.meta.url);
  test("Unpivot columns into rows", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "When the video says \"Four measurements per flower, each in its own column\"", () => videoSays(page, "Four measurements per flower, each in its own column"));
    await session.step(21, "And user picks \"Data > Unpivot...\" from the top menu", () => pickFromTopMenu(page, "Data > Unpivot..."), undefined, "Unpivot: turn columns into rows");
    await session.step(22, "Then Unpivot dialog should be visible", () => shouldBe(page, el("Unpivot dialog"), "visible"));
    await session.step(24, "When user enters \"Measure\" into \"Category column name\" input in Unpivot dialog", () => enterInto(page, "Measure", el("\"Category column name\" input in Unpivot dialog")), undefined, "Name the column that will hold the measurement names");
    await session.step(26, "And user selects \"Copy\" in 3rd input in Unpivot dialog", () => selectIn(page, "Copy", el("3rd input in Unpivot dialog")), undefined, "Copy: keep the column as it is in every new row");
    await session.step(27, "And user selects \"Copy\" in 8th input in Unpivot dialog", () => selectIn(page, "Copy", el("8th input in Unpivot dialog")));
    await session.step(29, "And user selects \"Merge\" in 4th input in Unpivot dialog", () => selectIn(page, "Merge", el("4th input in Unpivot dialog")), undefined, "Merge: fold the columns into Measure and Value");
    await session.step(30, "And user selects \"Merge\" in 5th input in Unpivot dialog", () => selectIn(page, "Merge", el("5th input in Unpivot dialog")));
    await session.step(31, "And user selects \"Merge\" in 6th input in Unpivot dialog", () => selectIn(page, "Merge", el("6th input in Unpivot dialog")));
    await session.step(32, "And user selects \"Merge\" in 7th input in Unpivot dialog", () => selectIn(page, "Merge", el("7th input in Unpivot dialog")));
    await session.step(33, "And user clicks on OK button in Unpivot dialog", () => clickOn(page, el("OK button in Unpivot dialog")));
    await session.step(34, "Then Unpivot dialog should be hidden", () => shouldBe(page, el("Unpivot dialog"), "hidden"));
    await session.step(35, "And the table should have 600 rows", () => rowCount(page, 600));
    await session.step(36, "And the table should have 4 columns", () => columnCount(page, 4));
    await session.step(37, "And the table should have a column \"Measure\"", () => hasColumn(page, "Measure"));
    await session.step(38, "And the value of \"Measure\" column in row 1 should be \"Sepal.Length\"", () => valueInRow(page, "Measure", 1, "Sepal.Length"));
    await session.step(39, "And the \"Value\" cell of row 1 should be displayed as \"5.10\"", () => displayedInRow(page, "Value", 1, "5.10"));
    await session.step(40, "When the video says \"One row per flower and measurement: 150 × 4 = 600\" at the bottom", () => videoSaysAt(page, "One row per flower and measurement: 150 × 4 = 600", "bottom"));
    await session.step(41, "And user hovers over the \"cell 1 of Measure\" area of grid", () => hoverArea(page, "cell 1 of Measure", el("grid")));
    await session.step(42, "Then no errors should have been logged", () => noErrors(page));
  });
});
