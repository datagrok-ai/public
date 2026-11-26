/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/11-categorize-numbers.feature
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
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/11-categorize-numbers.feature", import.meta.url);
  test("Categorize numbers into ranges", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(15, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "When user picks \"Data > Categorize...\" from the top menu", () => pickFromTopMenu(page, "Data > Categorize..."), undefined, "Categorize: put numbers into named ranges");
    await session.step(20, "Then Categorize Data dialog should be visible", () => shouldBe(page, el("Categorize Data dialog"), "visible"));
    await session.step(21, "When user selects \"AGE\" in \"Categorize\" input in Categorize Data dialog", () => selectIn(page, "AGE", el("\"Categorize\" input in Categorize Data dialog")));
    await session.step(23, "And user clicks on \"Bin...\" link in Categorize Data dialog", () => clickOn(page, el("\"Bin...\" link in Categorize Data dialog")), undefined, "Bin: ranges of equal width");
    await session.step(24, "And user enters \"4\" into \"Number of bins\" input in Edit dialog", () => enterInto(page, "4", el("\"Number of bins\" input in Edit dialog")));
    await session.step(25, "And user clicks on OK button in Edit dialog", () => clickOn(page, el("OK button in Edit dialog")));
    await session.step(26, "Then Edit dialog should be hidden", () => shouldBe(page, el("Edit dialog"), "hidden"));
    await session.step(28, "And Categorize Data dialog should contain text \"Matching: 5849\"", () => shouldContainText(page, el("Categorize Data dialog"), "Matching: 5849"));
    await session.step(30, "When user clicks on \"SAVE COLUMN\" button in Categorize Data dialog", () => clickOn(page, el("\"SAVE COLUMN\" button in Categorize Data dialog")), undefined, "SAVE COLUMN adds each row's range as a new column");
    await session.step(31, "Then the table should have a column \"category\"", () => hasColumn(page, "category"));
    await session.step(32, "And \"category\" column should be color-coded categorically", () => colorCodedCategorically(page, "category"));
    await session.step(33, "When user clicks on \"Close\" icon in Categorize Data dialog", () => clickOn(page, el("\"Close\" icon in Categorize Data dialog")));
    await session.step(34, "Then Categorize Data dialog should be hidden", () => shouldBe(page, el("Categorize Data dialog"), "hidden"));
    await session.step(35, "When the video says \"Each row's age range, color-coded\" at the bottom", () => videoSaysAt(page, "Each row's age range, color-coded", "bottom"));
    await session.step(36, "And user hovers over the \"header category\" area of grid", () => hoverArea(page, "header category", el("grid")));
    await session.step(37, "Then no errors should have been logged", () => noErrors(page));
  });
});
