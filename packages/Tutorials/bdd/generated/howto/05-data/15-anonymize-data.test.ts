/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/15-anonymize-data.feature
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
import {clickOn, hoverOver, shouldBe, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/15-anonymize-data.feature", import.meta.url);
  test("Anonymize data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the \"DEMOG\" column of the grid is 90 px wide", () => gridColumnWidth(page, "DEMOG", 90));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When user hovers over the \"cell 1 of USUBJID\" area of grid", () => hoverArea(page, "cell 1 of USUBJID", el("grid")), undefined, "Subject IDs and real ages: not for sharing as they are");
    await session.step(24, "And user picks \"Data > Anonymize...\" from the top menu", () => pickFromTopMenu(page, "Data > Anonymize..."), undefined, "Anonymize: hide real values before you share");
    await session.step(25, "Then Anonymize Data dialog should be visible", () => shouldBe(page, el("Anonymize Data dialog"), "visible"));
    await session.step(27, "When user toggles the \"USUBJID\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "USUBJID", el("Anonymize Data dialog")), undefined, "Tick the columns to anonymize");
    await session.step(28, "And user toggles the \"AGE\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "AGE", el("Anonymize Data dialog")));
    await session.step(29, "And user toggles the \"WEIGHT\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "WEIGHT", el("Anonymize Data dialog")));
    await session.step(31, "And user toggles the \"DEMOG\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "DEMOG", el("Anonymize Data dialog")), undefined, "DEMOG repeats the age: anonymize it too");
    await session.step(33, "And user hovers over \"Number randomization factor\" input in Anonymize Data dialog", () => hoverOver(page, el("\"Number randomization factor\" input in Anonymize Data dialog")), undefined, "Numbers get random noise: 0.1 is up to about 25%");
    await session.step(34, "And user clicks on OK button in Anonymize Data dialog", () => clickOn(page, el("OK button in Anonymize Data dialog")));
    await session.step(35, "Then Anonymize Data dialog should be hidden", () => shouldBe(page, el("Anonymize Data dialog"), "hidden"));
    await session.step(36, "And every value of \"USUBJID\" column should match \"^USUBJID \\d+$\"", () => everyValueMatches(page, "USUBJID", "^USUBJID \\d+$"));
    await session.step(37, "And every value of \"DEMOG\" column should match \"^DEMOG \\d+$\"", () => everyValueMatches(page, "DEMOG", "^DEMOG \\d+$"));
    await session.step(38, "When the video says \"IDs and DEMOG became codes; ages and weights got noise\" at the bottom", () => videoSaysAt(page, "IDs and DEMOG became codes; ages and weights got noise", "bottom"));
    await session.step(39, "And user hovers over the \"cell 1 of AGE\" area of grid", () => hoverArea(page, "cell 1 of AGE", el("grid")));
    await session.step(40, "Then no errors should have been logged", () => noErrors(page));
  });
});
