/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/09-split-a-column.feature
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
import {clickOn, enterInto, hoverOver, selectIn, shouldBe, shouldContainText, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, clickViewTab, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/09-split-a-column.feature", import.meta.url);
  test("Split a column by a separator or a pattern", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens titanic dataset", () => openDataset(page, ds("titanic")));
    await session.step(20, "And the grid shows only the columns \"pclass, survived, name, sex, age\"", () => gridShowsColumns(page, "pclass, survived, name, sex, age"));
    await session.step(21, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "When the video says \"DEMOG holds three values in one text: age, race, sex\"", () => videoSays(page, "DEMOG holds three values in one text: age, race, sex"));
    await session.step(27, "And user clicks on the \"cell 1 of DEMOG\" area of grid", () => clickArea(page, "cell 1 of DEMOG", el("grid")), undefined, "Click the column to split");
    await session.step(29, "And user picks \"Data > Split Column...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column..."), undefined, "Split Column: one text column into several");
    await session.step(30, "Then Split Column dialog should be visible", () => shouldBe(page, el("Split Column dialog"), "visible"));
    await session.step(31, "And Column input in Split Column dialog should contain text \"DEMOG\"", () => shouldContainText(page, el("Column input in Split Column dialog"), "DEMOG"));
    await session.step(33, "When user enters \" \" into \"Split by\" input in Split Column dialog", () => enterInto(page, " ", el("\"Split by\" input in Split Column dialog")), undefined, "Split by: here a space");
    await session.step(34, "And user clicks on OK button in Split Column dialog", () => clickOn(page, el("OK button in Split Column dialog")));
    await session.step(35, "Then Split Column dialog should be hidden", () => shouldBe(page, el("Split Column dialog"), "hidden"));
    await session.step(36, "And the table should have a column \"col3\"", () => hasColumn(page, "col3"));
    await session.step(37, "And the value of \"col1\" column in row 1 should be \"53\"", () => valueInRow(page, "col1", 1, "53"));
    await session.step(38, "And the value of \"col2\" column in row 1 should be \"C\"", () => valueInRow(page, "col2", 1, "C"));
    await session.step(39, "And the value of \"col3\" column in row 1 should be \"F\"", () => valueInRow(page, "col3", 1, "F"));
    await session.step(40, "When the video says \"Three new columns: age, race and sex\" at the bottom", () => videoSaysAt(page, "Three new columns: age, race and sex", "bottom"));
    await session.step(41, "And user hovers over the \"header col1\" area of grid", () => hoverArea(page, "header col1", el("grid")));
    await session.step(42, "When user clicks on the tab of the \"titanic\" view", () => clickViewTab(page, "titanic"));
    await session.step(44, "And user hovers over the \"cell 1 of name\" area of grid", () => hoverArea(page, "cell 1 of name", el("grid")), undefined, "In titanic, each name holds a surname, a title and given names");
    await session.step(46, "And user picks \"Data > Split Column by RegExp...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column by RegExp..."), undefined, "No one separator splits these: split by a pattern");
    await session.step(47, "Then Split Column by RegExp dialog should be visible", () => shouldBe(page, el("Split Column by RegExp dialog"), "visible"));
    await session.step(48, "When user selects \"name\" in \"Column\" input in Split Column by RegExp dialog", () => selectIn(page, "name", el("\"Column\" input in Split Column by RegExp dialog")));
    await session.step(50, "And user enters \"(.+), (\\w+)\\. (.+)\" into \"RegExp\" input in Split Column by RegExp dialog", () => enterInto(page, "(.+), (\\w+)\\. (.+)", el("\"RegExp\" input in Split Column by RegExp dialog")), undefined, "Each (group) in the pattern becomes a column");
    await session.step(52, "And user hovers over \"Preview result columns\" dialog", () => hoverOver(page, el("\"Preview result columns\" dialog")), undefined, "The preview shows the parts it found");
    await session.step(53, "And user clicks on OK button in Split Column by RegExp dialog", () => clickOn(page, el("OK button in Split Column by RegExp dialog")));
    await session.step(54, "Then Split Column by RegExp dialog should be hidden", () => shouldBe(page, el("Split Column by RegExp dialog"), "hidden"));
    await session.step(56, "And the table should have a column \"r (3)\"", () => hasColumn(page, "r (3)"));
    await session.step(57, "And the value of \"r\" column in row 1 should be \"Allen\"", () => valueInRow(page, "r", 1, "Allen"));
    await session.step(58, "And the value of \"r (2)\" column in row 1 should be \"Miss\"", () => valueInRow(page, "r (2)", 1, "Miss"));
    await session.step(59, "And the value of \"r (3)\" column in row 1 should be \"Elisabeth Walton\"", () => valueInRow(page, "r (3)", 1, "Elisabeth Walton"));
    await session.step(60, "When the video says \"Surname, title and given names, each in a column\" at the bottom", () => videoSaysAt(page, "Surname, title and given names, each in a column", "bottom"));
    await session.step(61, "And user hovers over the \"cell 1 of r (2)\" area of grid", () => hoverArea(page, "cell 1 of r (2)", el("grid")));
    await session.step(62, "Then no errors should have been logged", () => noErrors(page));
  });
});
