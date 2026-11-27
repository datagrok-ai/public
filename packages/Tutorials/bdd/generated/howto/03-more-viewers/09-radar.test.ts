/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/09-radar.feature
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
import {check, clickOn, expand, hoverOver, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/09-radar.feature", import.meta.url);
  test("Compare rows across many columns with a radar", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "When the video says \"A radar shows a row as a shape across several columns\"", () => videoSays(page, "A radar shows a row as a shape across several columns"));
    await session.step(21, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(22, "And user types \"radar\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "radar", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(23, "And user clicks on first \"Radar\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Radar\" card in \"Add Viewer\" dialog")));
    await session.step(24, "Then radar viewer should be visible", () => shouldBe(page, el("radar viewer"), "visible"));
    await session.step(26, "When user hovers over radar viewer", () => hoverOver(page, el("radar viewer")), undefined, "Pick the columns for the axes in the settings");
    await session.step(27, "And user clicks on settings icon of radar viewer", () => clickOn(page, el("settings icon of radar viewer")));
    await session.step(28, "And user expands \"Value\" category in context panel", () => expand(page, el("\"Value\" category in context panel")));
    await session.step(29, "And user clicks on \"...\" button in \"Values\" property in context panel", () => clickOn(page, el("\"...\" button in \"Values\" property in context panel")));
    await session.step(30, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(31, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(32, "And user toggles the \"curb.weight\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "curb.weight", el("\"Select columns...\" dialog")));
    await session.step(33, "And user toggles the \"eng.size\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "eng.size", el("\"Select columns...\" dialog")));
    await session.step(34, "And user toggles the \"horsepower\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "horsepower", el("\"Select columns...\" dialog")));
    await session.step(35, "And user toggles the \"city.mpg\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "city.mpg", el("\"Select columns...\" dialog")));
    await session.step(36, "And user toggles the \"price\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "price", el("\"Select columns...\" dialog")));
    await session.step(37, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(38, "Then \"Values\" property of radar viewer should be \"curb.weight, eng.size, horsepower, city.mpg, price\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "curb.weight, eng.size, horsepower, city.mpg, price"));
    await session.step(40, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: write the values of the current row on the axes");
    await session.step(41, "And user checks \"Show Values\" property in context panel", () => check(page, el("\"Show Values\" property in context panel")));
    await session.step(42, "When the video says \"The shaded band is where most cars fall\" at the top", () => videoSaysAt(page, "The shaded band is where most cars fall", "top"));
    await session.step(44, "And user clicks on the \"cell 19 of model\" area of grid", () => clickArea(page, "cell 19 of model", el("grid")), undefined, "Click a car to draw its shape over the rest");
    await session.step(45, "Then row 19 should be current", () => currentRowIs(page, 19));
    await session.step(47, "When user clicks on the \"cell 28 of model\" area of grid", () => clickArea(page, "cell 28 of model", el("grid")), undefined, "Another car, another shape");
    await session.step(48, "Then row 28 should be current", () => currentRowIs(page, 28));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
