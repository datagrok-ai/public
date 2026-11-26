/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/01-box-plot.feature
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
import {check, clickOn, collapse, expand, hoverOver, selectIn, shouldBe, shouldContainText, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hoverArea, noErrors, propertyShouldBe, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/01-box-plot.feature", import.meta.url);
  test("Compare groups with a box plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(19, "When the video says \"A box plot compares a number across groups\"", () => videoSays(page, "A box plot compares a number across groups"));
    await session.step(20, "And user clicks on box plot icon on toolbox", () => clickOn(page, el("box plot icon on toolbox")));
    await session.step(21, "Then box plot viewer should be added to the open tableview", () => viewerAdded(page, "box plot"));
    await session.step(23, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(24, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(26, "When user picks \"HEIGHT\" in the \"Value\" column selector of box plot viewer", () => pickInColumnSelector(page, "HEIGHT", "Value", el("box plot viewer")), undefined, "Value: the number to compare");
    await session.step(28, "And user picks \"RACE\" in the \"Category 1\" column selector of box plot viewer", () => pickInColumnSelector(page, "RACE", "Category 1", el("box plot viewer")), undefined, "Category: the groups");
    await session.step(29, "When the video says \"Under the plot: min, max, average and median of each group\" at the top", () => videoSaysAt(page, "Under the plot: min, max, average and median of each group", "top"));
    await session.step(30, "And user hovers over the \"stats\" area of box plot viewer", () => hoverArea(page, "stats", el("box plot viewer")));
    await session.step(32, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")), undefined, "The p-value tells whether the groups really differ");
    await session.step(33, "Then tooltip should contain text \"Alexander and Govern\"", () => shouldContainText(page, el("tooltip"), "Alexander and Govern"));
    await session.step(35, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")), undefined, "Everything else is in the settings");
    await session.step(36, "And user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
    await session.step(37, "Then the context panel should show \"Box plot\"", () => contextPanelShows(page, "Box plot"));
    await session.step(39, "When user selects \"SEX\" in \"Category2\" property in context panel", () => selectIn(page, "SEX", el("\"Category2\" property in context panel")), undefined, "Data: split each group again by a second category");
    await session.step(40, "Then box plot viewer should have a \"category Asian | F\" area", () => hasArea(page, el("box plot viewer"), "category Asian | F"));
    await session.step(42, "When user expands \"Statistics\" category in context panel", () => expand(page, el("\"Statistics\" category in context panel")), undefined, "Statistics: choose what the table under the plot shows");
    await session.step(43, "And user checks \"Show Q1\" property in context panel", () => check(page, el("\"Show Q1\" property in context panel")));
    await session.step(44, "And user checks \"Show Q3\" property in context panel", () => check(page, el("\"Show Q3\" property in context panel")));
    await session.step(45, "Then box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
    await session.step(46, "When user collapses \"Statistics\" category in context panel", () => collapse(page, el("\"Statistics\" category in context panel")));
    await session.step(48, "When user expands \"Color\" category in context panel", () => expand(page, el("\"Color\" category in context panel")), undefined, "Color: color the dots by any column");
    await session.step(49, "And user selects \"AGE\" in \"Marker Color\" property in context panel", () => selectIn(page, "AGE", el("\"Marker Color\" property in context panel")));
    await session.step(50, "Then box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
    await session.step(52, "When user expands \"Style\" category in context panel", () => expand(page, el("\"Style\" category in context panel")), undefined, "Style: a violin shows the shape of each distribution");
    await session.step(53, "And user selects \"violin\" in \"Plot Style\" property in context panel", () => selectIn(page, "violin", el("\"Plot Style\" property in context panel")));
    await session.step(54, "Then \"Plot Style\" property of box plot viewer should be \"violin\"", () => propertyShouldBe(page, "Plot Style", el("box plot viewer"), "violin"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
});
