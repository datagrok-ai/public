/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/01-add-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dockedAtViewEdge} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/01-add-viewers.feature", import.meta.url);
  test("Add viewers to a table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When the video says \"Viewers are charts of the open table\"", () => videoSays(page, "Viewers are charts of the open table"));
    await session.step(22, "And user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "The toolbox: the common viewers, one click each");
    await session.step(23, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(24, "When user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")));
    await session.step(25, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(26, "When user clicks on \"histogram\" icon in toolbox", () => clickOn(page, el("\"histogram\" icon in toolbox")));
    await session.step(27, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(29, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")), undefined, "Add viewer opens the gallery of every viewer");
    await session.step(30, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(32, "When user clicks on \"Trend\" viewer tag", () => clickOn(page, el("\"Trend\" viewer tag")), undefined, "Tags group viewers by what they show");
    await session.step(33, "And user clicks on \"Line chart\" viewer card", () => clickOn(page, el("\"Line chart\" viewer card")));
    await session.step(34, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(35, "And line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
    await session.step(36, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(38, "And user types \"bar\" into viewer gallery search", () => typeInto(page, "bar", el("viewer gallery search")), undefined, "Or search by name");
    await session.step(39, "And user clicks on \"Bar chart\" viewer card", () => clickOn(page, el("\"Bar chart\" viewer card")));
    await session.step(40, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(42, "When user docks bar chart viewer to the bottom edge of the view", () => dockToViewEdge(page, el("bar chart viewer"), "bottom"), undefined, "Drag a viewer by its title to dock it elsewhere");
    await session.step(43, "Then bar chart viewer should be docked along the bottom edge of the view", () => dockedAtViewEdge(page, el("bar chart viewer"), "bottom"));
    await session.step(45, "When user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")), undefined, "The x on the title bar closes a viewer");
    await session.step(46, "And user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
    await session.step(47, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
    await session.step(48, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
