/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/02-scatter-plot.feature
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
import {clickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, paintedInColors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/02-scatter-plot.feature", import.meta.url);
  test("Add a scatter plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(14, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(15, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(16, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")), undefined, "Add viewer opens the gallery of all viewers");
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(21, "When user types \"scatter\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "scatter", el("viewer gallery search in \"Add Viewer\" dialog")), undefined, "Search the gallery, then click the viewer");
    await session.step(22, "And user clicks on first \"Scatter plot\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Scatter plot\" button in \"Add Viewer\" dialog")));
    await session.step(23, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(24, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(26, "When user picks \"AGE\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "AGE", "x", el("scatter plot viewer")), undefined, "Pick the columns right on the plot");
    await session.step(27, "And user picks \"WEIGHT\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("scatter plot viewer")));
    await session.step(29, "And user picks \"SEX\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "SEX", "color", el("scatter plot viewer")), undefined, "Color the points by any column");
    await session.step(30, "Then scatter plot viewer should be painted in at least 2 colors", () => paintedInColors(page, el("scatter plot viewer"), 2));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
});
