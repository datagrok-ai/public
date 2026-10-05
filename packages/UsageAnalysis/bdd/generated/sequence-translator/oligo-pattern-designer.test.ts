/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/oligo-pattern-designer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sequencetranslator.app.oligo-toolkit]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, shouldBe, shouldContainText, shouldHaveValue, switchOff} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Oligo Pattern: the default example pattern and the strand switches", () => {
  const session = feature(test, "features/sequence-translator/oligo-pattern-designer.feature", import.meta.url);
  test("The app opens on the default example, and saving over it is refused", {tag: ["@realizes:sequencetranslator.app.oligo-toolkit"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(22, "And Apps---Peptides tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Peptides tree node inside browse tree")));
    await session.step(23, "And Apps---Peptides---Oligo-Toolkit tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Peptides---Oligo-Toolkit tree node inside browse tree")));
    await session.step(24, "When user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree", () => doubleClickOn(page, el("Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree")));
    await session.step(25, "Then the \"Oligo Pattern\" view should be current", () => viewIsCurrent(page, "Oligo Pattern"));
    await session.step(26, "When user hovers over \"Load\" heading", () => hoverOver(page, el("\"Load\" heading")));
    await session.step(27, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(30, "Then Author input should contain the text \"(me)\"", () => shouldContainText(page, el("Author input"), "(me)"));
    await session.step(31, "And Pattern input should have the value \"<default example>\"", () => shouldHaveValue(page, el("Pattern input"), "<default example>"));
    await session.step(32, "And \"Translation example\" heading should be visible", () => shouldBe(page, el("\"Translation example\" heading"), "visible"));
    await session.step(33, "And \"Sense strand\" heading should be visible", () => shouldBe(page, el("\"Sense strand\" heading"), "visible"));
    await session.step(34, "And \"Anti sense\" heading should be visible", () => shouldBe(page, el("\"Anti sense\" heading"), "visible"));
    await session.step(35, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
    await session.step(36, "When user enters \"22\" into \"Sense strand length\" input", () => enterInto(page, "22", el("\"Sense strand length\" input")));
    await session.step(37, "Then Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
    await session.step(38, "When user clicks on Save button", () => clickOn(page, el("Save button")));
    await session.step(39, "Then a warning balloon containing \"Cannot save default pattern\" should have been shown", () => warningBalloonText(page, "Cannot save default pattern"));
    await session.step(40, "And Pattern input should have the value \"<default example>\"", () => shouldHaveValue(page, el("Pattern input"), "<default example>"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("Switching the antisense strand off hides its length and its example", {tag: ["@realizes:sequencetranslator.app.oligo-toolkit"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(22, "And Apps---Peptides tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Peptides tree node inside browse tree")));
    await session.step(23, "And Apps---Peptides---Oligo-Toolkit tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Peptides---Oligo-Toolkit tree node inside browse tree")));
    await session.step(24, "When user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree", () => doubleClickOn(page, el("Apps---Peptides---Oligo-Toolkit---Oligo-Pattern tree node inside browse tree")));
    await session.step(25, "Then the \"Oligo Pattern\" view should be current", () => viewIsCurrent(page, "Oligo Pattern"));
    await session.step(26, "When user hovers over \"Load\" heading", () => hoverOver(page, el("\"Load\" heading")));
    await session.step(27, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(44, "Then \"Anti sense length\" input should be visible", () => shouldBe(page, el("\"Anti sense length\" input"), "visible"));
    await session.step(45, "And \"Anti sense\" heading should be visible", () => shouldBe(page, el("\"Anti sense\" heading"), "visible"));
    await session.step(46, "When user enters \"10\" into \"Sense strand length\" input", () => enterInto(page, "10", el("\"Sense strand length\" input")));
    await session.step(47, "And user switches off \"Anti sense strand\" checkbox", () => switchOff(page, el("\"Anti sense strand\" checkbox")));
    await session.step(48, "Then \"Anti sense length\" input should be hidden", () => shouldBe(page, el("\"Anti sense length\" input"), "hidden"));
    await session.step(49, "And \"Anti sense\" heading should be hidden", () => shouldBe(page, el("\"Anti sense\" heading"), "hidden"));
    await session.step(50, "And \"Sense strand\" heading should be visible", () => shouldBe(page, el("\"Sense strand\" heading"), "visible"));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
