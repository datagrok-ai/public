/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-summary-five-rgroups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {callWith} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary over five R-group columns", () => {
  const session = feature(test, "features/guides/sar-summary-five-rgroups.feature", import.meta.url);
  test("six components on four thousand compounds", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And user opens five-rgroups dataset", () => openDataset(page, ds("five-rgroups")));
    await session.step(12, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(13, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:Compound"],["activity","column:pIC50"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Scaffold"],["fragmentColumns","columns:R1,R2,R3,R4,R5"],["columnAxis","R2"]]), [["table","table"],["molecules","column:Compound"],["activity","column:pIC50"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Scaffold"],["fragmentColumns","columns:R1,R2,R3,R4,R5"],["columnAxis","R2"]]);
    await session.step(27, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(28, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(29, "And user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(30, "Then tab panel should contain text \"Changing R2 moves pIC50 most: R2 > Scaffold > R1 > R5\"", () => shouldContainText(page, el("tab panel"), "Changing R2 moves pIC50 most: R2 > Scaffold > R1 > R5"));
    await session.step(31, "And tab panel should contain text \"R2 — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "R2 — offsets from the additive fit"));
    await session.step(32, "When user clicks on \"R5\" effects tab", () => clickOn(page, el("\"R5\" effects tab")));
    await session.step(33, "Then tab panel should contain text \"R5 — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "R5 — offsets from the additive fit"));
    await session.step(34, "When user clicks on \"R3\" effects tab", () => clickOn(page, el("\"R3\" effects tab")));
    await session.step(35, "Then tab panel should contain text \"nothing at this position\"", () => shouldContainText(page, el("tab panel"), "nothing at this position"));
    await session.step(36, "When user clicks on \"Method\" summary segment", () => clickOn(page, el("\"Method\" summary segment")));
    await session.step(37, "Then tab panel should contain text \"Fit quality\"", () => shouldContainText(page, el("tab panel"), "Fit quality"));
    await session.step(38, "And tab panel should contain text \"Where it holds best\"", () => shouldContainText(page, el("tab panel"), "Where it holds best"));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
