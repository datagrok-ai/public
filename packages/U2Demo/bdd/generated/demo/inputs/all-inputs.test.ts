/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/inputs/all-inputs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.property-form]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {asTallAs, leftAligned, noLayoutProblems} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {check, enterInto, selectIn, shouldBe, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("All inputs", () => {
  const session = feature(test, "features/demo/inputs/all-inputs.feature", import.meta.url);
  test("All inputs", {tag: ["@demo", "@journey", "@realizes:u2.dg.property-form"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(8, "Given user opens the \"All inputs\" demo page", () => openDemoPage(page, "All inputs"));
    enter(page, "U2 Demo");
    await run.scenario("The u2 editors stand in one column at the platform's row height", async () => {
      await session.step(11, "Then Name input in \"u2\" section should be as tall as Name input in \"Platform\" section", () => asTallAs(page, el("Name input in \"u2\" section"), el("Name input in \"Platform\" section")));
      await session.step(12, "And Stage input in \"u2\" section should be as tall as Stage input in \"Platform\" section", () => asTallAs(page, el("Stage input in \"u2\" section"), el("Stage input in \"Platform\" section")));
      await session.step(13, "And Replicates input in \"u2\" section should be as tall as Replicates input in \"Platform\" section", () => asTallAs(page, el("Replicates input in \"u2\" section"), el("Replicates input in \"Platform\" section")));
      await session.step(14, "And editor of Name input in \"u2\" section should be left-aligned with editor of Stage input in \"u2\" section", () => leftAligned(page, el("editor of Name input in \"u2\" section"), el("editor of Stage input in \"u2\" section")));
      await session.step(15, "And editor of Stage input in \"u2\" section should be left-aligned with editor of Replicates input in \"u2\" section", () => leftAligned(page, el("editor of Stage input in \"u2\" section"), el("editor of Replicates input in \"u2\" section")));
      await session.step(16, "And editor of Replicates input in \"u2\" section should be left-aligned with editor of Active checkbox in \"u2\" section", () => leftAligned(page, el("editor of Replicates input in \"u2\" section"), el("editor of Active checkbox in \"u2\" section")));
    });
    await run.scenario("The u2 checkbox row is as tall as the platform's", async () => {
      await session.step(19, "Then Active checkbox in \"u2\" section should be as tall as Active checkbox in \"Platform\" section", () => asTallAs(page, el("Active checkbox in \"u2\" section"), el("Active checkbox in \"Platform\" section")));
    });
    await run.scenario("The page shows no layout problems", async () => {
      await session.step(22, "Then demo page should have no layout problems", () => noLayoutProblems(page, el("demo page")));
    });
    await run.scenario("A name typed on either side shows on the other", async () => {
      await session.step(25, "When user enters \"Ibuprofen\" into Name input in \"u2\" section", () => enterInto(page, "Ibuprofen", el("Name input in \"u2\" section")));
      await session.step(26, "Then Name input in \"Platform\" section should have value \"Ibuprofen\"", () => shouldHaveValue(page, el("Name input in \"Platform\" section"), "Ibuprofen"));
      await session.step(27, "When user enters \"Aspirin\" into Name input in \"Platform\" section", () => enterInto(page, "Aspirin", el("Name input in \"Platform\" section")));
      await session.step(28, "Then Name input in \"u2\" section should have value \"Aspirin\"", () => shouldHaveValue(page, el("Name input in \"u2\" section"), "Aspirin"));
    });
    await run.scenario("A stage picked on either side shows on the other", async () => {
      await session.step(31, "When user selects \"Phase II\" in Stage input in \"u2\" section", () => selectIn(page, "Phase II", el("Stage input in \"u2\" section")));
      await session.step(32, "Then Stage input in \"Platform\" section should have value \"Phase II\"", () => shouldHaveValue(page, el("Stage input in \"Platform\" section"), "Phase II"));
      await session.step(33, "When user selects \"Phase I\" in Stage input in \"Platform\" section", () => selectIn(page, "Phase I", el("Stage input in \"Platform\" section")));
      await session.step(34, "Then Stage input in \"u2\" section should have value \"Phase I\"", () => shouldHaveValue(page, el("Stage input in \"u2\" section"), "Phase I"));
    });
    await run.scenario("Replicates entered on either side show on the other", async () => {
      await session.step(37, "When user enters \"7\" into Replicates input in \"u2\" section", () => enterInto(page, "7", el("Replicates input in \"u2\" section")));
      await session.step(38, "Then Replicates input in \"Platform\" section should have value \"7\"", () => shouldHaveValue(page, el("Replicates input in \"Platform\" section"), "7"));
      await session.step(39, "When user enters \"3\" into Replicates input in \"Platform\" section", () => enterInto(page, "3", el("Replicates input in \"Platform\" section")));
      await session.step(40, "Then Replicates input in \"u2\" section should have value \"3\"", () => shouldHaveValue(page, el("Replicates input in \"u2\" section"), "3"));
    });
    await run.scenario("Active toggled on either side shows on the other", async () => {
      await session.step(43, "When user unchecks Active checkbox in \"u2\" section", () => uncheck(page, el("Active checkbox in \"u2\" section")));
      await session.step(44, "Then Active checkbox in \"Platform\" section should be unchecked", () => shouldBe(page, el("Active checkbox in \"Platform\" section"), "unchecked"));
      await session.step(45, "When user checks Active checkbox in \"Platform\" section", () => check(page, el("Active checkbox in \"Platform\" section")));
      await session.step(46, "Then Active checkbox in \"u2\" section should be checked", () => shouldBe(page, el("Active checkbox in \"u2\" section"), "checked"));
    });
    run.finish();
  });
});
