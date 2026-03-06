/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {renderReactionToCanvas} from '../../rendering/rdkit-reaction-renderer';
import {NamedReaction, ReactionMode} from './types';
import {getReactionCategories} from './consts';
import {validateReactionSmarts, saveUserReaction, generateReactionId, loadAllReactions} from './reaction-storage';
import {runTransformationReaction} from './reactions';

/**
 * Opens a dialog for creating or editing a named reaction.
 * Supports pasting SMARTS directly or using molecule sketchers.
 */
export async function openReactionEditor(
  rdkit: RDModule,
  options: {
    /** Existing reaction to edit (null = create new) */
    preset?: NamedReaction;
    /** Default mode to pre-select and lock (e.g. when opened from two-component dialog) */
    defaultMode?: ReactionMode;
    /** Called after a reaction is saved */
    onSave?: (reaction: NamedReaction) => void;
    /** Called when the dialog is cancelled */
    onCancel?: () => void;
  } = {},
): Promise<void> {
  const {preset, defaultMode, onSave, onCancel} = options;
  const isEdit = !!preset;

  // ---- State ----
  let currentSmarts = preset?.reactionSmarts ?? '';
  let validationState = {valid: false, error: 'Enter a reaction SMARTS', numReactants: 0, numProducts: 0};
  let _okBtn: HTMLButtonElement | null = null;

  // ---- Inputs ----
  const nameInput = ui.input.string('Reaction Name', {
    value: preset?.name ?? '',
    nullable: false,
  });
  ui.tooltip.bind(nameInput.input, 'A descriptive name for this reaction, e.g. "Suzuki Coupling".');

  const allReactions = await loadAllReactions();
  const existingCategories = getReactionCategories(allReactions);
  const categoryInput = ui.input.choice('Category', {
    value: preset?.category ?? existingCategories[0] ?? 'Uncategorized',
    items: existingCategories,
    nullable: false,
  });
  ui.tooltip.bind(categoryInput.input, 'Grouping category. You can also type a new category name.');
  categoryInput.addOptions(ui.icons.add(() => {
    const catInput = ui.input.string('Category name', {value: '', nullable: false});
    ui.dialog('New Category')
      .add(catInput)
      .onOK(() => {
        const newCat = catInput.value?.trim();
        if (newCat) {
          existingCategories.push(newCat);
          (categoryInput as DG.ChoiceInput<string>).items = existingCategories;
          categoryInput.value = newCat;
        }
      })
      .show();
  }, 'Add new category'));

  const descriptionInput = ui.input.string('Description', {
    value: preset?.description ?? '',
    nullable: true,
  });
  ui.tooltip.bind(descriptionInput.input, 'A short description shown in tooltips when hovering this reaction.');

  const modeValue = preset?.mode ?? defaultMode ?? 'transformation';
  const modeInput = ui.input.choice('Mode', {
    value: modeValue,
    items: ['transformation', 'two-component'] as ReactionMode[],
    nullable: false,
  });
  ui.tooltip.bind(modeInput.input, 'Transformation: single column. Two-component: two columns of molecules.');
  if (defaultMode) {
    modeInput.enabled = false;
    ui.tooltip.bind(modeInput.input, `Mode is locked to "${defaultMode}".`);
  }

  const smartsInput = ui.input.textArea('Reaction SMARTS', {
    value: currentSmarts,
    nullable: false,
  });
  ui.tooltip.bind(smartsInput.input, 'Enter a valid reaction SMIRKS/SMARTS string (reactants>>products).');
  smartsInput.root.style.minHeight = '60px';

  // ---- Preview canvas ----
  const previewWidth = 500;
  const previewHeight = 130;
  const r = window.devicePixelRatio;
  const previewCanvas = ui.canvas(previewWidth * r, previewHeight * r);
  previewCanvas.style.width = `${previewWidth}px`;
  previewCanvas.style.height = `${previewHeight}px`;
  previewCanvas.style.border = '1px solid var(--grey-2)';
  previewCanvas.style.borderRadius = '4px';
  previewCanvas.style.marginTop = '8px';

  // ---- Validation label ----
  const validationDiv = ui.div([], {style: {padding: '4px 0', fontSize: '12px'}});

  function updateValidation() {
    const smarts = smartsInput.value?.trim() ?? '';
    currentSmarts = smarts;
    const ctx = previewCanvas.getContext('2d')!;
    ctx.clearRect(0, 0, previewCanvas.width, previewCanvas.height);

    if (!smarts) {
      validationState = {valid: false, error: 'Reaction SMARTS is required.', numReactants: 0, numProducts: 0};
      renderValidation();
      return;
    }

    validationState = validateReactionSmarts(rdkit, smarts) as typeof validationState;
    if (validationState.valid) {
      try {
        renderReactionToCanvas(rdkit, previewCanvas, smarts, previewCanvas.width, previewCanvas.height);
      } catch {/* failed to render preview */}
      // auto-detect mode (only when not locked)
      if (!defaultMode) {
        if (validationState.numReactants === 1)
          modeInput.value = 'transformation';
        else if (validationState.numReactants === 2)
          modeInput.value = 'two-component';
      }
    } else {
      ctx.save();
      ctx.scale(r, r);
      ctx.fillStyle = '#cc0000';
      ctx.font = '13px Roboto, sans-serif';
      ctx.fillText(validationState.error ?? 'Invalid reaction', 12, previewHeight / 2);
      ctx.restore();
    }
    renderValidation();
  }

  function renderValidation() {
    validationDiv.innerHTML = '';
    if (validationState.valid) {
      validationDiv.style.color = 'var(--green-1)';
      validationDiv.textContent = `Valid  |  Reactants: ${validationState.numReactants}  |  Products: ${validationState.numProducts}`;
    } else {
      validationDiv.style.color = '#cc0000';
      validationDiv.textContent = validationState.error ?? 'Invalid';
    }
    _okBtn?.classList.toggle('d4-disabled', !canSave());
  }

  smartsInput.onInput.subscribe(() => updateValidation());
  nameInput.onInput.subscribe(() => _okBtn?.classList.toggle('d4-disabled', !canSave()));
  categoryInput.onChanged.subscribe(() => _okBtn?.classList.toggle('d4-disabled', !canSave()));

  // ---- Test run section ----
  const testSmilesInput = ui.input.string('Test SMILES', {
    value: preset?.exampleReactants?.[0] ?? '',
    nullable: true,
  });
  ui.tooltip.bind(testSmilesInput.input, 'Paste a molecule SMILES to test the reaction.');

  const testResultDiv = ui.div([], {style: {minHeight: '80px', border: '1px solid var(--grey-2)', borderRadius: '4px', padding: '4px', marginTop: '4px'}});

  const testButton = ui.button('Run Test', async () => {
    const smarts = smartsInput.value?.trim();
    const testSmiles = testSmilesInput.value?.trim();
    if (!smarts || !testSmiles) {
      grok.shell.warning('Enter both SMARTS and a test molecule SMILES.');
      return;
    }
    try {
      const result = await runTransformationReaction(rdkit, smarts, [testSmiles], {removeSaltsAndWater: false});
      testResultDiv.innerHTML = '';
      if (result.products[0] && result.products[0] !== testSmiles) {
        const mol = grok.chem.drawMolecule(result.products[0], 200, 80);
        testResultDiv.append(ui.divH([
          ui.divV([ui.label('Product:'), mol]),
          ui.divV([ui.label('SMILES:'), ui.div([ui.element('code').append(document.createTextNode(result.products[0]))])]),
        ], {style: {gap: '12px', alignItems: 'center'}}));
      } else
        testResultDiv.append(ui.divText('No product generated.'));
    } catch (e: any) {
      testResultDiv.innerHTML = '';
      testResultDiv.append(ui.divText(`Error: ${e?.message ?? e}`, {style: {color: '#cc0000'}}));
    }
  }, 'Run the reaction on the test molecule to preview the product');

  // ---- Tags input ----
  const tagsInput = ui.input.string('Tags', {
    value: preset?.tags?.join(', ') ?? '',
    nullable: true,
  });
  ui.tooltip.bind(tagsInput.input, 'Comma-separated tags for searching (e.g. "palladium, cross-coupling").');

  function canSave(): boolean {
    return validationState.valid && !!nameInput.value?.trim() && !!categoryInput.value;
  }

  // ---- Layout ----
  const inputsSection = ui.divV([
    nameInput,
    categoryInput,
    descriptionInput,
    modeInput,
    tagsInput,
  ], {style: {paddingBottom: '8px'}});

  const smartsSection = ui.divV([
    ui.h3('Reaction SMARTS'),
    smartsInput,
    previewCanvas,
    validationDiv,
  ]);

  const testSection = ui.divV([
    ui.h3('Test Run'),
    ui.divH([testSmilesInput.root, testButton], {style: {alignItems: 'flex-end', gap: '8px'}}),
    testResultDiv,
  ]);

  const tabControl = ui.tabControl({
    'Reaction': ui.divV([inputsSection, smartsSection]),
    'Test': testSection,
  }, false);
  tabControl.root.style.minWidth = '550px';
  tabControl.root.style.minHeight = '500px';

  const dialog = ui.dialog(isEdit ? 'Edit Reaction' : 'New Reaction')
    .add(tabControl.root);

  // Use dialog's built-in OK/CANCEL buttons — rename OK to SAVE
  dialog.onOK(async () => {
    if (!canSave()) {
      grok.shell.warning('Please fix validation errors before saving.');
      return;
    }

    const reaction: NamedReaction = {
      id: preset?.id ?? generateReactionId(nameInput.value),
      name: nameInput.value.trim(),
      reactionSmarts: currentSmarts,
      mode: modeInput.value! as ReactionMode,
      category: categoryInput.value!,
      description: descriptionInput.value?.trim() || undefined,
      isUserDefined: true,
      author: DG.User.current().friendlyName,
      tags: tagsInput.value ? tagsInput.value.split(',').map((t) => t.trim()).filter((t) => t.length > 0) : undefined,
      exampleReactants: testSmilesInput.value ? [testSmilesInput.value.trim()] : undefined,
    };

    try {
      await saveUserReaction(reaction);
      grok.shell.info(`Reaction "${reaction.name}" saved.`);
      onSave?.(reaction);
    } catch (e: any) {
      grok.shell.error(`Failed to save: ${e?.message ?? e}`);
    }
  });

  _okBtn = dialog.getButton('OK');
  _okBtn.textContent = 'SAVE';

  // Initial validation (now _okBtn is set, renderValidation will toggle it)
  updateValidation();
  dialog.show({resizable: true});
}
