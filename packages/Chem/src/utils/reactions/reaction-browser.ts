/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {renderReactionToCanvas} from '../../rendering/rdkit-reaction-renderer';
import {NamedReaction, ReactionMode} from './types';
import {getReactionCategories} from './consts';
import {loadAllReactions, deleteUserReaction, invalidateReactionsCache} from './reaction-storage';
import {resolveReactionVariables} from './reactions';
import {openReactionEditor} from './reaction-editor';
import {Subject} from 'rxjs';

const CARD_WIDTH = 200;
const CARD_HEIGHT = 160;
const PREVIEW_WIDTH = 190;
const PREVIEW_HEIGHT = 80;

/**
 * A visual browsable catalog of reactions with search, category/mode filters,
 * card previews, and selection events.
 */
export class ReactionBrowser {
  private reactions: NamedReaction[] = [];
  private filteredReactions: NamedReaction[] = [];
  private selectedReaction: NamedReaction | null = null;
  private modeFilter: ReactionMode | 'all';

  // DOM
  readonly root: HTMLElement;
  private cardsContainer: HTMLElement;
  private detailsPanel: HTMLElement;
  private searchInput: DG.InputBase<string>;
  private categoryInput: DG.ChoiceInput<string | null>;

  // Events
  readonly onSelectionChanged = new Subject<NamedReaction | null>();
  readonly onReactionsChanged = new Subject<void>();

  constructor(
    private rdkit: RDModule,
    options: {
      /** Pre-filter to only show one mode */
      modeFilter?: ReactionMode | 'all';
      /** Whether to show the details panel at the bottom */
      showDetails?: boolean;
      /** Whether to show the "+ Add New" card */
      allowAdd?: boolean;
      /** Whether to show edit/delete actions */
      allowEdit?: boolean;
    } = {},
  ) {
    const {modeFilter = 'all', showDetails = true, allowAdd = true, allowEdit = true} = options;
    this.modeFilter = modeFilter;

    // ---- Search ----
    this.searchInput = ui.input.string('Search reaction:', {
      value: '',
      nullable: true,
      onValueChanged: () => this.applyFilters(),
    });
    this.searchInput.root.style.marginBottom = '0';
    this.searchInput.root.style.flexGrow = '1';
    this.searchInput.root.style.minWidth = '200px';
    ui.tooltip.bind(this.searchInput.input, 'Search reactions by name, description, tags, or category');
    (this.searchInput.input as HTMLInputElement).placeholder = 'Search reactions...';

    // ---- Category filter ----
    this.categoryInput = ui.input.choice('', {
      value: 'All',
      items: ['All'],
      nullable: false,
      onValueChanged: () => this.applyFilters(),
    }) as DG.ChoiceInput<string | null>;
    this.categoryInput.root.style.flexShrink = '0';
    ui.tooltip.bind(this.categoryInput.input, 'Filter reactions by category');

    // ---- Filter bar ----
    const filterBar = ui.divH([
      this.searchInput.root,
      this.categoryInput.root,
    ], {style: {gap: '8px', alignItems: 'center', paddingBottom: '8px', flexWrap: 'nowrap'}});

    // ---- Cards ----
    this.cardsContainer = ui.div([], {style: {
      display: 'flex', flexWrap: 'wrap', gap: '8px', overflow: 'auto',
      maxHeight: '340px', padding: '4px',
    }});

    // ---- Details ----
    this.detailsPanel = ui.div([], {style: {
      borderTop: '1px solid var(--grey-2)', paddingTop: '6px', marginTop: '6px',
      display: showDetails ? 'block' : 'none',
    }});

    // ---- Root ----
    this.root = ui.divV([filterBar, this.cardsContainer, this.detailsPanel]);
    this.root.style.minWidth = '450px';

    // Initial load
    this.reload(allowAdd, allowEdit);
  }

  /** Get the currently selected reaction. */
  get selected(): NamedReaction | null {
    return this.selectedReaction;
  }

  /** Reload reactions from storage. */
  async reload(allowAdd = true, allowEdit = true): Promise<void> {
    try {
      this.reactions = await loadAllReactions();
    } catch (e) {
      console.error('Chem | ReactionBrowser: failed to load reactions', e);
      this.reactions = [];
    }

    // Update category dropdown
    const categories = getReactionCategories(this.reactions);
    (this.categoryInput as DG.ChoiceInput<string>).items = ['All', ...categories];
    this.categoryInput.value = 'All';

    this.applyFilters(allowAdd, allowEdit);
  }

  private applyFilters(allowAdd = true, allowEdit = true): void {
    const query = (this.searchInput.value ?? '').toLowerCase().trim();
    const category = this.categoryInput.value;

    this.filteredReactions = this.reactions.filter((r) => {
      if (this.modeFilter !== 'all' && r.mode !== this.modeFilter)
        return false;
      if (category && category !== 'All' && r.category !== category)
        return false;
      if (query) {
        const haystack = [r.name, r.description ?? '', r.category, ...(r.tags ?? [])].join(' ').toLowerCase();
        if (!haystack.includes(query))
          return false;
      }
      return true;
    });

    this.renderCards(allowAdd, allowEdit);
  }

  private renderCards(allowAdd: boolean, allowEdit: boolean): void {
    this.cardsContainer.innerHTML = '';

    for (const reaction of this.filteredReactions) {
      const card = this.createReactionCard(reaction, allowEdit);
      this.cardsContainer.append(card);
    }

    if (allowAdd) {
      const addCard = this.createAddNewCard();
      this.cardsContainer.append(addCard);
    }
  }

  private createReactionCard(reaction: NamedReaction, allowEdit: boolean): HTMLElement {
    // Preview canvas
    const r = window.devicePixelRatio;
    const canvas = ui.canvas(PREVIEW_WIDTH * r, PREVIEW_HEIGHT * r);
    canvas.style.width = `${PREVIEW_WIDTH}px`;
    canvas.style.height = `${PREVIEW_HEIGHT}px`;
    canvas.style.borderRadius = '3px';
    canvas.style.backgroundColor = 'white';

    // Resolve variables with defaults for card preview
    let previewSmarts = reaction.reactionSmarts;
    if (reaction.variables) {
      const defaults: Record<string, any> = {};
      for (const [key, varDef] of Object.entries(reaction.variables))
        defaults[key] = varDef.defaultValue;
      previewSmarts = resolveReactionVariables(previewSmarts, defaults);
    }

    // Render reaction preview
    try {
      renderReactionToCanvas(this.rdkit, canvas, previewSmarts, canvas.width, canvas.height);
    } catch {
      const ctx = canvas.getContext('2d')!;
      ctx.save();
      ctx.scale(r, r);
      ctx.fillStyle = '#cc0000';
      ctx.font = '11px Roboto, sans-serif';
      ctx.fillText('Preview error', 8, PREVIEW_HEIGHT / 2);
      ctx.restore();
    }

    // Labels
    const nameLabel = ui.label(reaction.name);
    nameLabel.style.fontWeight = '500';
    nameLabel.style.fontSize = '12px';
    nameLabel.style.overflow = 'hidden';
    nameLabel.style.textOverflow = 'ellipsis';
    nameLabel.style.whiteSpace = 'nowrap';
    nameLabel.style.maxWidth = `${CARD_WIDTH - 10}px`;

    const categoryLabel = ui.label(reaction.category);
    categoryLabel.style.fontSize = '10px';
    categoryLabel.style.color = 'var(--grey-4)';

    const sourceLabel = ui.label(reaction.isUserDefined ? 'User' : 'Default');
    sourceLabel.style.fontSize = '9px';
    sourceLabel.style.color = reaction.isUserDefined ? 'var(--blue-1)' : 'var(--grey-4)';

    // Card container
    const card = ui.divV([canvas, nameLabel, categoryLabel, sourceLabel], {style: {
      width: `${CARD_WIDTH}px`,
      minHeight: `${CARD_HEIGHT}px`,
      border: '1px solid var(--grey-2)',
      borderRadius: '6px',
      padding: '5px',
      cursor: 'pointer',
      transition: 'border-color 0.15s, box-shadow 0.15s',
      overflow: 'hidden',
    }});

    // Tooltip with full details
    const tooltipContent = [
      reaction.name,
      reaction.description ?? '',
      `Category: ${reaction.category}`,
      `Mode: ${reaction.mode}`,
      reaction.tags?.length ? `Tags: ${reaction.tags.join(', ')}` : '',
    ].filter((s) => s).join('\n');
    ui.tooltip.bind(card, tooltipContent);

    // Selection
    card.addEventListener('click', () => {
      this.selectReaction(reaction, allowEdit);
      // Highlight selected card
      this.cardsContainer.querySelectorAll('.reaction-card-selected').forEach((el) => {
        el.classList.remove('reaction-card-selected');
        (el as HTMLElement).style.borderColor = 'var(--grey-2)';
        (el as HTMLElement).style.boxShadow = 'none';
      });
      card.classList.add('reaction-card-selected');
      card.style.borderColor = 'var(--blue-1)';
      card.style.boxShadow = '0 0 0 2px var(--blue-1, #1976d2)20';
    });

    // Hover effect
    card.addEventListener('mouseenter', () => {
      if (!card.classList.contains('reaction-card-selected'))
        card.style.borderColor = 'var(--grey-3)';
    });
    card.addEventListener('mouseleave', () => {
      if (!card.classList.contains('reaction-card-selected'))
        card.style.borderColor = 'var(--grey-2)';
    });

    return card;
  }

  private createAddNewCard(): HTMLElement {
    const icon = ui.iconFA('plus', () => {}, 'Add a new reaction');
    icon.style.fontSize = '32px';
    icon.style.color = 'var(--grey-3)';

    const label = ui.label('Add New');
    label.style.color = 'var(--grey-4)';
    label.style.marginTop = '8px';

    const card = ui.divV([icon, label], {style: {
      width: `${CARD_WIDTH}px`,
      minHeight: `${CARD_HEIGHT}px`,
      border: '2px dashed var(--grey-2)',
      borderRadius: '6px',
      padding: '5px',
      cursor: 'pointer',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      transition: 'border-color 0.15s',
    }});

    ui.tooltip.bind(card, 'Create a new custom reaction');

    card.addEventListener('click', () => {
      openReactionEditor(this.rdkit, {
        defaultMode: this.modeFilter !== 'all' ? this.modeFilter : undefined,
        onSave: async () => {
          invalidateReactionsCache();
          await this.reload();
          this.onReactionsChanged.next();
        },
      });
    });

    card.addEventListener('mouseenter', () => {card.style.borderColor = 'var(--blue-1)';});
    card.addEventListener('mouseleave', () => {card.style.borderColor = 'var(--grey-2)';});

    return card;
  }

  private selectReaction(reaction: NamedReaction, allowEdit: boolean): void {
    this.selectedReaction = reaction;
    this.onSelectionChanged.next(reaction);
    this.renderDetails(reaction, allowEdit);
  }

  private renderDetails(reaction: NamedReaction, allowEdit: boolean): void {
    this.detailsPanel.innerHTML = '';

    // Row 1: Name (left) + action buttons (right)
    const nameEl = ui.label(reaction.name);
    nameEl.style.fontWeight = '600';
    nameEl.style.fontSize = '13px';
    nameEl.style.whiteSpace = 'nowrap';
    nameEl.style.overflow = 'hidden';
    nameEl.style.textOverflow = 'ellipsis';

    const actions = ui.divH([], {style: {gap: '6px', flexShrink: '0'}});
    if (allowEdit) {
      if (reaction.isUserDefined) {
        actions.append(ui.button('Edit', () => {
          openReactionEditor(this.rdkit, {
            preset: reaction,
            onSave: async () => {
              invalidateReactionsCache();
              await this.reload();
              this.onReactionsChanged.next();
            },
          });
        }, 'Edit this reaction'));

        actions.append(ui.button('Delete', async () => {
          if (!confirm(`Delete reaction "${reaction.name}"?`))
            return;
          try {
            await deleteUserReaction(reaction.id);
            invalidateReactionsCache();
            this.selectedReaction = null;
            this.detailsPanel.innerHTML = '';
            await this.reload();
            this.onReactionsChanged.next();
            this.onSelectionChanged.next(null);
            grok.shell.info(`Reaction "${reaction.name}" deleted.`);
          } catch (e: any) {
            grok.shell.error(`Failed to delete: ${e?.message ?? e}`);
          }
        }, 'Delete this user-defined reaction'));
      }

      actions.append(ui.button('Duplicate', () => {
        const copy: NamedReaction = {
          ...reaction,
          id: '',
          name: `${reaction.name} (copy)`,
          isUserDefined: true,
        };
        openReactionEditor(this.rdkit, {
          preset: copy,
          onSave: async () => {
            invalidateReactionsCache();
            await this.reload();
            this.onReactionsChanged.next();
          },
        });
      }, 'Create a copy of this reaction for editing'));
    }

    const headerRow = ui.divH([nameEl, actions], {style: {
      justifyContent: 'space-between', alignItems: 'center', gap: '8px',
    }});

    // Row 2: description + mode/category info, compact
    const infoParts: string[] = [];
    if (reaction.description) infoParts.push(reaction.description);
    infoParts.push(`Mode: ${reaction.mode}`);
    infoParts.push(`Category: ${reaction.category}`);
    if (reaction.author) infoParts.push(`Author: ${reaction.author}`);
    const infoEl = ui.divText(infoParts.join('  |  '));
    infoEl.style.fontSize = '11px';
    infoEl.style.color = 'var(--grey-4)';
    infoEl.style.whiteSpace = 'nowrap';
    infoEl.style.overflow = 'hidden';
    infoEl.style.textOverflow = 'ellipsis';

    // Row 3: SMARTS (truncated, single line)
    const smartsEl = ui.element('code');
    smartsEl.textContent = reaction.reactionSmarts;
    smartsEl.style.fontSize = '10px';
    smartsEl.style.display = 'block';
    smartsEl.style.padding = '2px 4px';
    smartsEl.style.backgroundColor = 'var(--grey-1)';
    smartsEl.style.borderRadius = '3px';
    smartsEl.style.whiteSpace = 'nowrap';
    smartsEl.style.overflow = 'hidden';
    smartsEl.style.textOverflow = 'ellipsis';
    ui.tooltip.bind(smartsEl, reaction.reactionSmarts);

    this.detailsPanel.append(ui.divV([headerRow, infoEl, smartsEl], {style: {gap: '2px'}}));
  }
}
