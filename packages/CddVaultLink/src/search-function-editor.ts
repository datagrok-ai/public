import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ANY_RUN_CHOICE, CDD_SEARCH_TYPES, CDDVaultSearchType, PROTOCOL_RUN_COND_CHOICES, ProtocolRunCond} from './constants';
import {Protocol} from './cdd-vault-api';
import {funcs} from './package-api';

const LAST_SEARCH_KEY = 'CDDVaultLink.lastSearch';

type StoredSearch = {
  protocolName?: string;
  runChoice?: string;
  runName?: string;
  structure?: string;
  searchType?: string;
  similarityThreshold?: number;
};

type CDDSearchParams = {
  structure?: string;
  structure_search_type?: CDDVaultSearchType;
  structure_similarity_threshold?: number;
  protocol?: number;
  run?: number;
}

export class SearchEditor {
  vaultId;
  protocols: {[key: string]: Protocol} = {};
  currentProtocolRuns: {[key: string]: number} = {};
  protocolNames: string[] = [];

  specificProtocolDiv = ui.divV([]);
  protocolListChoice: DG.InputBase = ui.input.choice('', {value: '', items: ['']});
  protocolListChoiceDiv = ui.div('', 'cdd-protocol-list-input');

  runAndReadoutDiv = ui.divV([], {style: {display: 'none'}});

  runChoice: DG.InputBase = ui.input.choice('', {
    value: PROTOCOL_RUN_COND_CHOICES[0], items: PROTOCOL_RUN_COND_CHOICES, onValueChanged: () => {
      this.runListChoiceDiv.style.display = this.runChoice.value === ProtocolRunCond.SPECIFIC ? 'flex' : 'none';
    },
  });
  runListChoiceDiv = ui.div();
  runsListChoice: DG.InputBase = ui.input.choice('', {value: '', items: ['']});

  structureInput = ui.input.molecule('Structure');
  structureSearchTypeChoice = ui.input.choice('Search', {
    value: CDD_SEARCH_TYPES[0], items: CDD_SEARCH_TYPES, onValueChanged: () => {
      this.similarityThreshold.root.style.display = this.structureSearchTypeChoice.value === CDDVaultSearchType.SIMILARITY ? 'flex' : 'none';
    },
  });
  similarityThreshold = ui.input.float('', {value: 0.7, min: 0, max: 1, step: 0.05, showSlider: true});

  accordion: DG.Accordion;

  initComplete: Promise<void>;
  hasRestoredSearch = false;

  constructor(vaultId: number) {
    this.vaultId = vaultId;
    this.accordion = ui.accordion(`cdd_vault_search_${this.vaultId}`);
    this.accordion.root.classList.add('cdd-search-vault-tab');
    this.structureInput.classList.add('cdd-molecule-input');
    this.runListChoiceDiv.append(this.runsListChoice.root);
    this.runAndReadoutDiv.append(ui.divH([
      this.runChoice.root,
      this.runListChoiceDiv,
    ], {style: {gap: '20px'}}));

    this.protocolListChoiceDiv.append(this.protocolListChoice.root);
    this.specificProtocolDiv.append(this.protocolListChoiceDiv);
    this.specificProtocolDiv.append(this.runAndReadoutDiv);

    this.accordion.addPane('Protocol', () => ui.divH([
      this.specificProtocolDiv,
    ]));

    this.similarityThreshold.root.style.display = 'none';
    this.accordion.addPane('Structure', () => ui.divV([
      ui.divH([this.structureSearchTypeChoice.root, this.similarityThreshold.root], {style: {gap: '20px'}}),
      this.structureInput,
    ]));

    this.initComplete = this.init();
  }

  async init() {
    try {
      const protocolsStr = await funcs.getProtocolsAsync(this.vaultId, 1);
      const protocols = protocolsStr !== '' ? JSON.parse(protocolsStr) as Protocol[] : [];
      if (protocols && protocols.length) {
        protocols.forEach((it) => this.protocols[it.name] = it);
        this.protocolNames = Object.keys(this.protocols);
        ui.empty(this.protocolListChoiceDiv);
        this.protocolListChoice = ui.input.choice('', {
          value: '', items: [''].concat(this.protocolNames), onValueChanged: () => {
            this.runAndReadoutDiv.style.display = !this.protocolListChoice!.value ? 'none' : 'flex';
            if (this.protocolListChoice!.value) {
              this.runChoice.value = ANY_RUN_CHOICE;
              this.updateRuns(this.protocolListChoice!.value);
            }
          },
        });
        this.protocolListChoiceDiv.append(this.protocolListChoice.root);
      }
      this.restoreLastSearch();
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
    }
  }

  private restoreLastSearch(): void {
    const raw = grok.userSettings.getValue(LAST_SEARCH_KEY, String(this.vaultId));
    if (!raw) return;
    let s: StoredSearch;
    try {
      s = JSON.parse(raw);
    } catch { return; }
    this.hasRestoredSearch = !!(s.protocolName || s.structure);

    if (s.protocolName && this.protocolNames.includes(s.protocolName)) {
      this.protocolListChoice.value = s.protocolName;
      this.runAndReadoutDiv.style.display = 'flex';
      this.updateRuns(s.protocolName);
      if (s.runChoice && (PROTOCOL_RUN_COND_CHOICES as readonly string[]).includes(s.runChoice)) {
        this.runChoice.value = s.runChoice;
        this.runListChoiceDiv.style.display = s.runChoice === ProtocolRunCond.SPECIFIC ? 'flex' : 'none';
        if (s.runName && s.runName in this.currentProtocolRuns)
          this.runsListChoice.value = s.runName;
      }
    }

    if (s.searchType && CDD_SEARCH_TYPES.includes(s.searchType as CDDVaultSearchType)) {
      this.structureSearchTypeChoice.value = s.searchType as CDDVaultSearchType;
      this.similarityThreshold.root.style.display = s.searchType === CDDVaultSearchType.SIMILARITY ? 'flex' : 'none';
    }
    if (typeof s.similarityThreshold === 'number')
      this.similarityThreshold.value = s.similarityThreshold;
    if (s.structure)
      this.structureInput.value = s.structure;
  }

  public saveLastSearch(): void {
    const payload: StoredSearch = {
      protocolName: this.protocolListChoice.value || undefined,
      runChoice: this.runChoice.value ?? undefined,
      runName: this.runChoice.value === ProtocolRunCond.SPECIFIC ? this.runsListChoice.value : undefined,
      structure: this.structureInput.value || undefined,
      searchType: this.structureSearchTypeChoice.value ?? undefined,
      similarityThreshold: this.similarityThreshold.value ?? undefined,
    };
    grok.userSettings.add(LAST_SEARCH_KEY, String(this.vaultId), JSON.stringify(payload));
  }

  updateRuns(protocolName: string) {
    this.currentProtocolRuns = {};
    this.protocols[protocolName].runs.forEach((it) => this.currentProtocolRuns[`${it.run_date} (${it.person})`] = it.id);
    const runsNames = Object.keys(this.currentProtocolRuns);
    ui.empty(this.runListChoiceDiv);
    this.runsListChoice = ui.input.choice('', {value: runsNames.length ? runsNames[0] : '', items: runsNames});
    this.runListChoiceDiv.append(this.runsListChoice.root);
    this.runListChoiceDiv.style.display = 'none';
  }

  public async reset(): Promise<void> {
    await this.initComplete;
    this.protocolListChoice.value = '';
    this.runAndReadoutDiv.style.display = 'none';
    this.runChoice.value = PROTOCOL_RUN_COND_CHOICES[0];
    this.runListChoiceDiv.style.display = 'none';
    this.structureInput.value = '';
    this.structureSearchTypeChoice.value = CDD_SEARCH_TYPES[0];
    this.similarityThreshold.value = 0.7;
    this.similarityThreshold.root.style.display = 'none';
    this.hasRestoredSearch = false;
    grok.userSettings.delete(LAST_SEARCH_KEY, String(this.vaultId));
  }

  public getEditor(): HTMLElement {
    return this.accordion.root;
  }

  public getParams(): CDDSearchParams {
    const isSimilarity = this.structureSearchTypeChoice.value === CDDVaultSearchType.SIMILARITY;
    return {
      structure: this.structureInput.value,
      structure_search_type: this.structureSearchTypeChoice.value ?? undefined,
      structure_similarity_threshold: isSimilarity ? (this.similarityThreshold.value ?? undefined) : undefined,
      protocol: this.protocolListChoice.value ? this.protocols[this.protocolListChoice.value].id : undefined,
      run: this.runChoice.value === ProtocolRunCond.SPECIFIC ? this.currentProtocolRuns[this.runsListChoice.value] : undefined,
    };
  }
}
