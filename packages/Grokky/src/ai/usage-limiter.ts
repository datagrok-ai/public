import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {_package} from '../package';

const USAGE_EVENT_TYPE = 'ai-request';
const SETTINGS_NAME = 'Grokky.UsageLimiter';
const SETTINGS_KEY_DATE = 'date';
const SETTINGS_KEY_COUNT = 'count';

interface LimitsConfig {
  defaultDailyLimit: number;
  groupLimits: Record<string, number>;
}

const DEFAULT_CONFIG: LimitsConfig = {
  defaultDailyLimit: 20,
  groupLimits: {
    'AI Folks': 1000,
  },
};

export class UsageLimiter {
  private static _instance: UsageLimiter | null = null;
  private config: LimitsConfig = DEFAULT_CONFIG;
  private todayCount = 0;
  private lastCountDate = '';
  private dailyLimit: number | null = null;

  private constructor() {}

  static getInstance(): UsageLimiter {
    UsageLimiter._instance ??= new UsageLimiter();
    return UsageLimiter._instance;
  }

  async init(): Promise<void> {
    this.loadCount();
    await this.resolveDailyLimit();
  }

  private today(): string { return dayjs().format('YYYY-MM-DD'); }

  private loadCount(): void {
    const saved = grok.userSettings.get(SETTINGS_NAME, true);
    const today = this.today();
    if (saved?.[SETTINGS_KEY_DATE] === today) {
      this.todayCount = Number(saved[SETTINGS_KEY_COUNT]) || 0;
      this.lastCountDate = today;
    } else {
      this.todayCount = 0;
      this.lastCountDate = today;
      this.saveCount();
    }
  }

  private saveCount(): void {
    grok.userSettings.put(SETTINGS_NAME, {
      [SETTINGS_KEY_DATE]: this.lastCountDate,
      [SETTINGS_KEY_COUNT]: String(this.todayCount),
    }, true);
  }

  private async resolveDailyLimit(): Promise<number> {
    if (this.dailyLimit != null)
      return this.dailyLimit;
    const userGroup = await grok.dapi.groups.find(grok.shell.user.group.id);
    const groups = [...userGroup.memberships, ...userGroup.adminMemberships].map((g: DG.Group) => g.friendlyName);
    let max = this.config.defaultDailyLimit;
    for (const [group, limit] of Object.entries(this.config.groupLimits)) {
      if (groups.includes(group) && limit > max)
        max = limit;
    }
    this.dailyLimit = max;
    return max;
  }

  async checkAndIncrement(feature: string, prompt: string, model?: string): Promise<void> {
    this.loadCount();
    const limit = await this.resolveDailyLimit();
    if (this.todayCount >= limit)
      throw new Error(`You've reached your daily AI request limit (${limit}). Contact your administrator for a higher quota.`);
    this.todayCount++;
    this.saveCount();
    _package.logger.usage(`${USAGE_EVENT_TYPE}:${feature} | ${model ?? 'unknown'} | ${prompt}`);
  }

  /** Checks the limit and increments. Returns true if allowed, shows warning and returns false if over limit. */
  async tryCheckAndIncrement(feature: string, prompt: string, model?: string): Promise<boolean> {
    try {
      await this.checkAndIncrement(feature, prompt, model);
      return true;
    }
    catch (e: any) {
      grok.shell.warning(e.message);
      return false;
    }
  }

  get remaining(): number { return Math.max(0, (this.dailyLimit ?? this.config.defaultDailyLimit) - this.todayCount); }
  get limit(): number { return this.dailyLimit ?? this.config.defaultDailyLimit; }
  get used(): number { return this.todayCount; }
}
