import * as Vue from 'vue';
import type {EntrySourceKind, MatchConfidence} from './types';

// chip styling lives in RunComparison.css (c2-comparison-badge + these modifiers)
export const sourceBadge = (kind: EntrySourceKind) =>
  <span class={`c2-comparison-badge c2-badge-${kind}`}>{kind}</span>;

export const confidenceBadge = (confidence: MatchConfidence) =>
  <span class={`c2-comparison-badge c2-confidence-${confidence}`}>{confidence}</span>;

export const renderListFilter = (value: Vue.Ref<string>, placeholder: string) => (
  <input
    type='text'
    placeholder={placeholder}
    value={value.value}
    onInput={(e: Event) => value.value = (e.target as HTMLInputElement).value}
    style={{
      border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 6px',
      outline: 'none', width: '160px',
    }}
  />
);
