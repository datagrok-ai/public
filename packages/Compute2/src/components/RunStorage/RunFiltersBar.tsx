import * as Vue from 'vue';
import {FilterDropdown, FilterOption} from '../Inspector/FilterDropdown';
import {RunFilters} from './use-run-filters';

/** Search, tag and author filters; every change is a new value in update:modelValue. */
export const RunFiltersBar = Vue.defineComponent({
  name: 'RunFiltersBar',
  props: {
    modelValue: {type: Object as Vue.PropType<RunFilters>, required: true},
    tagOptions: {type: Array as Vue.PropType<FilterOption[]>, required: true},
    authorOptions: {type: Array as Vue.PropType<FilterOption[]>, required: true},
    searchPlaceholder: {type: String, default: 'Search'},
  },
  emits: {
    'update:modelValue': (_value: RunFilters) => true,
  },
  setup(props, {emit}) {
    const update = (patch: Partial<RunFilters>) => emit('update:modelValue', {...props.modelValue, ...patch});
    const authorName = (id: string) => props.authorOptions.find((option) => option.value === id)?.label ?? id;

    return () => {
      const {search, tags, authors} = props.modelValue;
      return <div class='c2-runs-filters'>
        <input class='ui-input-editor c2-runs-search' placeholder={props.searchPlaceholder} value={search}
          onInput={(e: Event) => update({search: (e.target as HTMLInputElement).value})}/>
        { props.tagOptions.length > 0 &&
          <div class='c2-runs-tags'>
            <FilterDropdown options={props.tagOptions} modelValue={tags} placeholder='Tags'
              onUpdate:modelValue={(value: string[]) => update({tags: value})}/>
            { tags.map((item) =>
              <span key={item} class='c2-runs-tag' title='Click to stop filtering by this tag'
                onClick={() => update({tags: tags.filter((other) => other !== item)})}>{item} &times;</span>) }
          </div> }
        { (props.authorOptions.length > 1 || authors.length > 0) &&
          <div class='c2-runs-tags'>
            <FilterDropdown options={props.authorOptions} modelValue={authors} placeholder='Authors'
              onUpdate:modelValue={(value: string[]) => update({authors: value})}/>
            { authors.map((id) =>
              <span key={id} class='c2-runs-tag' title='Click to stop filtering by this author'
                onClick={() => update({authors: authors.filter((other) => other !== id)})}>
                {authorName(id)} &times;</span>) }
          </div> }
      </div>;
    };
  },
});
