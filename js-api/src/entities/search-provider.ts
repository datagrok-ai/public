/**
 * Search provider types.
 * @module entities/search-provider
 */

import {ViewBase} from "../views/view";


export type ViewSpecificSearchProvider = {
  isApplicable?: (s: string) => boolean;
  returnType?: string;
  search: (s: string, view?: ViewBase) => Promise<{priority?: number, results: any} | null>;
  getSuggestions?: (s: string) => {priority?: number, suggestionText: string, suggestionValue?: string}[] | null;
  onValueEnter?: (s: string, view?: ViewBase) => Promise<void>;
  name: string;
  description?: string;
  options?: {relatedViewName?: string, widgetHeight?: number, [key: string]: any}
}

/** the home view should be under the name home and rest should match the views that they are applicable to */
export type SearchProvider = {
  [view: string]: ViewSpecificSearchProvider | ViewSpecificSearchProvider[];
}
