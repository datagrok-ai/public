// Revvity Signals REST API - Response Interfaces

// Top-level response interface
export interface RevvityApiResponse<T = any, I = any> {
  links?: RevvityLinks;
  data?: RevvityData<T>;
  included?: RevvityIncluded<I>[];
  errors?: RevvityApiError[];
}

// Links object (can be extended as needed)
export interface RevvityLinks {
  self: string;
  [key: string]: any;
}

// Data object
export interface RevvityData<T = any> {
  type: string;
  id: string;
  links?: RevvityLinks;
  attributes?: T;
  relationships?: RevvityRelationships;
}

// Included object (same structure as data, but can be more generic)
export type RevvityIncluded<I = any> = RevvityData<I>;

// Relationships object (can be extended as needed)
export interface RevvityRelationships {
  [key: string]: {
    links?: RevvityLinks;
    data?: { type: string; id: string } | { type: string; id: string }[];
    meta?: any;
  };
}

// Error object
export interface RevvityApiError {
  status: string;
  code: string;
  title: string;
  detail: string;
} 