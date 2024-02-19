// Stemming-based search specific constants

/** Ratio specifying most common words */
export const RATIO = 0.05;

/** A big number. */
export const INFTY = 100000;

/** A small number. */
export const TINY = 0.000001;

/** Frequency for polar angle */
export const POLAR_FREQ = 10;

/** Minimal number of characters in a word. */
export const MIN_CHAR_COUNT = 1;

/** Default metric weight. */
export const DEFAULT_WEIGHT = 1;

/** Number of the closest elements found by stemming-based search */
export const CLOSEST_COUNT = 6;

/** Words that are skipped, when applying stemming-based search. */
export const STOP_WORDS = ['a', 'about', 'above', 'after', 'again', 'against',
  'all', 'am', 'an', 'and', 'any', 'are', 'as', 'at',
  'be', 'because', 'been', 'before', 'being', 'below', 'between', 'both', 'but', 'by',
  'could',
  'did', 'do', 'does', 'doing', 'down', 'during',
  'each',
  'few', 'for', 'from', 'further',
  'had', 'has', 'have', 'having', 'he', 'he\'d', 'he\'ll', 'he\'s', 'her', 'here',
  'here\'s', 'hers', 'herself', 'him', 'himself', 'his', 'how', 'how\'s',
  'i', 'i\'d', 'i\'ll', 'i\'m', 'i\'ve', 'if', 'in', 'into', 'is', 'it', 'it\'s', 'its', 'itself', 'i.e',
  'let\'s',
  'me', 'more', 'most', 'my', 'myself',
  'nor', 'nan', 'no', 'not',
  'of', 'on', 'once', 'only', 'or', 'other', 'ought', 'our', 'ours', 'ourselves', 'out', 'over', 'own',
  'same', 'she', 'she\'d', 'she\'ll', 'she\'s', 'should', 'so', 'some', 'such',
  'than', 'that', 'that\'s', 'the', 'their', 'theirs', 'them', 'themselves', 'then',
  'there', 'there\'s', 'these', 'they', 'they\'d', 'they\'ll',
  'they\'re', 'they\'ve', 'this', 'those', 'through', 'to', 'too',
  'under', 'until', 'up',
  'very',
  'was', 'we', 'we\'d', 'we\'ll', 'we\'re', 'we\'ve', 'were',
  'what', 'what\'s', 'when', 'when\'s', 'where', 'where\'s', 'which',
  'while', 'who', 'who\'s', 'whom', 'why', 'why\'s', 'with', 'would',
  'you', 'you\'d', 'you\'ll', 'you\'re', 'you\'ve', 'your', 'yours', 'yourself', 'yourselves'];

/** Separator characters. */
export const SEPARATORS = '[].,:;?!(){} \n\t#';

/** Delimeter of stemming-based search results */
export const DELIMETER = '________________________________';

/** Text semantic type */
export const TEXT_SEM_TYPE = 'Text';
