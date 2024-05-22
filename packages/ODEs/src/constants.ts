// Control constants
export const CONTROL_TAG = '#';
export const CONTROL_TAG_LEN = CONTROL_TAG.length;
export const DF_NAME = 'df';
const META = `${CONTROL_TAG}meta`;
export const MAX_LINE_CHART = 4;

/** Control expressions for the problem specifying */
export enum CONTROL_EXPR {
    NAME = `${CONTROL_TAG}name`,
    TAGS = `${CONTROL_TAG}tags`,
    DESCR = `${CONTROL_TAG}description`,
    DIF_EQ = `${CONTROL_TAG}equations`,
    EXPR = `${CONTROL_TAG}expressions`,
    ARG = `${CONTROL_TAG}argument`,
    INITS = `${CONTROL_TAG}inits`,
    CONSTS = `${CONTROL_TAG}constants`,
    PARAMS = `${CONTROL_TAG}parameters`,
    TOL = `${CONTROL_TAG}tolerance`,
    LOOP = `${CONTROL_TAG}loop`,
    UPDATE = `${CONTROL_TAG}update`,
    RUN_ON_OPEN = `${META}.runOnOpen`,
    RUN_ON_INPUT = `${META}.runOnInput`,
    OUTPUT = `${CONTROL_TAG}output`,
    COMMENT = `${CONTROL_TAG}comment`,
    EPS_SCALE = `${CONTROL_TAG}scale`,
};

/** Loop consts */
export enum LOOP {
  MIN_LINES_COUNT = 1,
  COUNT_IDX = 0,
  COUNT_NAME = '_count',
  MIN_COUNT = 1,
};

/** UPDATE consts */
export enum UPDATE {
  MIN_LINES_COUNT = 1,
  DURATION_IDX = 0,
  DURATION = '_duration',
};
