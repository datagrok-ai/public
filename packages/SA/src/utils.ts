// utils.ts

// Error messeges
enum ERROR_MSG {
  INCORRECT_SIZE_TYPE = 'Size must be an integer.',
  SIZE_MUST_BE_POSITIVE = 'Size must be positive.',
};

export function checkSize(size: any): void {
  // check if it's a number
  if (!Number.isInteger(size))
    throw new Error(ERROR_MSG.INCORRECT_SIZE_TYPE);
  
  // check positivity
  if (size < 1)
    throw new Error(ERROR_MSG.SIZE_MUST_BE_POSITIVE);
}