
/** If the condition is false, throws message(). */
export function assure(condition: boolean, errorMessage: string | (() => string)) {
  if (!condition)
    throw typeof errorMessage === 'string' ? errorMessage : errorMessage();
}


/**
 * Converts a number to Excel column letters (A-ZZ)
 * Examples:
 * 0 -> "A"
 * 25 -> "Z"
 * 26 -> "AA"
 * 51 -> "AZ"
 * 701 -> "ZZ"
 *
 * @param columnNumber - The zero-based column number to convert
 * @returns The Excel column letter(s)
 * @throws Error if the column number is negative or would require more than 2 letters
 */
export function getExcelColumnLetter(columnNumber: number): string {
  // Check for invalid input
  if (columnNumber < 0)
    throw new Error("Column number cannot be negative");

  if (columnNumber > 701) // ZZ = 26 * 26 + 25 = 701
    throw new Error("Column number too large - would require more than 2 letters");

  const firstCharCode = 'A'.charCodeAt(0);

  // For single letter (A-Z)
  if (columnNumber <= 25) {
    return String.fromCharCode(firstCharCode + columnNumber);
  }

  // For double letters (AA-ZZ)
  const firstLetter = Math.floor(columnNumber / 26);
  const secondLetter = columnNumber % 26;

  return String.fromCharCode(firstCharCode + firstLetter - 1) +
    String.fromCharCode(firstCharCode + secondLetter);
}
