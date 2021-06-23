<!-- TITLE: Text functions -->
<!-- SUBTITLE: -->

# Text functions

As parameters of the function, you can pass a regular string or a column name. To pass a regular string, it must be enclosed in double quotes, like `"my string"`. To pass a column cell, you can use the syntax `${columnName}`.

>Note that character and substring indexing in strings and lists starts at 0.

*Function List:*

- [Contains](#contains)
- [EndsWith](#endswith)
- [Length](#length)
- [ParseFloat](#parsefloat)
- [ParseInt](#parseint)
- [RegExpExtract](#regexpextract)
- [RegExpReplace](#rexpreplace)
- [ReplaceAll](#replaceall)
- [SplitString](#splitstring)
- [StartsWith](#startswith)
- [StrFind](#strfind)
- [StrLeft](#strleft)
- [StrRight](#strright)
- [StrRepeat](#strrepeat)
- [Substring](#substring)
- [ToLowerCase](#tolowercase)
- [ToUpperCase](#touppercase)
- [Trim](#trim)

## Contains(`s`, `sub`)

Checks if the string `s` contains a match of `sub`.

```javascript
/* Example of using the Contains */
Contains("Stormy Weather", "Sunny")      // Returns boolean false
Contains("Stormy Weather", "Weather")    // Returns boolean true
```

## EndsWith(`s`, `postfix`)

Checks if the string `s` ends with a match of `postfix`.

```javascript
/* Example of using the EndsWith */
EndsWith("White Christmas", "White")        // Returns boolean false
EndsWith("White Christmas", "Christmas")    // Returns boolean true
```

## Length(`s`)

Returns the length of the string `s`.

```javascript
/* Example of using the Length */
Length("Text")    // Returns integer 4
```

## ParseFloat(`s`)

Parse `s` as a, possibly signed, real literal and return its value.

```javascript
/* Example of using the ParseFloat */
ParseFloat("2025")      // Returns real -2025
ParseFloat("-012.78")   // Returns real -12.78
```

## ParseInt(`s`)

Parse `s` as a, possibly signed, integer literal and return its value.

```javascript
/* Example of using the ParseInt */
ParseInt("2025")    // Returns integer 2025
ParseInt("-012")     // Returns integer -12
```

## RegExpExtract(`s`, `pattern`, `n`)

Returns the `n`-th part of a string `s` that matches a regular expression `pattern`.

```javascript
/* Example of using the RegExpExtract */
RegExpExtract("Hello World!", "l+", 0)    // Returns string "ll"
RegExpExtract("Hello World!", "l+", 1)    // Returns string "l"
```

## RegExpReplace(`s`, `pattern`, `sub`)

Returns the string after replacing a substring `sub` in string `s` according to a regular expression `pattern`.

```javascript
/* Example of using the RegExpReplace */
RegExpReplace("Hello World!", "l+", "LL")    // Returns string "HeLLo WorLLd!"
```

## ReplaceAll(`s`, `from`, `replace`)

Replaces all substrings in string `s` that match `from` with `replace` and returns the result obtained.

```javascript
/* Example of using the ReplaceAll */
ReplaceAll("New York, New York", "York", "Orleans")    // Returns string "New Orleans, New Orleans"
ReplaceAll("Every", "", ".")                           // Returns string ".E.v.e.r.y."
ReplaceAll("moto", "o", "")                            // Returns string "mt"
```

## SplitString(`s`, `separator`, `i`)

Splits the string `s` at matches of `separator` and returns the `i`-th of the substring between the matches.

```javascript
/* Example of using the SplitString */
SplitString("Born to Be Wild", " ", 2)     // Returns string "Be"
SplitString("Born to Be Wild", "to", 1)    // Returns string " Be Wild"
SplitString("a,b,c,d", ",", 0)             // Returns string "a"
```

## StartsWith(`s`, `prefix`)

Checks if the string `s` starts with a match of `prefix`.

```javascript
/* Example of using the StartsWith */
StartsWith("White Christmas", "White")        // Returns boolean true
StartsWith("White Christmas", "Christmas")    // Returns boolean false
```

## StrFind(`s`, `sub`)

Returns the index of the first occurrence of the string `sub` in the string `s`.

If `sub` or `s` are empty then the function returns -1.

```javascript
/* Example of using the StrFind */
StrFind("Hello world!", "Hello")  // Returns integer 0
StrFind("Hello world!", "world")  // Returns integer 6
StrFind("Hello world!", "Car")    // Returns integer -1
StrFind("", "Moon")               // Returns integer -1
```

## StrLeft(`s`, `count`)

Returns the first `count` characters of the string `s`.

If `s` is empty then the function returns empty string.

`count` is of type integer. If `count` is negative, then the `count` number of characters will be removed from the right-hand side of the `s` string.

```javascript
/* Example of using the StrLeft */
StrLeft("Daddy", 1)      // Returns string "D"
StrLeft("Daddy", 10))    // Returns string "Daddy"
StrLeft("", 3)           // Returns empty string
StrLeft("Daddy", -3)     // Returns string "Da"
StrLeft("Daddy", -8))    // Returns empty string
```

## StrRight(`s`, `count`)

Returns the last `count` characters of the string `s`.

If `s` is empty then the function returns empty string.

`count` is of type integer. If `count` is negative, then the `count` number of characters will be removed from the left-hand side of the `s` string.

```javascript
/* Example of using the StrRight */
StrRight("Daddy", 1)      // Returns string "y"
StrRight("Daddy", 10))    // Returns string "Daddy"
StrRight("", 3)           // Returns empty string
StrRight("Daddy", -3)     // Returns string "dy"
StrRight("Daddy", -8))    // Returns empty string
```

## StrRepeat(`s`, `num`)

Returns a string consisting of `num` repetitions of string `s`.

```javascript
/* Example of using the StrRepeat */
StrRepeat("Chain", 2)    // Returns string "ChainChain"
```

## Substring(`s`, `start`, `end`)

Returns the substring of the string `s` from `start` index, inclusive, to `end`, exclusive.

```javascript
/* Example of using the Substring */
Substring("Alfa Romeo", 5, 10)    // Returns string "Romeo"
```

## ToLowerCase(`s`)

Converts all characters in this string `s` to lower case and returns the result obtained.

```javascript
/* Example of using the ToLowerCase */
ToLowerCase("ALPHABET")    // Returns string "alphabet"
```

## ToUpperCase(`s`)

Converts all characters in this string `s` to upper case and returns the result obtained.

```javascript
/* Example of using the ToUpperCase */
ToUpperCase("alphabet")    // Returns string "ALPHABET"
```

## Trim(`s`)

Returns the string without any leading and trailing whitespace of the string `s`.

```javascript
/* Example of using the Trim */
Trim("  My home.   ")    // Returns string "My home."
```
