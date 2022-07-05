# Basics of technical writing

## Language and grammar

Datagrok documentation should be clear and easy to understand:

* Write in US English with US grammar.
* [Be clear and concise](###Conciseness).
* [Use active voice](###active-voice): make it clear who's performing the action.
* [Put conditional clauses before instructions, not after](###clause-order).
* [Use the second person point of view](###pronouns).
* [Don't use words that minimize the difficulty of a task](###Words-that-minimize-difficulty).
* [Use sentence case for headers and bulleted lists](###capitalization).
* [Write compound words and words with prefixes without hyphens or spaces](###Compound-words-and-words-with-prefixes).
* [Contractions](###contractions) are OK because we write in an informal tone.
* [Spell out: acronyms on first use and numbers from zero to nine](###acronyms-and-numbers).
* For usage and spelling of specific words, see the [Word list](link).

### Conciseness

Conciseness is about using the fewest words while retaining completeness in meaning. You can achieve conciseness in several ways:

* Shorten words and phrases.

  | Wordy                    |Short         |
  |:-------------------------|:--------------|
  | presently                |now            |
  | provide an opportunity   |allow/permit   |
  | with regard to           |about          |

* Substitute large and complex words with short, simple words.
  
  | Complex           | Simple        |
  |:------------------|:--------------|
  | commence          | begin/start   |
  | demonstrate       | show          |
  | employ / utilize  | use           |
  
* Avoid the words ending in _-sion_ or _-tion_ (nominalizations). Replace them with active words.

  | Nominalization          | Active word       |
  |:------------------------|:------------------|
  | come to the conclusion  | conclude / decide |
  | implementation of       | implement         |
  | take into consideration | consider          |

* Omit cliches, idioms, old-fashioned, or overused words.

  | Instead of       | Use             |
  |:-----------------|:----------------|
  | at all times     | always          |
  | best-in-class    | superior        |

* Avoid word redundancies.

  |Instead of        |Use                    |
  |:-----------------|:----------------------|
  |like for example  | like _or_ for example |
  |the reason why    | why                   |
  |past experience   | experience            |

* Avoid sentences starting with _there is_, _here is_, _it is_, and similar (that is, sentences with no true subjects). Use this tips:

    1. First, identify an actor.
    1. Then, find the key action and make it a verb.
    1. Lastly, make the actor the subject of this verb.

  | Instead of                              |Use                     |
  |:----------------------------------------|:-----------------------|
  | There are many features in the platform. |The platform has many features. |
  | There are multiple options.              |You have several options.|

* Drop unnecessary words; use shorter alternatives if available.

  | Instead of                      |Use                     |
  |:--------------------------------|:-----------------------|
  | It is necessary for you to do A. | You must do A.          |
  | Follow all of the steps below.   | Follow the steps below. |

  _Tip:_ If a step is optional, start with the word _optional_ followed by a period. For example: _Step 3._ Optional. Change the size by dragging the bottom right corner.

* Simplify long, complicated sentences.

  * Look for and remove an introductory phrase that doesn't add value.

    |Instead of| Use|
    |:---|:---|
    | As discussed previously, you can embed a viewer as an iframe. | You can embed a viewer as an iframe.|
    | It is important to note that this feature is in beta. | Important: This feature is in beta.|

  * If a sentence contains multiple prepositions, rewrite the sentence to remove as many prepositions as you can.

    |Instead of| Use|
    |:---|:---|
    |These representations can be _of_ great importance _for_ the description _of_ monomers _of_ a decomposed macromolecule. | These representations are essential to describing monomers of a decomposed macromolecule.|

### Active voice

Avoid using passive voice. When you use passive voice, you imply that the doer of the action is unknown, unimportant, or cannot be named for any reason. Because technical writing is action-oriented, the use of active voice is preferred. Never end a sentence with a passive word. Place the subject and verd close to each other and put them at the front of a sentence.

| Instead of                                |Use                        |
|:------------------------------------------|:--------------------------|
|In this article, challenges in interactive data exploration are discussed. | This article addresses challenges in interactive data exploration.|  
|A dialog is displayed on clicking the icon. | When you click the icon, a dialog opens.|
| Viewers can be customized.                 | You can customize viewers. |

### Pronouns

Use the second person point of view (you) instead of a third person point of view (he/she). If you must use the third person, use plural pronouns (they, their) instead of single pronouns (he/she, his/her).

_Note_: Avoid using words that explicitly favor one gender, such as _guys_ or _businessman/businesswoman_.

### Words that minimize the difficulty

Avoid using words that minimize the difficulty involved in a task or operation, such as _just_ or _simply_. People may become frustrated when they do not find a step as straightforward or simple as it is implied to be.

### Capitalization

Use sentence case. Sentence case refers to a capitalization style in which only the first word and proper nouns or acronyms are uppercased:

* Capitalize names of third party organizations and products (such as _Slack_), methods, or methodologies (such as _Agile_).
* For headings, capitalize the first word only.
* When referring to specific user interface text, like a menu item, use the same capitalization that the user interface shows.
* For bulleted lists, always capitalize the first word, including the text that completes an explanatory or introductory text (ending with a colon). For example:

  You have several options:
  * First option
  * Second option

### Compound words and words with prefixes

We tend to use the closed form of compound words and words with prefixes; that is, write these words without a space or a hyphen (for example, _a_ [_dataset_](link) or _an_ [_open source_](link) _plugin_). For the list of commonly used terms, see [Word list](link).

### Contractions

You can use most types of contractions. Negation contractions (such as _isn't_ or _don't_) are even helpful because it's harder to misread _don't_ compared to _do not_. However, avoid contractions formed from nouns and verbs:

|Instead of                           |Use                             |
|:------------------------------------|:-------------------------------|
|The platform's fast and easy to use.  |The platform is fast and easy to use.          |

### Acronyms and numbers

If you use an acronym, spell it out on the first use on a page. You don't need to spell it out more than once on a page. If you can, avoid using acronyms in headings.

When using numbers in the text, spell out zero through nine, and use numbers for 10 and greater. For more guidance, see [Google Developer Documentation Style Guide](https://developers.google.com/style/numbers?hl=en).

### Clause order

When giving instructions, follow this clause order:

1. goal first
1. then location
1. then action

Mentioning the goal or circumstance first lets the user skip the instruction if it doesn't apply.

|Instead of                                                  |Use       |
|:-----------------------------------------------------------|:---------|
| Click **Remove** on the **Sidebar** to delete the document. | To delete the document, on the **Sidebar**, click **Remove**.|

### Paragraphs

The first sentence of a paragraph is a topic sentence describing that paragraph's central idea. State that idea upfront, then build on it. When reading technical documentation, users often look for a specific piece of information. When they know the topic, users can skip the paragraph if that paragraph's information is irrelevant to their needs.

Never assume users have read every word in the preceding paragraph. Don't start a paragraph with a pronoun or words like _but_ or _however_. Every paragraph should encapsulate an idea that can stand on its own.

Follow the 1:1 rule: one paragraph, one idea. Remove any sentences that neither clarify the idea in that paragraph, nor logically connect this paragraph and the next.

End the paragraph with a sentence that summarizes the topic sentence.

## Punctuation

Follow these guidelines for punctuation:

* Don't use semicolons or dashes. Use two sentences instead.
* In a list of three or more items, use Oxford commas before the final _and_ or _or_.
* Don't use curly quotes. Use straight quotes instead.
* For bulleted lists:
  * Don't add commas or semicolons to the ends of list items.
  * If a list item is a complete sentence (with a subject and a verb), add a period in the end.
  * Use no punctuation after bullets that are not complete sentences.
  * Be consistent. Don't mix sentences and individual items in a list.
  * Separate list items from introductory or explanatory text with a colon. Also, for numbered lists, italicize the numbered list items to make them visually stand out.

    For example:

    This procedure has two steps:
    * _Step 1._ Do the first thing.
    * _Step 2._ Do the second thing.

### Commas and appositives

An appositive is a noun or a phrase placed next to another noun or a phrase to modify it. Appositives can be restrictive and nonrestrictive. The following table summarises the use of commas with appositives:

|                |Restrictive             |Nonrestrictive                |
|:---------------|:-----------------------|:-----------------------------|
|Purpose         |Narrow down the meaning |Provide additional information|
|When removed    | Meaning changes        | Meaning doesn't change       |
|Comma usage     | Don't wrap in commas   | Wrap in commas               |
|Example         | Only users _who have preinstalled the package_ can access this feature.    | Data augmentation, _one of the platform capabilities_, is used to push insights to the users.|
