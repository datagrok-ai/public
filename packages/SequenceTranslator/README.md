# SequenceTranslator

SequenceTranslator is a [package](https://datagrok.ai/help/develop/develop#packages) for
the [Datagrok](https://datagrok.ai) platform, used to translate [oligonucleotide](https://en.wikipedia.org/wiki/Oligonucleotide) 
sequences between [different representations]().

Oligonucleotides are characterized by the sequence of [nucleotide modifications]() that make up the entire molecule. 
The length of the oligonucleotide is usually 13-25 nucleotides long (denoted by "-mer").

# Use Cases
| Name |  Steps  |
|---|---|
| Translate one sequence of types list | Paste sequence into the text field in `MAIN` tab of App |
| Translate many sequences of types list | 1. Open [SequenceTranslator](https://public.datagrok.ai/apps/Sequencetranslator/) <br/>2. Drag & drop an Excel or CSV file with sequences into Datagrok. The platform will detect columns with sequences if they match corresponding RegExp rules <br/>3. Right-click on the column header, then see the `Convert` menu <br/>4. Choose translation function <br/>This will add the result column to the right of the table|
| Translate sequences using new Axolabs pattern | 1. Go to `AXOLABS` tab <br/>2. Drag & drop an Excel or CSV file with sequences into Datagrok <br/>3. Choose your table in the `Table` field <br/>4. Choose `SS Column`, and, if needed, `AS Column` and `ID Column`(needed to add column with concatenated ID and pattern name) <br/>5. Press `Convert Sequences` button  <br/>This will add the result column(s) to the right of the table. Save pattern, if needed|
| Translate sequences using existing Axolabs pattern | 1. Go to `AXOLABS` tab <br/>2. Select your pattern in `Load Pattern` field <br/>3. Do steps #2-5 from previous use case |

# Sequence Representations 
| Name | Example | Regular expression |
|---|---|---|
| DNA nucleotides | `AGGTCTTCATGACTTCGGCC` | `^[ATGC]{10,}$` |
| RNA nucleotides | `UUCAACUGCUUACGUCUUU` | `^[AUGC]{10,}$` |
| BioSpring / Gapmers | `6*8*8*5*7*T*T*9*A*T*G*A*9*T*T*7*8*8*7*7` | `^[*56789ATGC]{30,}$` |
| GCRS / Gapmers | `moeAnpsmoeGnpsmoeGnpsmoeUnpsmoe5mCnpsTpsTpsCpsApsTpsGpsApsCpsTpsTpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC` | `^(?=.*moe)(?=.*5mC)(?=.*ps){30,}` |
| BioSpring / siRNA | `5*1*766354715274575*5*5` | `^[*1-8]{30,}$` |
| Axolabs / siRNA | `usUfscaaCfuGfcUfuAfcGfucususu` | `^[fsACGUacgu]{20,}$` |
| GCRS | `mUpsfUpsmCmAmAfCmUfGmCfUmUfAmCfGmUmCmUpsmUpsmU` | `^[fmpsACGU]{30,}$` |
| OP100 | `ug*aa*uu*ag*ag*ga*ga*cg*ga*cac` | `^[acgu*]{10,}$` |
| MM12 | `kgKHKGLHIELGJFHKiehK` | `^[IiJjKkLlEeFfGgHhQq]{10,}$` |
| ABI | `58877TTTACCACGT56788` | `^[5678ATGC]{10,}$` |

# Axolabs Nucleotide Modifications
| Name | Description | Symbols |
|---|---|---|
| `RNA` | RNA nucleotides | `A, C, G, U` |
| `DNA` | DNA nucleotides | `dA, dC, dG, dT` |
| `2'-Fluoro` | 2'-Fluoro nucleotides | `Af, Cf, Gf, Uf` |
| `2'-O-Methyl` | 2'-O-Methyl nucleotides | `a, c, g, u` |
| `2'-O-MOE` | 2'-O-MOE nucleotides (including 5-Methyl C) | `Am, Cm, Gm, Tm` |
| `GNA` | Glycol nucleic acid | `(GNA-A), (GNA-C), (GNA-G), (GNA-T)` |
| `LNA` | Locked nucleic acid (including 5-Methyl C) | `Ab, Cb, Gb, Tb` |
| `UNA` | Unlocked nucleotides | `Ao, Co, Go, Uo` |
| `A` | Adenine | `a` |
| `C` | Cytosine | `c` |
| `G` | Guanine | `g` |
| `U` | Uracil | `u` |
| `X-New` | | `X` |
| `Y-New` | | `Y` |
| `Z-New` | | `Z` |
| `InvAbasic` | Inverted abasic capped | `(invabasic)` |
| `5'-vinylps` | 5'-vinylphosphonate-2'-OMe-uridine | `(vinu)` |
| `InvAbasic(o)` | Inverted abasic capped (overhang) | `(invabasic)` |
| `2'-OMe-U(o)` | Nucleotide Uridine with 2’O-Methyl protection (overhang) | `mU` |

| Term |  Definition  |
|---|---|
| Pattern | Defines rules of translation. Contains: pattern name |
| PTO | Indicates whether oligonucleotide is phosphorothioated (ps linkage) |
| Sense strand (SS) | Contains the exact nucleotide sequence to the mRNA which encodes for a functional protein. Has the information that would be readable on the RNA, and that's called the coding side.
|
| Antisense strand (AS) | Non-coding DNA strand of a gene. A cell uses antisense DNA strand as a template for producing messenger RNA (mRNA) that directs the synthesis of a protein. |

Oligonucleotides are chemically synthesized. Chain assembly proceeds in the 3' to 5' direction by following a routine 
procedure referred to as a "synthetic cycle". 
Completion of a single synthetic cycle results in the addition of one nucleotide residue to the growing chain.

Uploaded the latest Sequences and Conversion  file , this includes for  information on 4 synthesizers that are listed 
on ABI, OP100,MerMade tab .

In th modification section on the right side of the screen you can select modification for each base in your input sequence
and check if PTO after the base is required.

Examples of supported converters:

| CMO Codes |  Example of ASO Gapmers  |
|---|---|
| Raw | AGGTCTTCATGACTTCGGCC |
| BioSpring | 6\*8\*8\*5\*7\*T\*T\*9\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |
| Axolabs | - |
| GCRS | moeAnpsmoeGnpsmoeGnpsmoeUnpsmoe5mCnpsTpsTpsCpsApsTpsGpsApsCpsTpsTpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC |

| CMO Codes |  Example of 2'-OMe and 2'-F modified siRNA  |
|---|---|
| Raw | UUCAACUGCUUACGUCUUU |
| BioSpring | 5\*1\*766354715274575\*5\*5 |
| Axolabs | usUfscaaCfuGfcUfuAfcGfucususu |
| GCRS | mUpsfUpsmCmAmAfCmUfGmCfUmUfAmCfGmUmCmUpsmUpsmU |

| Function | Input | Output | 
|---|---|---|
| classicToBioSpring | TTGTCCAGATGACTTCGGCC | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |
| classicToGCRS | TTGTCCAGATGACTTCGGCC | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | 
| bioSpringToClassic | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 | TTGTCCAGATGACTTCGGCC |
| bioSpringToGCRS | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC |
| GCRSToBioSpring | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |
| GCRSToClassic | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | TTGTCCAGATGACTTCGGCC | 

Tables of codes accordance:

| CMO Codes For ASO Gapmers | BioSpring Code | Axolabs Code | GCRS code |
|---|---|---|---|
| 2'MOE-5Me-rU | 5 | | moeT |
| 2'MOE-rA     | 6 | | moeA |
| 2'MOE-5Me-rC | 7 | | moe5mC |
| 2'MOE-rG     | 8 | | moeG |
| 5-Methyl-dC  | 9 | | 5mC |
| ps linkage   | * | | ps |
| dA           | A | | A |
| dC           | C | | C |
| dG           | G | | T |
| dT           | T | | G |

| CMO Codes For 2'-OMe and 2'-F modified siRNA  | BioSpring Code | Axolabs Code | GCRS code |
|---|---|---|---|
| 2'-fluoro-U | 1 | Uf | fU |
| 2'-fluoro-A | 2 | Af | fA |
| 2'-fluoro-C | 3 | Cf | fC |
| 2'-fluoro-G | 4 | Gf | fG |
| 2'OMe-rU    | 5 | u  | mU |
| 2'OMe-rA    | 6 | a  | mA |
| 2'OMe-rC    | 7 | c  | mC |
| 2'OMe-rG    | 8 | g  | mG |
| ps linkage  | * | s  | ps |

# Representations / synthesizers