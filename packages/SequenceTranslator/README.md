# SequenceTranslator

SequenceTranslator is a [package](https://datagrok.ai/help/develop/develop#packages) for
the [Datagrok](https://datagrok.ai) platform.

Examples of supported converters:

| Nucleotide type | Classic Code | BioSpring Code | Axolabs Code | Janssen GCRS code |
|---|---|---|---|---|
| ASO Gapmers                    | TTGTCCAGATGACTTCGGCC | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |-|moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC |
| 2'-OMe and 2'-F modified siRNA | GUUUAAUCAGAAGAGGAUU | 8\*5\*15261728264648252\*5\*5 | gsusUfuAfaUfcAfgAfaGfaGfgAfuAfsusu | mGpsmUpsfUmUfAmAfUmCfAmGfAmAfGmAfGmGfAmUfApsmUpsmU |

| Function | Input | Output | 
|---|---|---|
| classicToBioSpring | TTGTCCAGATGACTTCGGCC | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |
| classicToGCRS | TTGTCCAGATGACTTCGGCC | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | 
| bioSpringToClassic | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 | TTGTCCAGATGACTTCGGCC |
| bioSpringToGCRS | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC |
| GCRSToBioSpring | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | 5\*5\*8\*5\*7\*9\*A\*G\*A\*T\*G\*A\*9\*T\*T\*7\*8\*8\*7\*7 |
| GCRSToClassic | moeUnpsmoeUnpsmoeGnpsmoeUnpsmoe5mCps5mCpsApsGpsApsTpsGpsAps5mCpsTpsTnpsmoe5mCnpsmoeGnpsmoeGnpsmoe5mCnpsmoe5mC | TTGTCCAGATGACTTCGGCC | 

Tables of codes accordance:

| CMO Codes For ASO Gapmers | BioSpring Code | Axolabs Code | Janssen GCRS code |
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

| CMO Codes For 2'-OMe and 2'-F modified siRNA  | BioSpring Code | Axolabs Code | Janssen GCRS code |
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