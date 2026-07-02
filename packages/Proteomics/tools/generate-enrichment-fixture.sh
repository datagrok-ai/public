#!/bin/bash
# Engineered human Treatment-vs-Control DE result for the enrichment demo.
# UP = cell cycle / mitosis, DOWN = oxidative phosphorylation, plus a diverse
# non-significant background. Gene symbols are real HGNC; accessions real UniProt.
out=files/demo/enrichment-demo.csv
echo "Protein ID,Gene,log2FC,p-value,adj.p-value,significant" > "$out"

UP="CDK1 P06493
CCNB1 P14635
CCNB2 O95067
CDC20 Q12834
BUB1 O43683
BUB1B O60566
AURKA O14965
AURKB Q96GD4
PLK1 P53350
CCNA2 P20248
CDK2 P24941
MAD2L1 Q13257
TOP2A P11388
KIF11 P52732
CENPA P49450
CENPE Q02224
NDC80 O14777
CDC6 Q99741
CDC45 O75419
MCM2 P49736
MCM3 P25205
MCM4 P33991
MCM5 P33992
MCM6 Q14566
MCM7 P33993
CDC25C P30307
CDK4 P11802
CCNE1 P24864
CHEK1 O14757
WEE1 P30291"

DOWN="NDUFA9 Q16795
NDUFS1 P28331
NDUFS2 O75306
NDUFB10 O96000
NDUFV1 P49821
NDUFA4 O00483
NDUFB8 O95169
SDHA P31040
SDHB P21912
UQCRC1 P31930
UQCRC2 P22695
UQCRFS1 P47985
CYC1 P08574
COX4I1 P13073
COX5A P20674
COX5B P10606
COX6C P09669
COX7A2 P14406
COX6B1 P14854
ATP5F1A P25705
ATP5F1B P06576
ATP5F1C P36542
ATP5PO P48047
ATP5PB P24539
ATP5MC1 P05496"

BG="GAPDH P04406
ACTB P60709
TUBB P07437
ALB P02768
ENO1 P06733
PKM P14618
LDHA P00338
PGK1 P00558
ALDOA P04075
TPI1 P60174
GPI P06744
PGAM1 P18669
HSPA8 P11142
HSP90AA1 P07900
HSPA5 P11021
HSPD1 P10809
CALR P27797
CANX P27824
PDIA3 P30101
P4HB P07237
VIM P08670
FLNA P21333
ACTN1 P12814
TLN1 Q9Y490
CFL1 P23528
PFN1 P07737
GSN P06396
CAP1 Q01518
ANXA1 P04083
ANXA2 P07355
ANXA5 P08758
S100A4 P26447
YWHAZ P63104
YWHAE P62258
YWHAB P31946
PPIA P62937
PPIB P23284
SOD1 P00441
SOD2 P04179
CAT P04040
PRDX1 Q06830
PRDX2 P32119
GSTP1 P09211
TXN P10599
FTL P02792
FTH1 P02794
TF P02787
CTSD P07339
CTSB P07858
LGALS1 P09382
LGALS3 P17931
CD44 P16070
ITGB1 P05556
EEF1A1 P68104
EEF2 P13639
RPS3 P23396
RPL4 P36578
RPLP0 P05388
NPM1 P06748
NCL P19338
HNRNPK P61978
SRSF1 Q07955
G3BP1 Q13283
DDX5 P17844"

i=0
while read -r gene acc; do
  [ -z "$gene" ] && continue
  fc=$(awk -v i=$i 'BEGIN{printf "%.2f", 1.5 + (i%10)*0.18}')
  pexp=$(( 3 + i%4 ))
  p=$(awk -v e=$pexp 'BEGIN{printf "%.2e", 10^(-e)}')
  adj=$(awk -v e=$pexp 'BEGIN{printf "%.2e", 2*10^(-e)}')
  echo "$acc,$gene,$fc,$p,$adj,true" >> "$out"
  i=$((i+1))
done <<< "$UP"

i=0
while read -r gene acc; do
  [ -z "$gene" ] && continue
  fc=$(awk -v i=$i 'BEGIN{printf "%.2f", -(1.5 + (i%8)*0.18)}')
  pexp=$(( 3 + i%4 ))
  p=$(awk -v e=$pexp 'BEGIN{printf "%.2e", 10^(-e)}')
  adj=$(awk -v e=$pexp 'BEGIN{printf "%.2e", 2*10^(-e)}')
  echo "$acc,$gene,$fc,$p,$adj,true" >> "$out"
  i=$((i+1))
done <<< "$DOWN"

i=0
while read -r gene acc; do
  [ -z "$gene" ] && continue
  fc=$(awk -v i=$i 'BEGIN{printf "%.2f", -0.6 + (i%13)*0.09}')
  p=$(awk -v i=$i 'BEGIN{printf "%.2f", 0.10 + (i%9)*0.09}')
  adj=$(awk -v i=$i 'BEGIN{v=(0.10 + (i%9)*0.09)*1.05; if(v>0.99)v=0.99; printf "%.2f", v}')
  echo "$acc,$gene,$fc,$p,$adj,false" >> "$out"
  i=$((i+1))
done <<< "$BG"

echo "=== rows: $(($(wc -l < "$out") - 1)) (up=$(grep -c ',true$' "$out") sig, bg=$(grep -c ',false$' "$out") ns) ==="
echo "=== head ==="; head -4 "$out"
echo "=== down sample ==="; grep -E ',-' "$out" | head -3
echo "=== bg sample ==="; grep ',false$' "$out" | head -3
