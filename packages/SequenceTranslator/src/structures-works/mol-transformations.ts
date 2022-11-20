const PHOSHATE = `
Datagrok monomer library Nucleotides

  0  0  0  0  0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -1.5 0 0 0
M  V30 2 P 0 0 0 0
M  V30 3 O 0 1 0 0
M  V30 4 O 0 -1 0 0
M  V30 5 O 1.5 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  V30 BEGIN COLLECTION
M  V30 END COLLECTION
M  END`;

const THIOPHOSHATE = `
Datagrok monomer library Nucleotides

  0  0  0  0  0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -1.5 0 0 0
M  V30 2 P 0 0 0 0
M  V30 3 O 0 1 0 0
M  V30 4 S 0 -1 0 0
M  V30 5 O 1.5 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  V30 BEGIN COLLECTION
M  V30 END COLLECTION
M  END`;

const INVABASIC = `
Datagrok monomer library Nucleotides

  0  0  0  0  0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 1.0934 -2.1636 0 0
M  V30 2 C 1.8365 -1.4945 0 0 CFG=2
M  V30 3 C 2.8147 -1.7024 0 0
M  V30 4 C 3.3147 -0.8364 0 0
M  V30 5 O 2.6455 -0.0932 0 0
M  V30 6 C 1.732 -0.5 0 0 CFG=1
M  V30 7 C 0.866 0 0 0
M  V30 8 O 0.866 1 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 6 5
M  V30 6 1 2 6
M  V30 7 1 6 7 CFG=3
M  V30 8 1 7 8
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 2 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  END`;

const GALNAC = `
Datagrok monomer library Nucleotides                            

0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 111 113 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -20.7313 -0.7027 0 0
M  V30 2 C -19.3976 0.0673 0 0
M  V30 3 C -18.0638 -0.7027 0 0
M  V30 4 C -16.7303 0.0673 0 0
M  V30 5 N -15.3965 -0.7027 0 0
M  V30 6 C -14.0628 0.0673 0 0
M  V30 7 C -12.7293 -0.7027 0 0
M  V30 8 C -11.3955 0.0673 0 0
M  V30 9 C -10.062 -0.7027 0 0
M  V30 10 C -8.7283 0.0673 0 0
M  V30 11 N -7.3947 -0.7027 0 0
M  V30 12 O -18.0638 -2.2427 0 0
M  V30 13 O -14.0628 1.6073 0 0
M  V30 14 O -8.7283 1.6073 0 0
M  V30 15 C -5.8547 -0.7027 0 0
M  V30 16 C -5.8547 0.8373 0 0
M  V30 17 C -5.8547 -2.2427 0 0
M  V30 18 C -3.4848 -3.0127 0 0
M  V30 19 C -2.4544 -4.157 0 0
M  V30 20 C -0.948 -3.8368 0 0
M  V30 21 N 0.0824 -4.9813 0 0
M  V30 22 C 1.5888 -4.6612 0 0
M  V30 23 C 2.6192 -5.8056 0 0
M  V30 24 C 4.1256 -5.4855 0 0
M  V30 25 N 5.156 -6.6297 0 0
M  V30 26 C 6.6624 -6.3096 0 0
M  V30 27 C 7.6928 -7.4541 0 0
M  V30 28 C 9.1992 -7.1339 0 0
M  V30 29 C 10.2296 -8.2784 0 0
M  V30 30 C 11.736 -7.9583 0 0
M  V30 31 O 12.7664 -9.1027 0 0
M  V30 32 O -0.4722 -2.3723 0 0
M  V30 33 O 7.1382 -4.845 0 0
M  V30 34 C 14.2728 -8.7824 0 0
M  V30 35 C 15.3032 -9.9267 0 0
M  V30 36 C 16.8098 -9.6065 0 0
M  V30 37 C 17.2856 -8.1421 0 0
M  V30 38 C 16.2552 -6.9975 0 0
M  V30 39 O 14.7486 -7.3178 0 0
M  V30 40 C 16.7312 -5.5329 0 0
M  V30 41 O 18.7918 -7.8218 0 0
M  V30 42 O 17.8404 -10.751 0 0
M  V30 43 N 14.8274 -11.3914 0 0
M  V30 44 C 15.7325 -12.6372 0 0
M  V30 45 C 15.2567 -14.1018 0 0
M  V30 46 O 17.2537 -12.3963 0 0
M  V30 47 O 18.2628 -5.372 0 0
M  V30 48 O -4.9494 -3.4885 0 0
M  V30 49 C -4.521 0.0673 0 0
M  V30 50 C -1.9414 0.2026 0 0
M  V30 51 C -0.6077 -0.5674 0 0
M  V30 52 C 0.726 0.2026 0 0
M  V30 53 N 2.0596 -0.5674 0 0
M  V30 54 C 3.3933 0.2026 0 0
M  V30 55 C 4.7271 -0.5674 0 0
M  V30 56 C 6.0606 0.2026 0 0
M  V30 57 N 7.3943 -0.5674 0 0
M  V30 58 C 8.7281 0.2026 0 0
M  V30 59 C 10.0618 -0.5674 0 0
M  V30 60 C 11.3953 0.2026 0 0
M  V30 61 C 14.0628 0.2026 0 0
M  V30 62 O 15.3964 -0.5674 0 0
M  V30 63 O 0.726 1.7426 0 0
M  V30 64 O 8.7281 1.7426 0 0
M  V30 65 C 16.7301 0.2026 0 0
M  V30 66 C 18.0638 -0.5676 0 0
M  V30 67 C 19.3976 0.2026 0 0
M  V30 68 C 19.3974 1.7426 0 0
M  V30 69 C 18.0638 2.5126 0 0
M  V30 70 O 16.7301 1.7426 0 0
M  V30 71 C 18.064 4.0526 0 0
M  V30 72 O 20.7311 2.5126 0 0
M  V30 73 O 20.7313 -0.5674 0 0
M  V30 74 N 18.0638 -2.1076 0 0
M  V30 75 C 19.3096 -3.0127 0 0
M  V30 76 C 19.3096 -4.5527 0 0
M  V30 77 O 20.6818 -2.3135 0 0
M  V30 78 O 19.4709 4.6791 0 0
M  V30 79 O -3.1872 -0.7027 0 0
M  V30 80 C 12.7291 -0.5674 0 0
M  V30 81 C -3.919 3.2277 0 0
M  V30 82 C -2.4126 2.9076 0 0
M  V30 83 C -1.3822 4.0519 0 0
M  V30 84 N 0.1242 3.7317 0 0
M  V30 85 C 1.1546 4.8762 0 0
M  V30 86 C 2.661 4.5561 0 0
M  V30 87 C 3.6914 5.7005 0 0
M  V30 88 N 5.1978 5.3804 0 0
M  V30 89 C 6.2282 6.5248 0 0
M  V30 90 C 7.7346 6.2045 0 0
M  V30 91 C 8.765 7.349 0 0
M  V30 92 C 10.2714 7.0288 0 0
M  V30 93 C 11.3018 8.1733 0 0
M  V30 94 O 12.8082 7.8532 0 0
M  V30 95 O -1.858 5.5167 0 0
M  V30 96 O 5.7524 7.9894 0 0
M  V30 97 C 13.8386 8.9976 0 0
M  V30 98 C 15.345 8.6773 0 0
M  V30 99 C 16.3756 9.8219 0 0
M  V30 100 C 15.8996 11.2863 0 0
M  V30 101 C 14.3934 11.6068 0 0
M  V30 102 O 13.3628 10.4622 0 0
M  V30 103 C 13.9176 13.0714 0 0
M  V30 104 O 16.93 12.4308 0 0
M  V30 105 O 17.882 9.5018 0 0
M  V30 106 N 15.8208 7.2127 0 0
M  V30 107 C 17.2856 6.7367 0 0
M  V30 108 C 17.7614 5.2721 0 0
M  V30 109 O 18.3744 7.8257 0 0
M  V30 110 O 15.062 14.1018 0 0
M  V30 111 O -4.8241 1.9817 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 3 4
M  V30 3 1 6 7
M  V30 4 1 7 8
M  V30 5 1 8 9
M  V30 6 1 9 10
M  V30 7 1 1 2
M  V30 8 1 3 12
M  V30 9 1 4 5
M  V30 10 1 5 6
M  V30 11 2 6 13
M  V30 12 1 10 11
M  V30 13 1 11 15
M  V30 14 1 15 16
M  V30 15 1 15 17
M  V30 16 2 10 14
M  V30 17 1 18 19
M  V30 18 1 19 20
M  V30 19 1 22 23
M  V30 20 1 23 24
M  V30 21 1 26 27
M  V30 22 1 27 28
M  V30 23 1 28 29
M  V30 24 1 29 30
M  V30 25 2 26 33
M  V30 26 2 20 32
M  V30 27 1 20 21
M  V30 28 1 21 22
M  V30 29 1 24 25
M  V30 30 1 25 26
M  V30 31 1 30 31
M  V30 32 1 31 34
M  V30 33 1 35 36
M  V30 34 1 36 37
M  V30 35 1 37 38
M  V30 36 1 34 35
M  V30 37 1 38 39
M  V30 38 1 34 39
M  V30 39 1 38 40
M  V30 40 1 35 43
M  V30 41 1 43 44
M  V30 42 1 44 45
M  V30 43 2 44 46
M  V30 44 1 36 42
M  V30 45 1 37 41
M  V30 46 1 40 47
M  V30 47 1 18 48
M  V30 48 1 15 49
M  V30 49 1 50 51
M  V30 50 1 51 52
M  V30 51 1 54 55
M  V30 52 1 55 56
M  V30 53 1 58 59
M  V30 54 1 59 60
M  V30 55 2 58 64
M  V30 56 2 52 63
M  V30 57 1 52 53
M  V30 58 1 53 54
M  V30 59 1 56 57
M  V30 60 1 57 58
M  V30 61 1 61 62
M  V30 62 1 62 65
M  V30 63 1 66 67
M  V30 64 1 67 68
M  V30 65 1 68 69
M  V30 66 1 65 66
M  V30 67 1 69 70
M  V30 68 1 65 70
M  V30 69 1 69 71
M  V30 70 1 66 74
M  V30 71 1 74 75
M  V30 72 1 75 76
M  V30 73 2 75 77
M  V30 74 1 67 73
M  V30 75 1 68 72
M  V30 76 1 71 78
M  V30 77 1 50 79
M  V30 78 1 49 79
M  V30 79 1 60 80
M  V30 80 1 80 61
M  V30 81 1 81 82
M  V30 82 1 82 83
M  V30 83 1 85 86
M  V30 84 1 86 87
M  V30 85 1 89 90
M  V30 86 1 90 91
M  V30 87 1 91 92
M  V30 88 1 92 93
M  V30 89 2 89 96
M  V30 90 2 83 95
M  V30 91 1 83 84
M  V30 92 1 84 85
M  V30 93 1 87 88
M  V30 94 1 88 89
M  V30 95 1 93 94
M  V30 96 1 94 97
M  V30 97 1 98 99
M  V30 98 1 99 100
M  V30 99 1 100 101
M  V30 100 1 97 98
M  V30 101 1 101 102
M  V30 102 1 97 102
M  V30 103 1 101 103
M  V30 104 1 98 106
M  V30 105 1 106 107
M  V30 106 1 107 108
M  V30 107 2 107 109
M  V30 108 1 99 105
M  V30 109 1 100 104
M  V30 110 1 103 110
M  V30 111 1 81 111
M  V30 112 1 16 111
M  V30 113 1 17 48
M  V30 END BOND
M  V30 END CTAB
M  END`;

const GALNACPRIME = `
Datagrok monomer library Nucleotides          

0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 111 113 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 20.7313 0.7027 0 0
M  V30 2 C 19.3976 -0.0673 0 0
M  V30 3 C 18.0638 0.7027 0 0
M  V30 4 C 16.7303 -0.0673 0 0
M  V30 5 N 15.3965 0.7027 0 0
M  V30 6 C 14.0628 -0.0673 0 0
M  V30 7 C 12.7293 0.7027 0 0
M  V30 8 C 11.3955 -0.0673 0 0
M  V30 9 C 10.062 0.7027 0 0
M  V30 10 C 8.7283 -0.0673 0 0
M  V30 11 N 7.3947 0.7027 0 0
M  V30 12 O 18.0638 2.2427 0 0
M  V30 13 O 14.0628 -1.6073 0 0
M  V30 14 O 8.7283 -1.6073 0 0
M  V30 15 C 5.8547 0.7027 0 0
M  V30 16 C 5.8547 -0.8373 0 0
M  V30 17 C 5.8547 2.2427 0 0
M  V30 18 C 3.4848 3.0127 0 0
M  V30 19 C 2.4544 4.157 0 0
M  V30 20 C 0.948 3.8368 0 0
M  V30 21 N -0.0824 4.9813 0 0
M  V30 22 C -1.5888 4.6612 0 0
M  V30 23 C -2.6192 5.8056 0 0
M  V30 24 C -4.1256 5.4855 0 0
M  V30 25 N -5.156 6.6297 0 0
M  V30 26 C -6.6624 6.3096 0 0
M  V30 27 C -7.6928 7.4541 0 0
M  V30 28 C -9.1992 7.1339 0 0
M  V30 29 C -10.2296 8.2784 0 0
M  V30 30 C -11.736 7.9583 0 0
M  V30 31 O -12.7664 9.1027 0 0
M  V30 32 O 0.4722 2.3723 0 0
M  V30 33 O -7.1382 4.845 0 0
M  V30 34 C -14.2728 8.7824 0 0
M  V30 35 C -15.3032 9.9267 0 0
M  V30 36 C -16.8098 9.6065 0 0
M  V30 37 C -17.2856 8.1421 0 0
M  V30 38 C -16.2552 6.9975 0 0
M  V30 39 O -14.7486 7.3178 0 0
M  V30 40 C -16.7312 5.5329 0 0
M  V30 41 O -18.7918 7.8218 0 0
M  V30 42 O -17.8404 10.751 0 0
M  V30 43 N -14.8274 11.3914 0 0
M  V30 44 C -15.7325 12.6372 0 0
M  V30 45 C -15.2567 14.1018 0 0
M  V30 46 O -17.2537 12.3963 0 0
M  V30 47 O -18.2628 5.372 0 0
M  V30 48 O 4.9494 3.4885 0 0
M  V30 49 C 4.521 -0.0673 0 0
M  V30 50 C 1.9414 -0.2026 0 0
M  V30 51 C 0.6077 0.5674 0 0
M  V30 52 C -0.726 -0.2026 0 0
M  V30 53 N -2.0596 0.5674 0 0
M  V30 54 C -3.3933 -0.2026 0 0
M  V30 55 C -4.7271 0.5674 0 0
M  V30 56 C -6.0606 -0.2026 0 0
M  V30 57 N -7.3943 0.5674 0 0
M  V30 58 C -8.7281 -0.2026 0 0
M  V30 59 C -10.0618 0.5674 0 0
M  V30 60 C -11.3953 -0.2026 0 0
M  V30 61 C -14.0628 -0.2026 0 0
M  V30 62 O -15.3964 0.5674 0 0
M  V30 63 O -0.726 -1.7426 0 0
M  V30 64 O -8.7281 -1.7426 0 0
M  V30 65 C -16.7301 -0.2026 0 0
M  V30 66 C -18.0638 0.5676 0 0
M  V30 67 C -19.3976 -0.2026 0 0
M  V30 68 C -19.3974 -1.7426 0 0
M  V30 69 C -18.0638 -2.5126 0 0
M  V30 70 O -16.7301 -1.7426 0 0
M  V30 71 C -18.064 -4.0526 0 0
M  V30 72 O -20.7311 -2.5126 0 0
M  V30 73 O -20.7313 0.5674 0 0
M  V30 74 N -18.0638 2.1076 0 0
M  V30 75 C -19.3096 3.0127 0 0
M  V30 76 C -19.3096 4.5527 0 0
M  V30 77 O -20.6818 2.3135 0 0
M  V30 78 O -19.4709 -4.6791 0 0
M  V30 79 O 3.1872 0.7027 0 0
M  V30 80 C -12.7291 0.5674 0 0
M  V30 81 C 3.919 -3.2277 0 0
M  V30 82 C 2.4126 -2.9076 0 0
M  V30 83 C 1.3822 -4.0519 0 0
M  V30 84 N -0.1242 -3.7317 0 0
M  V30 85 C -1.1546 -4.8762 0 0
M  V30 86 C -2.661 -4.5561 0 0
M  V30 87 C -3.6914 -5.7005 0 0
M  V30 88 N -5.1978 -5.3804 0 0
M  V30 89 C -6.2282 -6.5248 0 0
M  V30 90 C -7.7346 -6.2045 0 0
M  V30 91 C -8.765 -7.349 0 0
M  V30 92 C -10.2714 -7.0288 0 0
M  V30 93 C -11.3018 -8.1733 0 0
M  V30 94 O -12.8082 -7.8532 0 0
M  V30 95 O 1.858 -5.5167 0 0
M  V30 96 O -5.7524 -7.9894 0 0
M  V30 97 C -13.8386 -8.9976 0 0
M  V30 98 C -15.345 -8.6773 0 0
M  V30 99 C -16.3756 -9.8219 0 0
M  V30 100 C -15.8996 -11.2863 0 0
M  V30 101 C -14.3934 -11.6068 0 0
M  V30 102 O -13.3628 -10.4622 0 0
M  V30 103 C -13.9176 -13.0714 0 0
M  V30 104 O -16.93 -12.4308 0 0
M  V30 105 O -17.882 -9.5018 0 0
M  V30 106 N -15.8208 -7.2127 0 0
M  V30 107 C -17.2856 -6.7367 0 0
M  V30 108 C -17.7614 -5.2721 0 0
M  V30 109 O -18.3744 -7.8257 0 0
M  V30 110 O -15.062 -14.1018 0 0
M  V30 111 O 4.8241 -1.9817 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 3 4
M  V30 3 1 6 7
M  V30 4 1 7 8
M  V30 5 1 8 9
M  V30 6 1 9 10
M  V30 7 1 1 2
M  V30 8 1 3 12
M  V30 9 1 4 5
M  V30 10 1 5 6
M  V30 11 2 6 13
M  V30 12 1 10 11
M  V30 13 1 11 15
M  V30 14 1 15 16
M  V30 15 1 15 17
M  V30 16 2 10 14
M  V30 17 1 18 19
M  V30 18 1 19 20
M  V30 19 1 22 23
M  V30 20 1 23 24
M  V30 21 1 26 27
M  V30 22 1 27 28
M  V30 23 1 28 29
M  V30 24 1 29 30
M  V30 25 2 26 33
M  V30 26 2 20 32
M  V30 27 1 20 21
M  V30 28 1 21 22
M  V30 29 1 24 25
M  V30 30 1 25 26
M  V30 31 1 30 31
M  V30 32 1 31 34
M  V30 33 1 35 36
M  V30 34 1 36 37
M  V30 35 1 37 38
M  V30 36 1 34 35
M  V30 37 1 38 39
M  V30 38 1 34 39
M  V30 39 1 38 40
M  V30 40 1 35 43
M  V30 41 1 43 44
M  V30 42 1 44 45
M  V30 43 2 44 46
M  V30 44 1 36 42
M  V30 45 1 37 41
M  V30 46 1 40 47
M  V30 47 1 18 48
M  V30 48 1 15 49
M  V30 49 1 50 51
M  V30 50 1 51 52
M  V30 51 1 54 55
M  V30 52 1 55 56
M  V30 53 1 58 59
M  V30 54 1 59 60
M  V30 55 2 58 64
M  V30 56 2 52 63
M  V30 57 1 52 53
M  V30 58 1 53 54
M  V30 59 1 56 57
M  V30 60 1 57 58
M  V30 61 1 61 62
M  V30 62 1 62 65
M  V30 63 1 66 67
M  V30 64 1 67 68
M  V30 65 1 68 69
M  V30 66 1 65 66
M  V30 67 1 69 70
M  V30 68 1 65 70
M  V30 69 1 69 71
M  V30 70 1 66 74
M  V30 71 1 74 75
M  V30 72 1 75 76
M  V30 73 2 75 77
M  V30 74 1 67 73
M  V30 75 1 68 72
M  V30 76 1 71 78
M  V30 77 1 50 79
M  V30 78 1 49 79
M  V30 79 1 60 80
M  V30 80 1 80 61
M  V30 81 1 81 82
M  V30 82 1 82 83
M  V30 83 1 85 86
M  V30 84 1 86 87
M  V30 85 1 89 90
M  V30 86 1 90 91
M  V30 87 1 91 92
M  V30 88 1 92 93
M  V30 89 2 89 96
M  V30 90 2 83 95
M  V30 91 1 83 84
M  V30 92 1 84 85
M  V30 93 1 87 88
M  V30 94 1 88 89
M  V30 95 1 93 94
M  V30 96 1 94 97
M  V30 97 1 98 99
M  V30 98 1 99 100
M  V30 99 1 100 101
M  V30 100 1 97 98
M  V30 101 1 101 102
M  V30 102 1 97 102
M  V30 103 1 101 103
M  V30 104 1 98 106
M  V30 105 1 106 107
M  V30 106 1 107 108
M  V30 107 2 107 109
M  V30 108 1 99 105
M  V30 109 1 100 104
M  V30 110 1 103 110
M  V30 111 1 81 111
M  V30 112 1 16 111
M  V30 113 1 17 48
M  V30 END BOND
M  V30 END CTAB
M  END`;

export function getNucleotidesMol(smilesCodes: string[]) {
  const molBlocks: string[] = [];

  for (let i = 0; i < smilesCodes.length - 1; i++) {
    smilesCodes[i] == 'OP(=O)(O)O' ? molBlocks.push(PHOSHATE) :
      smilesCodes[i] == 'OP(=O)(S)O' ? molBlocks.push(THIOPHOSHATE) :
        smilesCodes[i] == 'O[C@@H]1C[C@@H]O[C@H]1CO' ? molBlocks.push(rotateNucleotidesV3000(INVABASIC)) :
          smilesCodes[i] == 'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' ? molBlocks.push(GALNAC) :
            smilesCodes[i] == 'C(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)NC(=O)CCCC(=O)NCC(O)CO' ? molBlocks.push(GALNACPRIME) :
          molBlocks.push(rotateNucleotidesV3000(smilesCodes[i]));
  }

  return linkV3000(molBlocks);
}

export function linkStrandsV3000(strands:{senseStrands: string[], antiStrands: string[]}, useChirality: boolean = true) {
  let macroMolBlock = '\nDatagrok macromolecule handler\n\n';
  macroMolBlock += '  0  0  0  0  0  0            999 V3000\n';
  macroMolBlock += 'M  V30 BEGIN CTAB\n';
  let atomBlock = '';
  let bondBlock = '';
  let collectionBlock = '';
  const collection: number [] = [];
  let natom = 0;
  let nbond = 0;
  let xShift = 0;

  // if (twoChains && molBlocks.length > 1)
  //   molBlocks[1] = invertNucleotidesV3000(molBlocks[1]);

  if (strands.antiStrands.length > 0) {
    for(let i = 0; i < strands.antiStrands.length; i++) {
      strands.antiStrands[i] = invertNucleotidesV3000(strands.antiStrands[i]);
    }
  }

  let inverted = false;
  let molBlocks = strands.senseStrands.concat(strands.antiStrands);

  for (let i = 0; i < molBlocks.length; i++) {

    if (i >= strands.senseStrands.length && inverted == false) {
      inverted = true;
      xShift = 0;
    }


    molBlocks[i] = molBlocks[i].replaceAll('(-\nM  V30 ', '(')
      .replaceAll('-\nM  V30 ', '').replaceAll(' )', ')');
    const numbers = extractAtomsBondsNumbersV3000(molBlocks[i]);
    const coordinates = extractAtomDataV3000(molBlocks[i]);

    if (inverted) {
      const xShiftRight = Math.min(...coordinates.x) - xShift;
      const yShift = !inverted ? Math.min(...coordinates.y) - 1 : Math.max(...coordinates.y) + 10;
      for (let j = 0; j < coordinates.x.length; j++)
        coordinates.x[j] -= xShiftRight;
      for (let j = 0; j < coordinates.y.length; j++)
        coordinates.y[j] -= yShift;
    }

    let indexAtoms = molBlocks[i].indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
    indexAtoms = molBlocks[i].indexOf('\n', indexAtoms);
    let index = indexAtoms;
    let indexEnd = indexAtoms;

    for (let j = 0; j < numbers.natom; j++) {
      // if (coordinates.atomIndex[j] != 1 || i == 0 || twoChains) {
      //rewrite atom number
      index = molBlocks[i].indexOf('V30', index) + 4;
      indexEnd = molBlocks[i].indexOf(' ', index);
      const atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

      //rewrite coordinates
      index = molBlocks[i].indexOf(' ', index) + 1;
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf(' ', index);

      const totalShift = true ? 0 : xShift - coordinates.x[0];
      let coordinate = true ?
        Math.round(10000*coordinates.x[j])/10000 :
        Math.round(10000*(parseFloat(molBlocks[i].substring(index, indexEnd)) + totalShift))/10000;
      molBlocks[i] = molBlocks[i].slice(0, index) + coordinate + molBlocks[i].slice(indexEnd);

      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf(' ', index);
      coordinate = true ?
        Math.round(10000*coordinates.y[j])/10000 :
        Math.round(10000*(parseFloat(molBlocks[i].substring(index, indexEnd))))/10000;
      molBlocks[i] = molBlocks[i].slice(0, index) + coordinate + molBlocks[i].slice(indexEnd);

      index = molBlocks[i].indexOf('\n', index) + 1;
    }

    const indexAtomsEnd = molBlocks[i].indexOf('M  V30 END ATOM');
    atomBlock += molBlocks[i].substring(indexAtoms + 1, indexAtomsEnd);

    let indexBonds = molBlocks[i].indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
    indexBonds = molBlocks[i].indexOf('\n', indexBonds);
    index = indexBonds;
    indexEnd = indexBonds;

    for (let j = 0; j < numbers.nbond; j++) {
      //rewrite bond number
      index = molBlocks[i].indexOf('V30', index) + 4;
      indexEnd = molBlocks[i].indexOf(' ', index);
      const bondNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + nbond;
      molBlocks[i] = molBlocks[i].slice(0, index) + bondNumber + molBlocks[i].slice(indexEnd);

      //rewrite atom pair in bond
      index = molBlocks[i].indexOf(' ', index) + 1;
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf(' ', index);
      let atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = Math.min(molBlocks[i].indexOf('\n', index), molBlocks[i].indexOf(' ', index));
      atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

      index = molBlocks[i].indexOf('\n', index) + 1;
    }

    const indexBondEnd = molBlocks[i].indexOf('M  V30 END BOND');
    bondBlock += molBlocks[i].substring(indexBonds + 1, indexBondEnd);

    let indexCollection = molBlocks[i].indexOf('M  V30 MDLV30/STEABS ATOMS=('); // V3000 index for collections

    while (indexCollection != -1) {
      indexCollection += 28;
      const collectionEnd = molBlocks[i].indexOf(')', indexCollection);
      const collectionEntries = molBlocks[i].substring(indexCollection, collectionEnd).split(' ').slice(1);
      collectionEntries.forEach((e) => {
        collection.push(parseInt(e) + natom);
      });
      indexCollection = collectionEnd;
      indexCollection = molBlocks[i].indexOf('M  V30 MDLV30/STEABS ATOMS=(', indexCollection);
    }

    natom += true ? numbers.natom : numbers.natom - 1;
    nbond += numbers.nbond;
    xShift += Math.max(...coordinates.x) + 1;//twoChains ? 0 : coordinates.x[numbers.natom - 1] - coordinates.x[0];
  }

  const entries = 4;
  const collNumber = Math.ceil(collection.length / entries);

  //if (oclRender) {
    // collectionBlock += 'M  V30 MDLV30/STEABS ATOMS=(' + collection.length;

    // for (let j = 0; j < collection.length; j++)
    //   collectionBlock += ' ' + collection[j];

    // collectionBlock += ')\n';
  //} else {
    collectionBlock += 'M  V30 MDLV30/STEABS ATOMS=(' + collection.length + ' -\n';
    for (let i = 0; i < collNumber; i++) {
      collectionBlock += 'M  V30 ';
      const entriesCurrent = i + 1 == collNumber ? collection.length - (collNumber - 1)*entries : entries;
      for (let j = 0; j < entriesCurrent; j++) {
        collectionBlock += (j + 1 == entriesCurrent) ?
          (i == collNumber - 1 ? collection[entries*i + j] + ')\n' : collection[entries*i + j] + ' -\n') :
          collection[entries*i + j] + ' ';
      }
    }
  //}

  //generate file
  true? natom : natom++;
  macroMolBlock += 'M  V30 COUNTS ' + natom + ' ' + nbond + ' 0 0 0\n';
  macroMolBlock += 'M  V30 BEGIN ATOM\n';
  macroMolBlock += atomBlock;
  macroMolBlock += 'M  V30 END ATOM\n';
  macroMolBlock += 'M  V30 BEGIN BOND\n';
  macroMolBlock += bondBlock;
  macroMolBlock += 'M  V30 END BOND\n';
  if(useChirality){
    macroMolBlock += 'M  V30 BEGIN COLLECTION\n';
    macroMolBlock += collectionBlock;
    macroMolBlock += 'M  V30 END COLLECTION\n';
  } else
    macroMolBlock = macroMolBlock.replace(/ CFG=\d/g, ' ');

  macroMolBlock += 'M  V30 END CTAB\n';
  macroMolBlock += 'M  END';

  return macroMolBlock;
}

export function linkV3000(molBlocks: string[], useChirality: boolean = true) {
  let macroMolBlock = '\nDatagrok macromolecule handler\n\n';
  macroMolBlock += '  0  0  0  0  0  0            999 V3000\n';
  macroMolBlock += 'M  V30 BEGIN CTAB\n';
  let atomBlock = '';
  let bondBlock = '';
  let collectionBlock = '';
  const collection: number [] = [];
  let natom = 0;
  let nbond = 0;
  let xShift = 0;

  for (let i = 0; i < molBlocks.length; i++) {
    molBlocks[i] = molBlocks[i].replaceAll('(-\nM  V30 ', '(')
      .replaceAll('-\nM  V30 ', '').replaceAll(' )', ')');
    const numbers = extractAtomsBondsNumbersV3000(molBlocks[i]);
    const coordinates = extractAtomDataV3000(molBlocks[i]);

    let indexAtoms = molBlocks[i].indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
    indexAtoms = molBlocks[i].indexOf('\n', indexAtoms);
    let index = indexAtoms;
    let indexEnd = indexAtoms;

    for (let j = 0; j < numbers.natom; j++) {
      if (coordinates.atomIndex[j] != 1 || i == 0) {
        //rewrite atom number
        index = molBlocks[i].indexOf('V30', index) + 4;
        indexEnd = molBlocks[i].indexOf(' ', index);
        const atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
        molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

        //rewrite coordinates
        index = molBlocks[i].indexOf(' ', index) + 1;
        index = molBlocks[i].indexOf(' ', index) + 1;
        indexEnd = molBlocks[i].indexOf(' ', index);

        const totalShift = xShift - coordinates.x[0];
        let coordinate = 
          Math.round(10000*(parseFloat(molBlocks[i].substring(index, indexEnd)) + totalShift))/10000;
        molBlocks[i] = molBlocks[i].slice(0, index) + coordinate + molBlocks[i].slice(indexEnd);

        index = molBlocks[i].indexOf(' ', index) + 1;
        indexEnd = molBlocks[i].indexOf(' ', index);
        coordinate = 
          Math.round(10000*(parseFloat(molBlocks[i].substring(index, indexEnd))))/10000;
        molBlocks[i] = molBlocks[i].slice(0, index) + coordinate + molBlocks[i].slice(indexEnd);

        index = molBlocks[i].indexOf('\n', index) + 1;
      } else {
        index = molBlocks[i].indexOf('M  V30', index) - 1;
        indexEnd = molBlocks[i].indexOf('\n', index + 1);
        molBlocks[i] = molBlocks[i].slice(0, index) + molBlocks[i].slice(indexEnd);
      }
    }

    const indexAtomsEnd = molBlocks[i].indexOf('M  V30 END ATOM');
    atomBlock += molBlocks[i].substring(indexAtoms + 1, indexAtomsEnd);

    let indexBonds = molBlocks[i].indexOf('M  V30 BEGIN BOND'); // V3000 index for bonds
    indexBonds = molBlocks[i].indexOf('\n', indexBonds);
    index = indexBonds;
    indexEnd = indexBonds;

    for (let j = 0; j < numbers.nbond; j++) {
      //rewrite bond number
      index = molBlocks[i].indexOf('V30', index) + 4;
      indexEnd = molBlocks[i].indexOf(' ', index);
      const bondNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + nbond;
      molBlocks[i] = molBlocks[i].slice(0, index) + bondNumber + molBlocks[i].slice(indexEnd);

      //rewrite atom pair in bond
      index = molBlocks[i].indexOf(' ', index) + 1;
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = molBlocks[i].indexOf(' ', index);
      let atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);
      index = molBlocks[i].indexOf(' ', index) + 1;
      indexEnd = Math.min(molBlocks[i].indexOf('\n', index), molBlocks[i].indexOf(' ', index));
      atomNumber = parseInt(molBlocks[i].substring(index, indexEnd)) + natom;
      molBlocks[i] = molBlocks[i].slice(0, index) + atomNumber + molBlocks[i].slice(indexEnd);

      index = molBlocks[i].indexOf('\n', index) + 1;
    }

    const indexBondEnd = molBlocks[i].indexOf('M  V30 END BOND');
    bondBlock += molBlocks[i].substring(indexBonds + 1, indexBondEnd);

    let indexCollection = molBlocks[i].indexOf('M  V30 MDLV30/STEABS ATOMS=('); // V3000 index for collections

    while (indexCollection != -1) {
      indexCollection += 28;
      const collectionEnd = molBlocks[i].indexOf(')', indexCollection);
      const collectionEntries = molBlocks[i].substring(indexCollection, collectionEnd).split(' ').slice(1);
      collectionEntries.forEach((e) => {
        collection.push(parseInt(e) + natom);
      });
      indexCollection = collectionEnd;
      indexCollection = molBlocks[i].indexOf('M  V30 MDLV30/STEABS ATOMS=(', indexCollection);
    }

    natom += numbers.natom - 1;
    nbond += numbers.nbond;
    xShift += coordinates.x[numbers.natom - 1] - coordinates.x[0];
  }

  const entries = 4;
  const collNumber = Math.ceil(collection.length / entries);

  //if (oclRender) {
    // collectionBlock += 'M  V30 MDLV30/STEABS ATOMS=(' + collection.length;

    // for (let j = 0; j < collection.length; j++)
    //   collectionBlock += ' ' + collection[j];

    // collectionBlock += ')\n';
  //} else {
    collectionBlock += 'M  V30 MDLV30/STEABS ATOMS=(' + collection.length + ' -\n';
    for (let i = 0; i < collNumber; i++) {
      collectionBlock += 'M  V30 ';
      const entriesCurrent = i + 1 == collNumber ? collection.length - (collNumber - 1)*entries : entries;
      for (let j = 0; j < entriesCurrent; j++) {
        collectionBlock += (j + 1 == entriesCurrent) ?
          (i == collNumber - 1 ? collection[entries*i + j] + ')\n' : collection[entries*i + j] + ' -\n') :
          collection[entries*i + j] + ' ';
      }
    }
  //}

  //generate file
  natom++;
  macroMolBlock += 'M  V30 COUNTS ' + natom + ' ' + nbond + ' 0 0 0\n';
  macroMolBlock += 'M  V30 BEGIN ATOM\n';
  macroMolBlock += atomBlock;
  macroMolBlock += 'M  V30 END ATOM\n';
  macroMolBlock += 'M  V30 BEGIN BOND\n';
  macroMolBlock += bondBlock;
  macroMolBlock += 'M  V30 END BOND\n';
  if(useChirality){
    macroMolBlock += 'M  V30 BEGIN COLLECTION\n';
    macroMolBlock += collectionBlock;
    macroMolBlock += 'M  V30 END COLLECTION\n';
  } else
    macroMolBlock = macroMolBlock.replace(/ CFG=\d/g, ' ');

  macroMolBlock += 'M  V30 END CTAB\n';
  macroMolBlock += 'M  END';

  return macroMolBlock;
}

function rotateNucleotidesV3000(molecule: string) {
  // @ts-ignore
  let molBlock = molecule.includes('M  END') ? molecule : OCL.Molecule.fromSmiles(molecule).toMolfileV3();
  const coordinates = extractAtomDataV3000(molBlock);
  const natom = coordinates.atomIndex.length;

  const indexFivePrime = coordinates.atomIndex.indexOf(1);
  const indexThreePrime = coordinates.atomIndex.indexOf(natom);

  //fix 5 prime if inadequate
  if (natom > 8)
    fix5Prime(coordinates, indexFivePrime, indexThreePrime);

  const xCenter = (coordinates.x[indexThreePrime] + coordinates.x[indexFivePrime])/2;
  const yCenter = (coordinates.y[indexThreePrime] + coordinates.y[indexFivePrime])/2;

  //place to center
  for (let i = 0; i < natom; i++) {
    coordinates.x[i] -= xCenter;
    coordinates.y[i] -= yCenter;
  }

  let angle = 0;
  if (coordinates.x[indexFivePrime] == 0)
    angle = coordinates.y[indexFivePrime] > coordinates.y[indexThreePrime] ? Math.PI/2 : 3*Math.PI/2;
  else if (coordinates.y[indexFivePrime] == 0)
    angle = coordinates.x[indexFivePrime] > coordinates.x[indexThreePrime] ? Math.PI : 0;
  else {
    const derivative = coordinates.y[indexFivePrime]/coordinates.x[indexFivePrime];
    angle = derivative > 0 ? Math.PI - Math.atan(derivative) : Math.atan(derivative);
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < natom; i++) {
    const xAdd = coordinates.x[i];
    coordinates.x[i] = xAdd*cos - coordinates.y[i]*sin;
    coordinates.y[i] = xAdd*sin + coordinates.y[i]*cos;
  }

  //place to right
  const xShift = coordinates.x[indexFivePrime];
  for (let i = 0; i < natom; i++)
    coordinates.x[i] -= xShift;

  //rewrite molBlock
  let index = molBlock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = molBlock.indexOf('\n', index);
  let indexEnd = index;
  for (let i = 0; i < natom; i++) {
    index = molBlock.indexOf('V30', index) + 4;
    index = molBlock.indexOf(' ', index) + 1;
    index = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', indexEnd);

    molBlock = molBlock.slice(0, index) +
      coordinates.x[i] + ' ' + coordinates.y[i] +
      molBlock.slice(indexEnd);

    index = molBlock.indexOf('\n', index) + 1;
  }

  return molBlock;
}

function invertNucleotidesV3000(molecule: string) {
  // @ts-ignore
  let molBlock = molecule.includes('M  END') ? molecule : OCL.Molecule.fromSmiles(molecule).toMolfileV3();
  const coordinates = extractAtomDataV3000(molBlock);
  const natom = coordinates.atomIndex.length;

  const xCenter = (Math.max(...coordinates.x) + Math.min(...coordinates.x))/2;
  const yCenter = (Math.max(...coordinates.y) + Math.min(...coordinates.y))/2;

  //place to center
  for (let i = 0; i < natom; i++) {
    coordinates.x[i] -= xCenter;
    coordinates.y[i] -= yCenter;
  }

  const angle = Math.PI;

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < natom; i++) {
    const xAdd = coordinates.x[i];
    coordinates.x[i] = xAdd*cos - coordinates.y[i]*sin;
    coordinates.y[i] = xAdd*sin + coordinates.y[i]*cos;
  }

  //place back
  const yShift = Math.max(...coordinates.y);
  for (let i = 0; i < natom; i++) {
    coordinates.x[i] += xCenter;
    coordinates.y[i] -= yShift;
  }

  //rewrite molBlock
  let index = molBlock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = molBlock.indexOf('\n', index);
  let indexEnd = index;
  for (let i = 0; i < natom; i++) {
    index = molBlock.indexOf('V30', index) + 4;
    index = molBlock.indexOf(' ', index) + 1;
    index = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', index) + 1;
    indexEnd = molBlock.indexOf(' ', indexEnd);

    molBlock = molBlock.slice(0, index) +
      coordinates.x[i] + ' ' + coordinates.y[i] +
      molBlock.slice(indexEnd);

    index = molBlock.indexOf('\n', index) + 1;
  }

  return molBlock;
}

function fix5Prime(coordinates: {atomIndex: number[], atomType: string[], x: number[], y: number[]},
  indexFivePrime: number, indexThreePrime: number) {
  const indexFivePrimeNeighbour = indexFivePrime + 1;
  const xShift = coordinates.x[indexFivePrimeNeighbour];
  const yShift = coordinates.y[indexFivePrimeNeighbour];
  const base3PrimeX = coordinates.x[indexThreePrime] - xShift;
  const base3PrimeY = coordinates.y[indexThreePrime] - yShift;
  const base5PrimeX = coordinates.x[indexFivePrime] - xShift;
  const base5PrimeY = coordinates.y[indexFivePrime] - yShift;

  const rotated5PrimeX = base5PrimeX*Math.cos(Math.PI*2/3) - base5PrimeY*Math.sin(Math.PI*2/3);
  const rotated5PrimeY = base5PrimeX*Math.sin(Math.PI*2/3) + base5PrimeY*Math.cos(Math.PI*2/3);

  const dx = base5PrimeX - base3PrimeX;
  const dy = base5PrimeY - base3PrimeY;
  const dxRotated = rotated5PrimeX - base3PrimeX;
  const dyRotated = rotated5PrimeY - base3PrimeY;

  if (Math.sqrt(dyRotated*dyRotated + dxRotated*dxRotated) >= Math.sqrt(dy*dy + dx*dx)) {
    coordinates.x[indexFivePrime] = rotated5PrimeX + xShift;
    coordinates.y[indexFivePrime] = rotated5PrimeY + yShift;
  }
}

function extractAtomsBondsNumbersV3000(molBlock: string): {natom: number, nbond: number} {
  molBlock = molBlock.replaceAll('\r', ''); //equalize old and new sdf standards
  let index = molBlock.indexOf('COUNTS') + 7; // V3000 index for atoms and bonds number
  let indexEnd = molBlock.indexOf(' ', index);

  const atomsNumber = parseInt(molBlock.substring(index, indexEnd));
  index = indexEnd + 1;
  indexEnd = molBlock.indexOf(' ', index);
  const bondsNumber = parseInt(molBlock.substring(index, indexEnd));

  return {natom: atomsNumber, nbond: bondsNumber};
}

function extractAtomDataV3000(molBlock: string) {
  const numbers = extractAtomsBondsNumbersV3000(molBlock);
  let index = molBlock.indexOf('M  V30 BEGIN ATOM'); // V3000 index for atoms coordinates
  index = molBlock.indexOf('\n', index);
  let indexEnd = index;

  const indexes: number[] = Array(numbers.natom);
  const types: string[] = Array(numbers.natom);
  const x: number[] = Array(numbers.natom);
  const y: number[] = Array(numbers.natom);

  for (let i = 0; i < numbers.natom; i++) {
    index = molBlock.indexOf('V30', index) + 4;
    indexEnd = molBlock.indexOf(' ', index);
    indexes[i] = parseInt(molBlock.substring(index, indexEnd));

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    types[i] = molBlock.substring(index, indexEnd);

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    x[i] = parseFloat(molBlock.substring(index, indexEnd));

    index = indexEnd + 1;
    indexEnd = molBlock.indexOf(' ', index);
    y[i] = parseFloat(molBlock.substring(index, indexEnd));

    index = molBlock.indexOf('\n', index) + 1;
  }

  return {atomIndex: indexes, atomType: types, x: x, y: y};
}
