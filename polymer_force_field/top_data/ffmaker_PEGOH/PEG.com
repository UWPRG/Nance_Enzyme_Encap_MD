--Link1--
%chk=molecule
#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) NoSymm iop(6/42=6) opt

remark line goes here

0   1
    H   -6.1020000000       -0.3640000000        0.0250000000     
    O   -5.3540000000        0.2570000000        0.0330000000     
    C   -4.1660000000       -0.5230000000       -0.0090000000     
    H   -4.1450000000       -1.1950000000        0.8550000000     
    H   -4.1710000000       -1.1390000000       -0.9140000000     
    C   -2.9650000000        0.4030000000        0.0030000000     
    H   -2.9890000000        1.0320000000        0.9010000000     
    H   -2.9980000000        1.0680000000       -0.8690000000     
    O   -1.7780000000       -0.3760000000       -0.0190000000     
    C   -0.6150000000        0.4440000000       -0.0080000000     
    H   -0.6070000000        1.0960000000       -0.8900000000     
    H   -0.6120000000        1.0860000000        0.8810000000     
    C    0.6160000000       -0.4420000000       -0.0080000000     
    H    0.6060000000       -1.0950000000       -0.8900000000     
    H    0.6110000000       -1.0850000000        0.8810000000     
    O    1.7760000000        0.3770000000       -0.0190000000     
    C    2.9650000000       -0.4040000000        0.0030000000     
    H    2.9850000000       -1.0320000000        0.9020000000     
    H    2.9980000000       -1.0690000000       -0.8680000000     
    C    4.1660000000        0.5220000000       -0.0090000000     
    H    4.1520000000        1.2060000000        0.8460000000     
    H    4.1770000000        1.1340000000       -0.9170000000     
    O    5.3540000000       -0.2580000000        0.0330000000     
    H    6.1020000000        0.7650000000        0.0250000000     

