%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt=(Maxstep=10) freq  int=(grid=99590) scrf=(cpcm, solvent=water) EmpiricalDispersion=GD3

solv sp

0 1
C                5.439445           -0.653483            0.021986
N                4.223876            0.076776           -0.257536
H                4.325994            0.879467           -0.869342
C                3.049091           -0.233640            0.293810
O                2.888076           -1.188052            1.042741
C                1.898663            0.717053           -0.029565
N                0.673376            0.014235            0.276801
H                0.767270           -0.781595            0.898637
C                6.602892            0.274382           -0.322209
O                6.442215            1.245789           -1.050895
N                7.786334           -0.073534            0.187231
H                7.888819           -0.890736            0.779299
C                9.010445            0.631486           -0.119571
C               10.161631           -0.338492            0.150464
O                9.998687           -1.317424            0.874032
N               11.329366           -0.012262           -0.415259
H               11.357858            0.788150           -1.029997
H                5.460996           -0.891295            1.091760
H                8.992426            0.901047           -1.181856
C               -1.665715           -0.586554            0.081547
N               -2.880878            0.131108           -0.230347
H               -2.775284            0.920201           -0.858978
C               -4.061954           -0.175142            0.310770
O               -4.227223           -1.114835            1.077823
C               -5.213108            0.761289           -0.049188
N               -6.439185            0.059993            0.255894
H               -6.350779           -0.722287            0.895223
C               -7.608730            0.362787           -0.315206
O               -7.751384            1.292545           -1.100382
C               -0.501929            0.339168           -0.265489
O               -0.653562            1.287066           -1.024943
H               -1.658673           -0.802586            1.156074
H               -5.195794            0.971800           -1.124680
C               12.551749           -0.762721           -0.221935
H               12.384308           -1.500641            0.562663
H               12.841931           -1.280742           -1.141146
H               13.360665           -0.092978            0.080567
C               -8.768628           -0.565973            0.033204
H               -8.758534           -0.783342            1.107439
H               -9.908761            0.906779           -0.925892
C              -11.175736           -0.180208            0.286591
C              -12.360067            0.646636           -0.145098
H              -12.094798            1.477380           -0.802470
H              -13.069637           -0.003893           -0.664443
H              -12.859242            1.036768            0.745415
O              -11.291157           -1.101985            1.106107
N               -9.991747            0.136544           -0.274245
C               -8.650152           -1.880551           -0.745809
H               -8.627526           -1.684845           -1.822133
H               -7.741408           -2.421591           -0.467136
H               -9.513456           -2.510681           -0.519172
C               -5.098721            2.078031            0.725191
H               -5.087251            1.886849            1.802368
H               -4.185124            2.613706            0.452888
H               -5.956842            2.711036            0.487463
C               -1.560585           -1.900385           -0.699078
H               -1.536139           -1.703465           -1.774979
H               -0.656492           -2.448886           -0.420739
H               -2.428689           -2.523871           -0.472930
C                2.029760            2.015648            0.772176
H                2.056364            1.800595            1.844538
H                2.941470            2.553541            0.497986
H                1.171140            2.657362            0.561112
C                5.525716           -1.951912           -0.786262
H                5.536398           -1.733750           -1.858276
H                6.430086           -2.511420           -0.531608
H                4.657206           -2.574934           -0.560355
C                9.161967            1.902572            0.721575
H                8.300028            2.551936            0.551084
H                9.211878            1.652110            1.785390
H               10.067987            2.448959            0.445408
H                1.902163            0.951032           -1.100236

