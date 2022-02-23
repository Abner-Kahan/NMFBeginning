%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt freq  int=(grid=99590) scrf=(cpcm, solvent=water) EmpiricalDispersion=GD3

solv sp

0 1
C	-12.7984860	-1.0064860	2.0958270
C	-11.6121510	-1.1530850	1.1666880
N	-11.8816880	-1.5138140	-0.1255020
C	-10.8082700	-1.9365750	-1.0296710
C	-9.7513050	-0.8430420	-1.2429020
N	-10.1823050	0.4420590	-1.2840670
C	-9.2502750	1.5489940	-1.4881140
C	-8.1944700	1.6570020	-0.3693380
N	-8.5948640	1.3542250	0.8922940
C	-7.6444470	1.3826240	2.0034270
C	-6.5706160	0.2849720	1.8779850
N	-6.9967840	-0.9405750	1.4813640
C	-6.0462680	-2.0229390	1.2402250
C	-5.0388560	-1.6652060	0.1276750
N	-5.5398940	-1.0759800	-0.9846800
C	-4.6589860	-0.6758800	-2.0788820
C	-3.6207590	0.3732860	-1.6319430
N	-4.0694120	1.3952430	-0.8665250
C	-3.1504190	2.4179650	-0.3744110
C	-2.0549370	1.8225880	0.5326090
N	-2.4501220	0.9157370	1.4557790
C	-1.4766050	0.2760790	2.3355030
C	-0.4309430	-0.5382690	1.5476360
N	-0.8922580	-1.3152080	0.5406570
C	0.0279470	-2.0931400	-0.2820430
C	1.0253580	-1.1971520	-1.0437600
N	0.5257650	-0.0811090	-1.6217920
C	1.4005050	0.8479220	-2.3293010
C	2.4730840	1.4561100	-1.4029930
N	2.0643670	1.8683320	-0.1808640
C	3.0194590	2.4126870	0.7791690
C	4.1048370	1.3863590	1.1659900
N	3.6861580	0.1232660	1.4165520
C	4.6459230	-0.9233080	1.7488580
C	5.6639680	-1.1612580	0.6158190
N	5.1782610	-1.2023310	-0.6458370
C	6.0689910	-1.4104720	-1.7823470
C	7.0704080	-0.2543640	-1.9833070
N	6.6148390	0.9999720	-1.7564280
C	7.4963000	2.1512240	-1.9090800
C	8.6101220	2.2213910	-0.8415430
N	8.2962700	1.7346600	0.3915430
C	9.2389090	1.8133940	1.4970920
C	10.1876110	0.6018790	1.6319730
N	9.7395410	-0.5900200	1.1397700
C	10.4664730	-1.8274790	1.3963960
C	10.7369070	-2.6753300	0.1373600
N	10.3588800	-2.1503300	-1.0519530
O	-10.4513670	-0.9443860	1.5417910
C	-11.3948250	-2.3846610	-2.3742460
O	-8.5675320	-1.1532270	-1.4014820
C	-10.0138490	2.8700280	-1.6304060
O	-7.0603370	2.0611740	-0.6403050
C	-8.3825070	1.2832090	3.3436830
O	-5.3887090	0.5271540	2.1390070
C	-6.7890090	-3.3237600	0.9168710
O	-3.8389210	-1.9224720	0.2683020
C	-5.4820980	-0.1573300	-3.2631720
O	-2.4412370	0.2693070	-1.9871570
C	-3.9251720	3.5188070	0.3588510
O	-0.8799960	2.1895410	0.4123400
C	-2.1903640	-0.6067730	3.3654600
O	0.7652200	-0.4798490	1.8574160
C	-0.7516870	-2.9850270	-1.2549070
O	2.2173580	-1.5193470	-1.1216140
C	0.5724110	1.9493460	-3.0009850
O	3.6415240	1.5648260	-1.7938920
C	2.2863190	2.9225270	2.0250420
O	5.2887210	1.7313830	1.2529870
C	3.9125380	-2.2234080	2.0987600
O	6.8601180	-1.3278780	0.8826030
C	5.2509030	-1.6320000	-3.0601800
O	8.2223290	-0.4939250	-2.3628550
C	6.6785820	3.4484430	-1.8853520
O	9.6891240	2.7443340	-1.1172380
C	8.4947400	2.0152710	2.8259290
O	11.2701610	0.7462750	2.2025130
O	11.2859910	-3.7780440	0.2446240
H	-12.7225470	-0.0527800	2.6255370
H	-12.8253090	-1.8021110	-0.3535480
H	-10.2516610	-2.7651900	-0.5770830
H	-11.1631400	0.6338920	-1.1237880
H	-8.6687910	1.3583670	-2.3953280
H	-9.4574390	0.8341580	1.0279990
H	-7.0847590	2.3199550	1.9511380
H	-7.9743770	-1.0816830	1.2513760
H	-5.4259870	-2.1569520	2.1317210
H	-6.5441280	-0.9519790	-1.0876760
H	-4.0548510	-1.5344330	-2.3903320
H	-5.0597310	1.4810930	-0.6491190
H	-2.5977560	2.8448610	-1.2178160
H	-3.4345770	0.6764890	1.5540230
H	-0.8875890	1.0458660	2.8450510
H	-1.8912140	-1.3728310	0.3525340
H	0.6614350	-2.7084000	0.3656080
H	-0.4720100	0.1182210	-1.5771370
H	1.9782390	0.3009680	-3.0819330
H	1.0814520	1.8085270	0.0791610
H	3.5796610	3.2289700	0.3110140
H	2.6946780	-0.1071780	1.3851810
H	5.2615270	-0.5987620	2.5946700
H	4.1768460	-1.1167510	-0.8109640
H	6.7050950	-2.2801830	-1.5868250
H	5.6338620	1.1555120	-1.5336750
H	8.0401190	2.0650030	-2.8552150
H	7.3408470	1.4443650	0.5817580
H	9.9057810	2.6594310	1.3100830
H	8.7569310	-0.6712890	0.8840000
H	11.4513000	-1.5268310	1.7671180
C	9.7698780	-2.6856520	2.4668030
H	9.8720340	-1.2624000	-1.0902140
H	-7.6588760	1.2672750	4.1636720
H	-9.0368390	2.1504980	3.4812830
H	-8.9900480	0.3734430	3.3958710
H	-6.0680960	-4.1241840	0.7278790
H	-7.4184630	-3.6234580	1.7619950
H	-7.4184650	-3.2063930	0.0277270
H	-10.5878890	-2.6861140	-3.0471820
H	-12.0573740	-3.2468450	-2.2379100
H	-11.9610760	-1.5783950	-2.8538250
H	-4.8127830	0.1374510	-4.0764390
H	-6.1528930	-0.9396130	-3.6345680
H	-6.0813690	0.7125910	-2.9730420
H	-9.3055200	3.6899180	-1.7746010
H	-10.6806180	2.8392680	-2.4989840
H	-10.6070660	3.0867970	-0.7342190
H	-3.2292440	4.2817940	0.7192900
H	-4.6446430	3.9940950	-0.3160600
H	-4.4695750	3.1082940	1.2159310
H	-1.4517580	-1.0863430	4.0143520
H	-2.8604360	-0.0034250	3.9869930
H	-2.7813330	-1.3842420	2.8699450
H	5.9272190	-1.7828210	-3.9066560
H	4.6120660	-2.5157640	-2.9578930
H	4.6133340	-0.7657490	-3.2676590
H	1.2370960	2.6405880	-3.5274980
H	-0.1245480	1.5141690	-3.7250770
H	-0.0040450	2.5107810	-2.2582190
H	4.6406040	-3.0027090	2.3427540
H	3.2585030	-2.0715420	2.9640700
H	3.3008140	-2.5620470	1.2559200
H	7.3508580	4.3036470	-1.9982860
H	5.9488180	3.4585730	-2.7021170
H	6.1371970	3.5463600	-0.9377870
H	3.0106110	3.3261250	2.7387250
H	1.5810960	3.7159680	1.7553610
H	1.7299080	2.1119820	2.5075280
H	9.2163440	2.0551030	3.6467810
H	7.9227470	2.9487540	2.8037300
H	7.7935590	1.1932560	3.0117640
H	-0.0519520	-3.5625320	-1.8660400
H	-1.3919350	-3.6821370	-0.7039720
H	-1.3825220	-2.3814850	-1.9160580
H	10.3139710	-3.6248740	2.5944400
H	9.7467420	-2.1453350	3.4193940
H	8.7392580	-2.9083730	2.1687550
H	-12.7495020	-1.8020450	2.8471950
H	-13.7627920	-1.0629210	1.5824330
C	10.5562850	-2.8661530	-2.2989980
H	10.3391920	-3.9300530	-2.1582130
H	9.8818420	-2.4397580	-3.0453930
H	11.5921540	-2.7814570	-2.6544020

