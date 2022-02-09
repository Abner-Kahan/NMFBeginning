%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt int=(grid=99590) freq EmpiricalDispersion=GD3

solv sp

0 1
C	5.5747440	4.4546890	-0.6325680
N	4.2983350	5.0167640	-0.2093540
H	4.2161200	6.0459300	-0.1552700
C	3.2279840	4.2538440	0.0719710
O	3.1963370	3.0167060	-0.0358670
C	2.0088070	5.0153050	0.6095030
N	0.8069630	4.5287510	-0.0603900
H	0.6885800	3.5148380	-0.1666900
C	6.6850240	4.9837720	0.2872910
O	6.7204580	6.1935610	0.5852810
N	7.5871120	4.0703770	0.6856500
H	7.4460010	3.0864760	0.4263300
C	8.7411570	4.3920350	1.5193920
C	10.0143600	3.8718450	0.8302570
O	10.1094970	2.6834010	0.4820590
N	10.9961750	4.7760510	0.6631880
H	10.8558060	5.7489690	0.9692170
H	5.5112620	3.3689700	-0.5459850
H	8.7847320	5.4785710	1.6144810
C	1.7115220	-4.6381810	0.4197330
C	2.8871240	-5.3744870	-0.2396310
O	2.9511800	-6.6211310	-0.2034260
N	3.8017980	-4.5738540	-0.8112560
H	3.6741960	-3.5595870	-0.7408550
C	5.0406560	-5.0438280	-1.4217390
C	6.2229590	-4.4812940	-0.6122420
O	6.2354170	-3.2844000	-0.2676930
N	7.1976220	-5.3562500	-0.3390940
H	7.1056440	-6.3321080	-0.6816130
N	0.4677040	-5.0576690	-0.2218800
H	0.3451670	-6.0605250	-0.4591920
C	-0.5378180	-4.2010830	-0.4577060
O	-0.5044530	-2.9859160	-0.1936100
C	-1.7810360	-4.8192090	-1.1235850
N	-2.9724820	-4.3186150	-0.4471360
H	-3.1213030	-3.3045790	-0.4503610
H	1.8174040	-3.5634910	0.2635410
H	5.0445090	-6.1343970	-1.3797020
H	-1.7649270	-5.9050720	-1.0125530
C	-1.3564520	4.6692830	-1.1711460
N	-2.6111220	5.2679430	-0.7380080
H	-2.6921050	6.2975280	-0.7460900
C	-3.6677220	4.5411500	-0.3355590
O	-3.7223880	3.3002460	-0.3452300
C	-4.8275730	5.3736750	0.2204440
N	-6.0934190	4.7493850	-0.1461840
H	-6.1464570	3.7257570	-0.1993860
C	-7.2142040	5.4697560	-0.2989510
O	-7.2840970	6.7097980	-0.1555260
C	-0.1956580	5.3458380	-0.4282440
O	-0.2031190	6.5800480	-0.2443960
H	-1.3810740	3.6075030	-0.9214310
H	-4.8007090	6.3664300	-0.2283710
C	-4.9753640	-4.3287720	0.9280490
C	-3.8279930	-5.0878090	0.2452200
O	-3.7224560	-6.3277730	0.3681860
N	-6.2538230	-4.7985670	0.3989710
H	-6.4125500	-5.8230410	0.3655550
C	-7.2744550	-3.9611260	0.1507850
O	-7.2022810	-2.7225390	0.2479090
C	-8.5863150	-4.6357890	-0.2931720
N	-9.6999380	-4.0871410	0.4695120
H	-10.0186800	-3.1471830	0.2261170
C	-10.1007000	-4.6155120	1.6494440
O	-9.6181410	-5.6507310	2.1371320
H	-4.8964770	-3.2663060	0.6934550
H	-8.5312960	-5.7067960	-0.0805080
C	-11.2291280	-3.8801130	2.3498270
H	-11.4848390	-2.9385300	1.8572580
H	-12.1076640	-4.5338350	2.3741780
H	-10.9308140	-3.6901690	3.3853180
C	12.2639190	4.4093970	0.0499330
H	12.7952100	3.6615000	0.6507750
H	12.1074050	3.9864030	-0.9487380
H	12.8820370	5.3066930	-0.0330900
C	8.3788840	-4.9873730	0.4352020
H	8.5431980	-3.9152410	0.3146060
C	-8.4555230	4.6787170	-0.7215430
H	-8.4452530	3.6993160	-0.2409510
H	-9.4939120	6.4115460	-0.1286700
C	-10.7961390	4.8130020	0.0415960
C	-11.9087240	5.7183500	0.5386930
H	-11.6212620	6.7727200	0.5867560
H	-12.7723870	5.6106320	-0.1249060
H	-12.2134550	5.3831880	1.5351150
O	-10.9670850	3.5897640	-0.0720830
N	-9.6255280	5.4153020	-0.2647970
C	9.5874460	-5.7667570	-0.1029840
H	10.6050690	-3.9938540	-0.2126500
C	11.9396330	-5.5682830	-0.7919290
H	12.7598150	-5.2869650	-0.1205760
H	11.8446690	-6.6557530	-0.8092610
H	12.1803930	-5.2172040	-1.8033420
O	9.5417620	-6.9955690	-0.2665280
N	10.6825760	-5.0060010	-0.3293390
C	8.1880470	-5.3073460	1.9298700
H	7.3192390	-4.7645710	2.3159810
H	8.0314280	-6.3809230	2.0751370
H	9.0712520	-5.0030760	2.5028490
C	-8.4577810	4.4739780	-2.2499710
H	-8.4017600	5.4347470	-2.7732850
H	-7.6079070	3.8533500	-2.5490240
H	-9.3800580	3.9648590	-2.5459760
C	-4.6762320	5.5109150	1.7488400
H	-4.6773560	4.5246300	2.2246890
H	-3.7382250	6.0185880	1.9956500
H	-5.5044800	6.1014440	2.1515080
C	5.1446410	-4.5862020	-2.8854290
H	4.2966010	-4.9715420	-3.4603670
H	5.1458390	-3.4926130	-2.9489640
H	6.0699620	-4.9626990	-3.3332380
C	1.6970700	-4.9415340	1.9280050
H	0.8603240	-4.4175880	2.4007550
H	1.5906940	-6.0160850	2.1046630
H	2.6282820	-4.6036420	2.3968450
C	-1.1478860	4.8366680	-2.6886610
H	-1.1361880	5.8972860	-2.9609770
H	-0.2011880	4.3819160	-2.9989190
H	-1.9626390	4.3442730	-3.2290510
C	-1.8021090	-4.4562570	-2.6175890
H	-2.6896000	-4.8833260	-3.0950880
H	-1.8194660	-3.3692350	-2.7496210
H	-0.9117610	-4.8528210	-3.1162390
C	1.9223590	4.8265570	2.1358980
H	1.8201440	3.7651090	2.3840760
H	2.8246250	5.2124700	2.6226500
H	1.0561030	5.3679390	2.5287240
C	-4.8935040	-4.5322420	2.4515050
H	-5.0090970	-5.5885600	2.7119280
H	-3.9291250	-4.1806920	2.8359010
H	-5.6908610	-3.9632930	2.9397690
C	5.8988220	4.8446620	-2.0864390
H	5.9360340	5.9339740	-2.1907220
H	6.8642140	4.4264470	-2.3911490
H	5.1251970	4.4538670	-2.7551160
C	-8.8054100	-4.4103670	-1.7949120
H	-8.8778240	-3.3402750	-2.0161280
H	-7.9719670	-4.8270920	-2.3694010
H	-9.7308190	-4.8978470	-2.1174530
C	8.6126870	3.7562830	2.9138480
H	7.7066410	4.1185460	3.4098890
H	8.5608740	2.6660050	2.8325680
H	9.4782010	4.0136400	3.5346130
H	2.1092140	6.0773070	0.3843180
O	1.4334530	8.7090080	0.0739940
H	1.2066510	9.2537310	-0.6909330
H	0.8391310	7.9223920	0.0184630
O	3.9971380	7.8268540	-0.0916420
H	4.5507870	8.1467600	0.6319460
H	3.1003240	8.2290530	0.0170500
O	-5.5162620	8.7923840	0.1768770
H	-5.9153190	9.5273840	-0.3056870
H	-6.1558630	8.0458160	0.1012390
O	-2.9965350	8.0933920	-0.6681050
H	-2.3026770	8.3690350	-0.0549320
H	-3.8546590	8.4330600	-0.3204360
O	0.2511830	-7.7621940	-0.8315850
H	-0.4780090	-8.1608510	-0.2999800
H	1.0969280	-7.9446010	-0.3951410
O	-2.0265740	-8.4952490	0.4453350
H	-2.6002100	-7.6970810	0.4313520
H	-2.5420430	-9.1911500	0.0177960
O	7.0168270	-7.9622600	-1.2855680
H	7.8713700	-8.1309220	-0.8533210
H	6.2991970	-8.4189710	-0.7936290
O	4.6349210	-8.8192770	-0.2182180
H	4.0748390	-8.0123350	-0.1900540
H	4.2065890	-9.4026420	-0.8579270
O	8.1839970	8.3877580	0.9794820
H	8.0603670	8.8834960	0.1598720
H	7.6312570	7.5743170	0.8830350
O	10.6492170	7.4657170	1.4696370
H	10.8893640	7.7787090	2.3504110
H	9.7668940	7.8745250	1.2676590
O	-7.7178970	-7.4728170	2.7995510
H	-8.3916850	-6.7937470	2.5699780
H	-8.2111420	-8.2096770	3.1805910
O	-6.3400600	-7.5875640	0.4704720
H	-5.3838210	-7.4933550	0.6217950
H	-6.7665530	-7.7091000	1.3518060
C	-5.7423480	0.0995470	0.6239290
N	-4.5156400	0.5553180	-0.0188210
H	-4.3590870	1.5620070	-0.1102930
C	-3.5629110	-0.2955720	-0.4481870
O	-3.6344930	-1.5305360	-0.3396370
C	-2.3544960	0.3588150	-1.1297390
N	-1.1379180	-0.2060070	-0.5561650
H	-1.0770750	-1.2237520	-0.4427090
C	-6.9447010	0.7652820	-0.0571660
O	-6.9387490	1.9855150	-0.3063060
N	-7.9911350	-0.0444090	-0.3122600
H	-7.8952650	-1.0519990	-0.1326450
C	-9.2181290	0.4468490	-0.9320420
C	-10.4327830	-0.1322900	-0.1922360
O	-10.5509440	-1.3544400	-0.0080520
N	-11.3579420	0.7693450	0.1890670
H	-11.1777150	1.7657710	0.0601680
H	-5.7999750	-0.9822330	0.5010740
H	-9.2129810	1.5344000	-0.8386660
C	1.0679940	-0.1944890	0.4708800
N	2.3401630	0.2890740	-0.0478390
H	2.5215440	1.2961880	-0.0287310
C	3.2998650	-0.5418090	-0.5005290
O	3.1923970	-1.7789260	-0.5224940
C	4.5616320	0.1452280	-1.0392540
N	5.7328050	-0.5050860	-0.4620260
H	5.7607120	-1.5313160	-0.4371880
C	6.7917810	0.1896590	-0.0049800
O	6.8729490	1.4310470	-0.0470050
C	-0.0830400	0.5572030	-0.2099080
O	-0.0370920	1.7894980	-0.3733310
H	0.9984810	-1.2587340	0.2424430
H	4.5745820	1.1907040	-0.7288500
C	-12.5913130	0.3727950	0.8506650
H	-12.9711100	-0.5520820	0.4076530
H	-12.4382310	0.1988760	1.9238450
H	-13.3311950	1.1677430	0.7251530
C	7.8950180	-0.6567250	0.6410420
H	7.8379120	-1.6780060	0.2636240
H	9.3422180	0.8922080	0.3147030
C	10.2556300	-0.9252900	0.0181570
C	11.5571210	-0.2228090	-0.3251510
H	11.4878830	0.8627280	-0.2257800
H	12.3443690	-0.6090520	0.3301520
H	11.8286850	-0.4838750	-1.3538730
O	10.1831050	-2.1671350	0.0391920
N	9.2003670	-0.1197650	0.2799500
C	7.6858110	-0.6851470	2.1688690
H	7.7164550	0.3271970	2.5868050
H	6.7190370	-1.1369600	2.4152530
H	8.4759520	-1.2825660	2.6343460
C	4.5699020	0.0676320	-2.5772440
H	4.5529850	-0.9753200	-2.9103450
H	3.6951310	0.5797610	-2.9933430
H	5.4731940	0.5479690	-2.9664290
C	0.9703600	0.0087490	1.9949160
H	1.0439380	1.0716460	2.2502660
H	0.0172120	-0.3758820	2.3732860
H	1.7838270	-0.5289320	2.4922450
C	-2.4315600	0.1275500	-2.6498140
H	-2.4446450	-0.9430690	-2.8791930
H	-3.3398030	0.5813560	-3.0618170
H	-1.5632800	0.5819050	-3.1374970
C	-5.7523470	0.4502930	2.1235630
H	-5.6869270	1.5340330	2.2694050
H	-6.6719630	0.0884670	2.5958720
H	-4.8990760	-0.0232000	2.6193970
C	-9.2910350	0.0660160	-2.4217610
H	-8.4274910	0.4783720	-2.9530450
H	-9.2922820	-1.0218780	-2.5408830
H	-10.2047950	0.4638920	-2.8773510
H	-2.3493410	1.4305180	-0.9271800

