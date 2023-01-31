%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt freq  int=(grid=99590) scrf=(cpcm, solvent=water) EmpiricalDispersion=GD3

solv sp

0 1
N	-4.4195420	-1.3810110	0.1554850
H	-4.2035200	-0.4558820	-0.2117680
C	-2.2193030	-2.0173030	-0.6941040
C	-3.4826040	-2.3627490	0.1288520
O	-3.6140920	-3.4723080	0.6468190
N	-1.0179860	-2.4595260	0.0028800
H	-0.5530980	-3.3051950	-0.3124840
C	-0.4334940	-1.7702820	0.9997670
O	-0.8528210	-0.6740690	1.4099530
C	0.8300670	-2.4418710	1.5847050
N	2.0087490	-2.0320620	0.8187500
H	2.4746680	-1.1705540	1.0981250
C	2.3223500	-2.6503820	-0.3494090
O	1.6763520	-3.6030150	-0.7969740
H	-2.1580990	-0.9342290	-0.8322690
H	0.7282920	-3.5232300	1.4498240
C	3.6215300	-2.1896320	-1.0555640
C	-2.2996360	-2.7128040	-2.0586170
H	-2.3755790	-3.7975400	-1.9250160
H	-3.1860440	-2.3748410	-2.6039130
H	-1.4115970	-2.4881320	-2.6589860
C	1.0084730	-2.1154560	3.0646790
H	1.1883280	-1.0472460	3.2085510
H	0.1105400	-2.3935630	3.6244580
H	1.8620200	-2.6703350	3.4676090
C	-1.4668230	1.9667780	-0.6526580
N	-2.5412100	2.9419300	-0.8224430
H	-2.2676550	3.8994030	-1.0064390
C	-3.8467900	2.5748830	-0.8794360
O	-4.2159210	1.4015160	-0.7438360
C	-4.8479900	3.6937070	-1.1085620
C	-0.3375670	2.6409140	0.1302630
O	0.1961040	3.6810790	-0.2890330
N	0.0250020	2.0209030	1.2778670
H	-0.3358400	1.0766720	1.4615860
C	1.1958530	2.4615470	2.0381110
C	2.4183510	1.7523430	1.4246330
O	2.8999340	0.7114580	1.8917840
N	2.8836770	2.3267090	0.2835090
H	2.2925320	3.0358530	-0.1491760
H	-1.8734130	1.1356290	-0.0740200
H	1.2956360	3.5379410	1.8673020
C	3.8591260	1.6286320	-0.5470700
C	-0.9353470	1.4602440	-2.0045860
H	-0.5529860	2.2992170	-2.5950780
H	-0.1176880	0.7468310	-1.8646780
H	-1.7472540	0.9773060	-2.5583240
C	1.0349340	2.1688740	3.5257140
H	0.1656670	2.7068940	3.9164190
H	0.8949320	1.0995980	3.7044810
H	1.9256700	2.4881850	4.0757020
H	-4.3837000	4.6792490	-1.2103910
C	3.2730910	0.2798850	-1.0224180
O	2.2113170	0.2294010	-1.6388250
N	4.0211600	-0.8205110	-0.7193700
H	4.7884450	-0.6924430	-0.0740440
H	-5.4231370	3.4696870	-2.0122770
H	-5.5483410	3.7112850	-0.2677830
C	-5.6916390	-1.5551290	0.8352500
H	-5.8414510	-2.6208200	1.0205130
H	-5.7082340	-1.0279520	1.7983290
H	-6.5077270	-1.1712930	0.2131000
C	3.5511550	-2.4343350	-2.5682550
H	2.8130480	-1.7741760	-3.0278730
H	4.5310050	-2.2479110	-3.0202710
H	3.2577070	-3.4697850	-2.7533310
H	4.4040680	-2.8382630	-0.6403840
H	4.7515850	1.4454630	0.0642450
C	4.2215300	2.4915040	-1.7584820
H	4.6197120	3.4585270	-1.4337470
H	4.9790780	1.9938260	-2.3718690
H	3.3348850	2.6580900	-2.3788760
