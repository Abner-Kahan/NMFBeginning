%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt int=(grid=99590) freq EmpiricalDispersion=GD3

solv sp

0 1
C	5.8296540	-2.0309880	0.7106340
N	4.5721920	-2.6818400	0.3664090
H	4.6440000	-3.6645580	0.1201060
C	3.4104480	-2.0084020	0.2319320
O	3.3003610	-0.7952200	0.4613340
C	2.2219850	-2.8331600	-0.2761640
N	0.9997870	-2.2008040	0.1931210
H	1.0098040	-1.1917370	0.3515960
C	6.9472020	-2.7505790	-0.0562750
O	6.8577330	-3.9535330	-0.3303860
N	8.0118730	-1.9679380	-0.3507920
H	7.9547020	-0.9596050	-0.1815990
C	9.1853170	-2.4829160	-1.0398890
C	10.4321870	-1.9151390	-0.3462590
O	10.5923250	-0.6979610	-0.2165630
N	11.3323900	-2.8387860	0.0834630
H	11.0959450	-3.8151560	-0.0275470
H	5.7778970	-0.9878400	0.3967900
H	9.1477590	-3.5731010	-0.9516190
C	-2.2201830	1.9023910	0.9245580
C	-3.3698560	2.7375350	0.3440070
O	-3.2297140	3.9471070	0.1240660
N	-4.5197880	2.0560380	0.1438860
H	-4.5422580	1.0414350	0.2579500
C	-5.7220890	2.7141220	-0.3419210
C	-6.9262740	1.8505790	0.0540830
O	-6.8101120	0.6267750	0.2240850
N	-8.1008150	2.5010660	0.1609590
H	-8.1809030	3.4934740	-0.0448080
N	-0.9785660	2.5767470	0.5704440
H	-1.0696100	3.5719160	0.3836660
C	0.1805810	1.9218880	0.3621060
O	0.3059640	0.6979660	0.5240410
C	1.3447730	2.7822340	-0.1434540
N	2.5884290	2.1281460	0.2301710
H	2.5811810	1.1142740	0.3493130
H	-2.2258460	0.9082140	0.4748070
H	-5.7892480	3.6960180	0.1394690
H	1.3210200	3.7657930	0.3386890
C	-1.3358260	-2.0952850	0.8705710
N	-2.5560490	-2.7446700	0.4109730
H	-2.4553560	-3.7278740	0.1737120
C	-3.7016100	-2.0743450	0.1739070
O	-3.8353450	-0.8620030	0.3982440
C	-4.8342570	-2.9002580	-0.4473640
N	-6.0978980	-2.2587930	-0.1216400
H	-6.0987970	-1.2518810	0.0506510
C	-7.2624240	-2.9447900	-0.1539470
O	-7.3529050	-4.1425000	-0.4531360
C	-0.1561880	-2.8938510	0.2996060
O	-0.2790080	-4.0900240	0.0073170
H	-1.3090620	-1.0744150	0.4857570
H	-4.8442220	-3.9070010	-0.0150560
C	4.9603280	1.9838880	0.7581460
C	3.7515670	2.8150950	0.3071740
O	3.8594440	4.0247950	0.0737050
N	6.1574100	2.6664340	0.2872850
H	6.0528000	3.6639720	0.1237430
C	7.2909780	2.0146890	-0.0399140
O	7.4191330	0.7870660	0.1026250
C	8.4001940	2.8695720	-0.6543100
N	9.6896300	2.2795590	-0.3425930
H	9.7638660	1.2667650	-0.2712880
C	10.8059900	3.0641990	-0.2679910
O	10.7671790	4.2869180	-0.4249730
H	4.9141730	0.9955780	0.2975230
H	8.3879640	3.8728890	-0.2157090
C	12.1019800	2.3307860	0.0398870
H	11.9755410	1.2464880	0.0848780
H	12.8367040	2.5888710	-0.7292250
H	12.4875330	2.6994420	0.9964500
C	12.5798490	-2.4921180	0.7493030
H	12.6429400	-1.4037020	0.7958510
H	12.6046510	-2.8937380	1.7690080
H	13.4414230	-2.8775620	0.1920450
C	-9.3684570	1.8234740	0.4005370
H	-9.3041660	0.8088120	0.0050940
C	-8.5019540	-2.1362590	0.2550880
H	-8.4373030	-1.1322690	-0.1684100
H	-9.5050590	-3.7953250	-0.5400540
C	-10.8155630	-2.1921530	-0.6099590
C	-11.9026180	-3.0512280	-1.2313130
H	-11.6124440	-4.0993450	-1.3481310
H	-12.7976930	-2.9929670	-0.6041330
H	-12.1585210	-2.6373850	-2.2118980
O	-10.9965180	-0.9812950	-0.3987030
N	-9.6585110	-2.8185970	-0.3071900
C	-10.4521390	2.6110370	-0.3520580
H	-11.4755340	0.8596630	-0.5960930
C	-12.6756320	2.4506500	-1.3895550
H	-13.6016690	2.1880560	-0.8647730
H	-12.5603900	3.5362730	-1.3948280
H	-12.7483540	2.0986540	-2.4263750
O	-10.3213090	3.8248490	-0.5501040
N	-11.5246540	1.8717180	-0.7171560
C	-9.7049800	1.7546800	1.9023450
H	-8.9002720	1.2426740	2.4387900
H	-9.8224630	2.7611400	2.3181760
H	-10.6354850	1.1994940	2.0580010
C	-8.5874980	-2.0197780	1.7894680
H	-8.5839430	-3.0119730	2.2538550
H	-7.7403570	-1.4436090	2.1733810
H	-9.5120120	-1.5040900	2.0652190
C	-4.6210300	-3.0240190	-1.9719680
H	-4.5846870	-2.0323590	-2.4354270
H	-3.6869610	-3.5499450	-2.1977670
H	-5.4510520	-3.5895360	-2.4037350
C	-5.6975660	2.9179810	-1.8722610
H	-4.8165200	3.5072190	-2.1403920
H	-5.6519270	1.9518500	-2.3861050
H	-6.5906050	3.4534430	-2.2126540
C	-2.3810070	1.7558850	2.4500760
H	-1.5244650	1.2109520	2.8578920
H	-2.4351300	2.7388600	2.9304600
H	-3.2923520	1.1967960	2.6834420
C	-1.2536990	-2.0470530	2.4088270
H	-1.2237170	-3.0590790	2.8270150
H	-0.3559640	-1.5067280	2.7249570
H	-2.1302580	-1.5269670	2.8062960
C	1.2283860	2.9822020	-1.6703480
H	2.0792600	3.5737410	-2.0185470
H	1.2308090	2.0150000	-2.1842010
H	0.3060810	3.5137220	-1.9289670
C	2.2716340	-2.9363040	-1.8160350
H	2.2501290	-1.9385910	-2.2671700
H	3.1803040	-3.4508260	-2.1475950
H	1.4048010	-3.5043690	-2.1642210
C	4.9577600	1.8151690	2.2900280
H	4.9430890	2.7912040	2.7870640
H	4.0820310	1.2411320	2.6080200
H	5.8585580	1.2768570	2.5999590
C	6.1076070	-2.0808670	2.2252550
H	6.2004670	-3.1167820	2.5689610
H	7.0341980	-1.5468200	2.4592630
H	5.2853570	-1.6020020	2.7656950
C	8.1712810	3.0050720	-2.1766010
H	8.1729230	2.0198070	-2.6559260
H	7.2167340	3.4979210	-2.3932860
H	8.9782020	3.6084030	-2.6014330
C	9.1966730	-2.0904780	-2.5277940
H	8.3036710	-2.4880130	-3.0195160
H	9.2072930	-1.0010740	-2.6311840
H	10.0819160	-2.4939910	-3.0328950
H	2.2599010	-3.8449940	0.1422980

