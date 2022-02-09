%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt int=(grid=99590) freq EmpiricalDispersion=GD3

solv sp

0 1
C	5.6026840	-4.4732850	0.8427680
N	4.3529680	-5.1188900	0.4649630
H	4.3980260	-6.1235700	0.3265780
C	3.2352110	-4.4230160	0.1667030
O	3.1618550	-3.1895630	0.2656800
C	2.0559230	-5.2512190	-0.3531460
N	0.8270390	-4.5831150	0.0393160
H	0.8278700	-3.5627020	0.1054230
C	6.7492500	-5.2567720	0.1925970
O	6.6791940	-6.4811760	0.0362840
N	7.8163540	-4.4927810	-0.1445330
H	7.7512350	-3.4756170	-0.0441650
C	9.0316910	-5.0562350	-0.7088840
C	10.2243570	-4.5247100	0.0991460
O	10.4271900	-3.3107690	0.2161680
N	11.0282960	-5.4725300	0.6462500
H	10.7511560	-6.4394160	0.5446790
H	5.5879030	-3.4517640	0.4619650
H	8.9456390	-6.1427150	-0.6105690
C	-2.2656740	-0.3090410	0.6271810
C	-3.4595960	0.4368250	0.0173990
O	-3.4042200	1.6579140	-0.2001280
N	-4.5582650	-0.3166250	-0.1918920
H	-4.5070250	-1.3302470	-0.0706690
C	-5.8054410	0.2528260	-0.6931440
C	-6.9825120	-0.4962020	-0.0566080
O	-6.9523530	-1.7325600	0.0736950
N	-8.0366710	0.2694810	0.2856170
H	-7.9583570	1.2900870	0.2193870
N	-1.0271280	0.2577470	0.1036960
H	-0.9816450	1.2706260	-0.0331800
C	0.0763400	-0.4920230	-0.0867760
O	0.1335530	-1.7109970	0.1507830
C	1.2734890	0.2468360	-0.6971300
N	2.5083900	-0.3257820	-0.1701080
H	2.5401520	-1.3338570	-0.0052860
H	-2.3069340	-1.3614870	0.3429080
H	-5.8398860	1.3042750	-0.4063690
H	1.2353260	1.3000020	-0.4147170
C	-1.5100060	-4.4131900	0.6735630
N	-2.7340810	-5.0783120	0.2497850
H	-2.6577800	-6.0792330	0.0932300
C	-3.8594990	-4.3961660	-0.0417530
O	-3.9527580	-3.1674160	0.0962600
C	-5.0179060	-5.2226020	-0.6089040
N	-6.2627520	-4.5590820	-0.2568350
H	-6.2558260	-3.5441160	-0.1258700
C	-7.4333820	-5.2336090	-0.2232230
O	-7.5487630	-6.4385840	-0.4843740
C	-0.3315730	-5.2651970	0.1892080
O	-0.4545440	-6.4815040	-0.0008820
H	-1.4655640	-3.4279620	0.2067260
H	-5.0261340	-6.2204900	-0.1566290
C	4.8107930	-0.3244340	0.6276440
C	3.6207920	0.4172230	0.0059900
O	3.6879690	1.6296450	-0.2525350
N	6.0557410	0.2477400	0.1239650
H	6.1047910	1.2609920	-0.0116210
C	7.1573500	-0.5031460	-0.0631590
O	7.1966750	-1.7284920	0.1525960
C	8.3767060	0.2309270	-0.6311240
N	9.5932640	-0.3629950	-0.0879970
H	9.6332560	-1.3760910	0.0225630
C	10.6997860	0.3731870	0.1691180
O	10.7459860	1.6029640	-0.0100820
H	4.7762960	-1.3756410	0.3383860
H	8.3453770	1.2779280	-0.3272310
C	11.8969150	-0.3832360	0.7142090
H	11.7388670	-1.4631500	0.7351210
H	12.7687970	-0.1441860	0.0972150
H	12.1012950	-0.0193910	1.7271190
C	12.2072900	-5.1651110	1.4433870
H	12.3803470	-4.0892260	1.3874660
H	12.0598410	-5.4457210	2.4930850
H	13.0852120	-5.6915770	1.0530120
C	-9.2570620	-0.2922550	0.8561230
H	-9.2754330	-1.3552520	0.6094460
C	-8.6452600	-4.3999620	0.2171220
H	-8.5938960	-3.4082030	-0.2370540
H	-9.6982270	-6.0682530	-0.4931520
C	-11.0002550	-4.4547430	-0.5104320
C	-12.1282140	-5.3094870	-1.0604990
H	-11.8537510	-6.3612410	-1.1828850
H	-12.9860490	-5.2369700	-0.3846740
H	-12.4351540	-4.9026440	-2.0292980
O	-11.1607650	-3.2417710	-0.2999000
N	-9.8324980	-5.0865480	-0.2698060
C	-10.4696640	0.4137740	0.2327830
H	-11.2961560	-1.4121650	-0.2270840
C	-12.6701290	0.1241850	-0.7973300
H	-13.4143860	-0.6758990	-0.8221810
H	-13.0440240	0.9513010	-0.1864810
H	-12.5174400	0.4999460	-1.8172260
O	-10.5441040	1.6517190	0.2070240
N	-11.4380300	-0.4010890	-0.2295360
C	-9.2944840	-0.1235650	2.3859310
H	-8.4295720	-0.6208430	2.8362840
H	-9.2735240	0.9378990	2.6534360
H	-10.2069530	-0.5668120	2.8003770
C	-8.6577650	-4.2396490	1.7502730
H	-8.6336650	-5.2180650	2.2425560
H	-7.7934310	-3.6535310	2.0756270
H	-9.5686330	-3.7151830	2.0549720
C	-4.8510070	-5.3780900	-2.1362220
H	-4.8234760	-4.3958020	-2.6196580
H	-3.9266330	-5.9140860	-2.3785240
H	-5.6969610	-5.9469240	-2.5312710
C	-5.9042800	0.1485360	-2.2280000
H	-5.0622480	0.6767260	-2.6860530
H	-5.8818520	-0.8992260	-2.5462400
H	-6.8344010	0.6045100	-2.5836600
C	-2.3318620	-0.2048740	2.1647900
H	-1.4762740	-0.7260190	2.6051320
H	-2.3085600	0.8437250	2.4801780
H	-3.2514840	-0.6645970	2.5417130
C	-1.4505890	-4.2351800	2.2032580
H	-1.4402300	-5.2075090	2.7077800
H	-0.5513250	-3.6788900	2.4851210
H	-2.3265600	-3.6723850	2.5398170
C	1.2084370	0.1431420	-2.2348980
H	2.0687510	0.6581440	-2.6729380
H	1.2240240	-0.9053790	-2.5515490
H	0.2936350	0.6113960	-2.6127640
C	2.1678340	-5.4191260	-1.8838080
H	2.1713300	-4.4411320	-2.3765330
H	3.0858150	-5.9520780	-2.1558300
H	1.3115460	-5.9962120	-2.2435790
C	4.7259820	-0.2301460	2.1652300
H	4.7394600	0.8163680	2.4882550
H	3.8057580	-0.6992460	2.5281820
H	5.5803740	-0.7478100	2.6121780
C	5.7908110	-4.4318900	2.3710480
H	5.8352510	-5.4449340	2.7852560
H	6.7165260	-3.9052410	2.6251720
H	4.9521320	-3.8992030	2.8301290
C	8.3434130	0.1638330	-2.1723290
H	8.3529780	-0.8767470	-2.5156610
H	7.4436680	0.6531480	-2.5595980
H	9.2202430	0.6772510	-2.5782860
C	9.2039790	-4.6861970	-2.1918710
H	8.3515920	-5.0606630	-2.7666590
H	9.2624850	-3.5994590	-2.3085740
H	10.1210320	-5.1257480	-2.6008820
C	5.6318090	5.1729240	0.3425060
N	4.4075710	4.5111830	-0.0760440
H	4.4140500	3.4929010	-0.1603770
C	3.2525410	5.1938640	-0.2448400
O	3.1245190	6.4085560	-0.0491330
C	2.0830990	4.3446350	-0.7571220
N	0.8518750	5.0096750	-0.3555380
H	0.9261360	6.0112720	-0.2012790
C	6.8149030	4.3187400	-0.1254310
O	6.7144460	3.0853960	-0.2302550
N	7.9572970	4.9865370	-0.3714650
H	8.0306950	5.9886980	-0.2178220
C	9.2120440	4.3176200	-0.6874850
C	10.3405680	5.1430720	-0.0520090
O	10.2161480	6.3632990	0.1062790
N	11.4416660	4.4225410	0.2644590
H	11.3846840	3.4063750	0.1834180
H	5.6638410	6.1572810	-0.1369650
H	2.1163750	3.3575780	-0.2937810
H	9.1939910	3.3180170	-0.2526610
C	-1.4463960	5.1661300	0.4571630
N	-2.6854790	4.5386670	0.0291540
H	-2.6991850	3.5237310	-0.0883160
C	-3.8371180	5.2454630	-0.0705210
O	-3.9397740	6.4479650	0.1965500
C	-5.0278370	4.4448330	-0.6083010
N	-6.2437250	5.0037410	-0.0323850
H	-6.2336310	6.0006450	0.1577260
C	-7.3082520	4.2407290	0.2894640
O	-7.3315760	3.0158330	0.0861860
C	-0.2795870	4.3321110	-0.0799550
O	-0.3743010	3.1014400	-0.2106620
H	-1.4175790	6.1752960	0.0320660
H	-4.9470590	3.4035790	-0.2939180
C	12.6336650	5.0305130	0.8304500
H	12.5065250	6.1147780	0.8143620
H	12.7887680	4.7092560	1.8682830
H	13.5211780	4.7628740	0.2447340
C	-8.4844700	4.9516410	0.9667310
H	-8.4366830	6.0289730	0.7823290
H	-9.9135360	3.4906150	0.3452220
C	-10.6271700	5.3694540	-0.1642070
C	-11.8757410	4.7252790	-0.7487280
H	-11.8663650	3.6348450	-0.6706280
H	-12.7498380	5.1273910	-0.2259860
H	-11.9580830	5.0224750	-1.7991540
O	-10.4546010	6.5906200	-0.1807590
N	-9.7293970	4.4934900	0.3791160
C	-8.4308070	4.6840130	2.4836880
H	-8.5168580	3.6109140	2.6864060
H	-7.4892600	5.0434340	2.9149250
H	-9.2602170	5.2014700	2.9753120
C	-5.0482200	4.5011430	-2.1485310
H	-5.1013960	5.5381050	-2.4971650
H	-4.1457140	4.0367990	-2.5593350
H	-5.9212430	3.9580950	-2.5237590
C	-1.3460680	5.2784410	1.9937130
H	-1.3636510	4.2838040	2.4518490
H	-0.4216690	5.7858470	2.2915790
H	-2.1958930	5.8577250	2.3647610
C	2.1724380	4.1733440	-2.2860500
H	2.1925380	5.1478840	-2.7859180
H	3.0773700	3.6188290	-2.5528230
H	1.3031690	3.6123060	-2.6426250
C	5.6972950	5.3691830	1.8724590
H	5.6774400	4.4006380	2.3833790
H	6.6112750	5.8993270	2.1621610
H	4.8356260	5.9602630	2.1945880
C	9.4279440	4.1937920	-2.2074560
H	8.5918850	3.6473880	-2.6555880
H	9.4928860	5.1832340	-2.6728970
H	10.3524800	3.6455800	-2.4143040
H	2.0580890	-6.2447180	0.1081090
