%nproc=24
%mem=64GB
#P PBE1PBE/6-31+G(d,p) opt int=(grid=99590) freq EmpiricalDispersion=GD3

solv sp

0 1
C               -6.577414            5.848026            0.812462
N               -5.392049            6.641741            0.593856
H               -5.490126            7.664838            0.560540
C               -4.202570            6.099354            0.316693
O               -3.970268            4.884903            0.324815
C               -3.132531            7.103386           -0.086096
N               -1.835201            6.610809            0.319530
H               -1.672903            5.601110            0.313396
C               -7.681666            6.343832           -0.117769
O               -7.714198            7.534218           -0.471272
N               -8.589365            5.427039           -0.457786
H               -8.436852            4.451901           -0.181668
C               -9.755778            5.711411           -1.268640
C              -10.969951            5.099849           -0.572773
O              -10.959061            3.914353           -0.219790
N              -12.022024            5.904488           -0.402953
H              -11.961685            6.888713           -0.680297
H               -6.341283            4.808587            0.575052
H               -9.859998            6.796862           -1.335220
C                3.233490           -6.331929            0.625134
C                4.531791           -6.877036            0.039506
O                4.732789           -8.104301           -0.003639
N                5.397431           -5.949473           -0.368687
H                5.136345           -4.963904           -0.303133
C                6.702211           -6.237415           -0.922557
C                7.734060           -5.389676           -0.181043
O                7.475755           -4.221527            0.148195
N                8.905920           -5.979732            0.034565
H                9.027620           -6.963230           -0.253097
N                2.120682           -6.967040           -0.056699
H                2.220194           -7.941472           -0.382734
C                0.935816           -6.367015           -0.173289
O                0.688434           -5.231576            0.257212
C               -0.127727           -7.156273           -0.935612
N               -1.407296           -6.936526           -0.297358
H               -1.711490           -5.965734           -0.193014
H                3.170994           -5.256542            0.440515
H                6.912109           -7.297728           -0.773762
H                0.093213           -8.224958           -0.896668
C                0.447212            6.783996            1.096979
N                1.573090            7.629564            0.773134
H                1.464306            8.650005            0.829513
C                2.754631            7.120702            0.417336
O                3.007889            5.909474            0.384609
C                3.782815            8.149553           -0.029408
N                5.109441            7.625996            0.211565
H                5.250946            6.613498            0.183870
C                6.151932            8.433322            0.404857
O                6.088822            9.675746            0.377162
C               -0.829299            7.445686            0.596430
O               -0.911901            8.683013            0.516294
H                0.576591            5.827347            0.585450
H                3.661634            9.052914            0.572153
C               -3.489275           -7.442779            0.818834
C               -2.198068           -7.909325            0.154431
O               -1.924447           -9.122800            0.081191
N               -4.616106           -8.063974            0.146241
H               -4.549307           -9.067194           -0.084799
C               -5.807133           -7.463071            0.072183
O               -6.027802           -6.309516            0.465846
C               -6.918639           -8.287688           -0.582977
N               -8.162214           -8.059794            0.114069
H               -8.630423           -7.167691           -0.043487
C               -8.552624           -8.805913            1.162405
O               -7.938473           -9.807170            1.555071
H               -3.590087           -6.361792            0.699587
H               -6.679108           -9.352583           -0.506005
C               -9.829806           -8.375752            1.837219
H              -10.247587           -7.462142            1.409025
H              -10.554532           -9.190898            1.756951
H               -9.622130           -8.226194            2.900210
C              -13.233797            5.419731            0.218800
H              -13.695244            4.629580           -0.383043
H              -13.022941            5.007152            1.210599
H              -13.932037            6.252230            0.315615
C               10.002231           -5.288283            0.682509
H                9.990081           -4.243371            0.361906
C                7.469812            7.742913            0.724753
H                7.591173            6.878618            0.065661
H                8.262790            9.659375            0.415765
C                9.793530            8.288195            0.218495
C               10.803467            9.366891           -0.076186
H               10.382163           10.374837           -0.062134
H               11.608397            9.302425            0.660611
H               11.241249            9.174627           -1.059439
O               10.117719            7.094714            0.242721
N                8.533791            8.684979            0.459448
C               11.305117           -5.950155            0.245079
H               12.085094           -4.096600           -0.081866
C               13.615659           -5.542650           -0.424026
H               14.376191           -4.944551            0.086388
H               13.722388           -6.589971           -0.136023
H               13.772092           -5.462397           -1.506320
O               11.419088           -7.178996            0.219065
N               12.297257           -5.092744           -0.047041
C                9.877021           -5.342753            2.206359
H                8.927805           -4.890961            2.506964
H                9.908838           -6.379610            2.553376
H               10.694217           -4.785912            2.675647
C                7.482659            7.249663            2.174101
H                7.298241            8.076507            2.867024
H                6.721126            6.478871            2.313826
H                8.461499            6.814620            2.392388
C                3.564312            8.491876           -1.506478
H                3.651690            7.590407           -2.120793
H                2.572927            8.928698           -1.652659
H                4.313042            9.218487           -1.831397
C                6.748847           -5.916922           -2.417412
H                5.998159           -6.510881           -2.945876
H                6.541801           -4.855198           -2.586126
H                7.735510           -6.156550           -2.823942
C                3.193721           -6.591727            2.130919
H                2.254251           -6.207262            2.536831
H                3.262718           -7.663862            2.334933
H                4.025021           -6.083316            2.629872
C                0.344236            6.519774            2.602121
H                0.193703            7.455603            3.149056
H               -0.487213            5.841535            2.814266
H                1.270882            6.053711            2.947753
C               -0.160052           -6.695405           -2.393529
H               -0.928794           -7.248880           -2.939357
H               -0.388056           -5.626153           -2.451095
H                0.807929           -6.875278           -2.869933
C               -3.203429            7.351552           -1.595932
H               -3.031258            6.418038           -2.140308
H               -4.183782            7.749481           -1.873205
H               -2.438338            8.077178           -1.884086
C               -3.467158           -7.800296            2.305154
H               -3.412632           -8.883798            2.439685
H               -2.610228           -7.333491            2.801786
H               -4.385522           -7.436180            2.773320
C               -7.052227            5.927232            2.265272
H               -7.294605            6.959396            2.537695
H               -7.936773            5.300964            2.412642
H               -6.257208            5.569985            2.925471
C               -7.050800           -7.887921           -2.050179
H               -7.293598           -6.823917           -2.137403
H               -6.113348           -8.075035           -2.580938
H               -7.846586           -8.467586           -2.525580
C               -9.607747            5.115802           -2.667995
H               -8.728485            5.540072           -3.160222
H               -9.490722            4.029938           -2.603457
H              -10.491257            5.333915           -3.276055
H               -3.315005            8.041706            0.440811
O               -2.901407           10.507614            0.111592
H               -2.710195           11.192963            0.757931
H               -2.194882            9.832190            0.225568
O               -5.451972            9.490928            0.403696
H               -6.016513            9.633978           -0.362125
H               -4.595943            9.936802            0.231972
O                4.052921           11.496459            0.008212
H                4.352885           12.304740            0.432373
H                4.784762           10.850790            0.116136
O                1.611117           10.491751            0.818912
H                0.862269           10.749227            0.273951
H                2.405359           10.953759            0.480305
O                2.383618           -9.622025           -0.883537
H                1.739521          -10.144421           -0.361176
H                3.240029           -9.652076           -0.438537
O                0.275214          -10.730407            0.451039
H               -0.491821          -10.133107            0.329021
H               -0.037839          -11.607256            0.213362
O                9.128489           -8.631907           -0.822740
H                9.978784           -8.704633           -0.368397
H                8.468701           -9.193165           -0.371014
O                6.862831           -9.858236            0.167721
H                6.126458           -9.214072            0.120963
H                6.580346          -10.613803           -0.354478
O               -9.389333            9.609924           -0.384985
H               -9.233777           10.034511            0.462892
H               -8.775871            8.841707           -0.424412
O              -11.805134            8.655831           -1.051125
H              -12.078967            9.079184           -1.868196
H              -10.958078            9.084354           -0.783162
O               -5.562639          -10.988437            2.117366
H               -6.399289          -10.520644            1.904416
H               -5.827943          -11.787921            2.578576
O               -4.197051          -10.816869           -0.237210
H               -3.285353          -10.544363           -0.050314
H               -4.599359          -11.061936            0.621808
C               -5.101826           -3.367854            0.624365
N               -3.999449           -2.607506            0.077696
H               -4.114705           -1.599850           -0.041482
C               -2.830347           -3.180385           -0.222108
O               -2.594703           -4.383256           -0.046817
C               -1.795812           -2.268592           -0.863939
N               -0.479140           -2.680621           -0.426680
H               -0.310557           -3.668874           -0.220857
C               -6.390937           -2.886722           -0.025225
O               -6.548407           -1.692227           -0.315743
N               -7.337495           -3.816453           -0.209196
H               -7.120590           -4.797778           -0.012339
C               -8.626843           -3.464615           -0.774184
C               -9.701478           -4.344998           -0.151378
O               -9.572822           -5.575159           -0.096156
N              -10.795487           -3.711387            0.279538
H              -10.843914           -2.693517            0.219268
H               -4.938291           -4.420695            0.387845
H               -8.821681           -2.418669           -0.526086
C                1.814015           -2.359774            0.265870
N                2.921128           -1.557663           -0.202071
H                2.793470           -0.550551           -0.310350
C                4.113737           -2.112056           -0.443980
O                4.340534           -3.318428           -0.286342
C                5.181992           -1.176913           -0.988736
N                6.479417           -1.637758           -0.542443
H                6.615724           -2.635203           -0.350504
C                7.494966           -0.781732           -0.384276
O                7.400129            0.429550           -0.635938
C                0.517760           -1.798304           -0.295040
O                0.398827           -0.594030           -0.563035
H                1.954330           -3.378481           -0.100919
H                5.016159           -0.173895           -0.590276
C              -11.915646           -4.444651            0.822825
H              -12.223345           -5.240186            0.137307
H              -11.658759           -4.907364            1.782228
H              -12.746185           -3.753124            0.972425
C                8.771303           -1.376492            0.189476
H                8.884628           -2.398907           -0.177820
H                9.801162            0.395990           -0.412777
C               11.123151           -1.173941           -0.366892
C               12.252012           -0.287251           -0.824134
H               11.941578            0.745864           -0.989186
H               13.043154           -0.321194           -0.069603
H               12.663570           -0.706208           -1.747030
O               11.315099           -2.372244           -0.111224
N                9.909988           -0.606243           -0.259570
C                8.678683           -1.430876            1.718002
H                8.522562           -0.429280            2.132463
H                7.856314           -2.083064            2.024046
H                9.611690           -1.838760            2.116032
C                5.090476           -1.109113           -2.515737
H                5.218132           -2.103855           -2.954289
H                4.120467           -0.705099           -2.820163
H                5.878044           -0.450907           -2.892258
C                1.749318           -2.409715            1.795242
H                1.579048           -1.410322            2.209231
H                0.946374           -3.077780            2.117985
H                2.697023           -2.795436            2.181259
C               -1.933876           -2.317237           -2.388043
H               -1.800025           -3.339484           -2.755906
H               -2.921282           -1.957098           -2.692216
H               -1.172357           -1.674837           -2.838422
C               -5.197288           -3.214767            2.144200
H               -5.370729           -2.169029            2.418762
H               -6.013916           -3.827400            2.537090
H               -4.261775           -3.549897            2.600598
C               -8.639330           -3.625288           -2.295021
H               -7.854801           -3.001270           -2.731659
H               -8.463007           -4.668722           -2.572272
H               -9.604537           -3.311998           -2.705398
H               -1.961583           -1.241327           -0.533474
C                5.409901            2.952769           -0.841051
N                4.122705            3.345531           -0.309950
H                3.968024            4.324630           -0.067991
C                3.128633            2.460252           -0.173280
O                3.227783            1.275197           -0.515553
C                1.868909            2.983388            0.499278
N                0.735105            2.216612            0.034957
H                0.866608            1.227103           -0.186266
C                6.488889            3.795193           -0.177311
O                6.268852            4.967860            0.159249
N                7.675929            3.194946           -0.039951
H                7.767224            2.200436           -0.268762
C                8.824799            3.871443            0.523954
C               10.059385            3.477035           -0.278128
O               10.259401            2.301525           -0.602146
N               10.915890            4.465413           -0.559586
H               10.686996            5.421095           -0.293500
H                5.574323            1.900710           -0.597935
H                8.654595            4.947496            0.442823
C               -1.552648            1.857176           -0.648578
N               -2.847693            2.291085           -0.176115
H               -3.004197            3.285730           -0.008710
C               -3.845769            1.417470           -0.005158
O               -3.734120            0.209081           -0.250686
C               -5.129290            1.988853            0.574594
N               -6.240747            1.152951            0.178138
H               -6.084866            0.154934            0.018203
C               -7.466507            1.668550            0.039083
O               -7.723190            2.858866            0.274837
C               -0.480740            2.764552           -0.065194
O               -0.725266            3.938901            0.246803
H               -1.383648            0.836967           -0.297297
H               -5.289224            2.987592            0.162459
C               -8.525573            0.712936           -0.485738
H               -8.377977           -0.264886           -0.021204
H              -10.012533            2.195967           -0.093244
C              -10.843754            0.329739            0.101043
C              -12.179491            0.919832            0.468250
H              -12.182244            2.011001            0.444637
H              -12.930479            0.531161           -0.225086
H              -12.446207            0.568608            1.469588
O              -10.687578           -0.898474            0.011385
N               -9.837264            1.191186           -0.113129
C               -8.382377            0.541527           -2.001576
H               -8.479983            1.504681           -2.513499
H               -7.412527            0.098927           -2.242800
H               -9.168750           -0.129014           -2.358805
C               -5.015496            2.111610            2.096981
H               -4.814602            1.135345            2.549662
H               -4.211746            2.803996            2.361879
H               -5.957255            2.499019            2.495284
C               -1.469965            1.852472           -2.177750
H               -1.601140            2.863506           -2.576636
H               -0.503459            1.458231           -2.503377
H               -2.259521            1.209083           -2.575493
C                2.022435            2.910396            2.021625
H                2.208236            1.880173            2.341948
H                2.851742            3.542793            2.351341
H                1.102224            3.266293            2.492848
C                5.477210            3.113014           -2.362057
H                5.336217            4.160406           -2.647420
H                6.443832            2.765859           -2.738569
H                4.687929            2.510979           -2.820181
C                9.024394            3.498602            1.994296
H                8.132552            3.771145            2.565421
H                9.193455            2.421953            2.093374
H                9.885078            4.029523            2.412983
H                1.716595            4.026876            0.216907
C               12.155842            4.198449           -1.252236
H               12.675920            5.143724           -1.414023
H               12.797945            3.532175           -0.666222
H               11.966206            3.719241           -2.217959
