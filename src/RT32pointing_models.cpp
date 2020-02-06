#include "cpeds-math.h"
#include <math.h>
#include <RT32pointing_models.h>


/*****************************************************************************************************************/
void cpeds_RT32_Model4(double AZ, double ZD, double *dAZ, double *dZD) {
	
	//const double p[13] = {
	//-0.032367, -0.049374, 0.051222, -0.054595, -0.000161, 0.002850,
	//0.001370, 0.020659, -0.012753, 0.008571, 0.004812, -0.006485, 0.028568 };
	
	const double p[13] = { -3.2367E-02,-4.9374E-02, 5.1222E-02,-5.4595E-02,
			-1.6113E-04, 2.8503E-03, 1.3701E-03, 2.0659E-02,-1.2753E-02,
			8.5710E-03, 4.8125E-03,-6.4854E-03, 2.8568E-02 };
	
	double pi180 = M_PI/180.;
	
	AZ = AZ * pi180;
	double Alt = (90.-ZD) * M_PI/180.;
	
	double xi = p[4]*pi180;
	double zeta = p[5]*pi180;
	double sigma = p[2]*pi180;
	double beta = p[3]*pi180;
	double sh = sin(Alt);
	double ch = cos(Alt);
	
	double sinxi = sin(xi);
	double sinzeta = sin(zeta);
	double se = sqrt(sinxi*sinxi + sinzeta*sinzeta);
	
	if (ch < 0.) { se = -se; }
	
	double ce = sqrt(1.-se*se);
	double alfa = atan2(sin(zeta), sin(xi));
	double AT = atan2(ch*sin(alfa-AZ), sh*se-ch*ce*cos(alfa-AZ)) -
			atan2(sin(alfa), -ce*cos(alfa));
	
	double hT = asin(ce*sh + se*ch*cos(alfa-AZ));
	double AT_Az = AT - AZ;
	if (ch < 0.) {
		AT_Az = M_PI - AT_Az;
		hT = M_PI - hT;
	}
	
	*dAZ = cpeds_RT4_arcsin((sin(sigma)*sin(hT) + sin(beta)) /
			(cos(hT)*cos(sigma)));
	
	double hb = atan2(sin(hT)*cos(sigma) + cos(hT)*sin(sigma)*sin(*dAZ),
			cos(hT)*cos(*dAZ));
	
	*dAZ = fmod((*dAZ + AT_Az)*180./M_PI, 360.) +
			p[0] + 
			p[8]*sin(2.*AZ) +
			p[9]*cos(2.*AZ) +
			p[11]*sin(3.*AZ)*sin(Alt) + 
			p[12]*cos(Alt)*cos(AZ/4.);
	
	*dZD = (hb-Alt) * 180./M_PI +
			p[1] +
			p[6]*ch +
			p[7]*sh +
			p[10]*sin(2.*AZ);
	*dZD = -(*dZD);
}
/***************************************************************************************/
void cpeds_RT32_Model4e(double AZ, double ZD, double *dAZ, double *dZD) {
	const double p[9] = { 
			3.4144712346E-03, -9.1488156466E-04, -2.2984989887E-03, 1.8860129779E-02, 
			-4.3813651976E-02, -1.0356952146E-02, 8.2631028153E-03, -7.4978342202E-03,
			5.0519549549E-03
	}; // deg
	const double q[7] = { 
			8.2681943479E-02, 4.2711366808E-05, 2.7049249587E-04, -8.7588060472E-04, 
			-3.4899753854E-02, -4.1424115143E-03, 3.6971967552E-03
	}; // deg
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	double MA=0,MZ=0;
	double A=AZ*PI180;
	double Z=ZD*PI180;
	double sinA=sin(A);
	double cosA=cos(A);
	double sin2A=sin(2.0*A);
	double cos2A=cos(2.0*A);
	double sinZ=sin(Z);
	double cosZ=cos(Z);
	
	// Model4e
	MA+=p[0]+( p[1]*sinA-p[2]*cosA+p[3] )/tan(Z) + p[4]/sinZ;
	MA+=p[5]*sin2A + p[6]*cos2A + p[7]*sin(3.0*A)*cosZ + p[8]*cos(A/4)*sinZ;
	
	MZ+=q[0]+q[1]*cosA+q[2]*sinA+q[3]*sinZ+q[4]*cosZ+q[5]*sin2A+q[6]*cos2A;
	
	*dAZ=MA;
	*dZD=MZ;
}
/***************************************************************************************/
/***************************************************************************************/
void cpeds_RT32_Model5(double AZ, double ZD, double *dAZ, double *dZD) {
	const double deltaA_A_phi_T[51][3] = {
			{9.1464679913E-06,0.0,0.0},
			{7.4812314678E-04,-1.0146440002E+00,3.5964783888E+02},
			{2.1099373994E-04,2.1395266021E+00,1.7982391944E+02},
			{7.5789593125E-04,3.7969801449E-01,1.1988261296E+02},
			{5.6486694523E-04,-2.5560328871E-01,8.9911959720E+01},
			{2.9594479076E-04,2.9213364814E+00,7.1929567776E+01},
			{1.5902554032E-03,2.0405774586E+00,5.9941306480E+01},
			{5.1419002198E-04,-2.6892085396E+00,5.1378262697E+01},
			{2.4555827774E-04,-2.6879404830E+00,4.4955979860E+01},
			{7.0900164890E-04,-2.4479912145E+00,3.9960870987E+01},
			{1.1665254860E-03,-1.8508451646E+00,3.5964783888E+01},
			{2.1747761860E-04,5.0597093916E-01,3.2695258080E+01},
			{3.4186615678E-04,9.6981714032E-01,2.9970653240E+01},
			{2.2048263759E-04,-1.0939918206E+00,2.7665218375E+01},
			{3.4358437983E-03,-2.6928190457E+00,2.5689131349E+01},
			{4.6753622473E-04,-2.8479383897E+00,2.3976522592E+01},
			{7.7632348788E-05,-1.4847739440E+00,2.2477989930E+01},
			{6.8024460573E-05,-2.2613095138E+00,2.1155755228E+01},
			{1.5780593336E-04,1.8409436599E+00,1.9980435493E+01},
			{1.5954792635E-04,-2.4045974145E+00,1.8928833625E+01},
			{1.5467645909E-04,-2.6642158778E+00,1.7982391944E+01},
			{1.8160142754E-04,2.3099465958E-01,1.7126087566E+01},
			{1.0035780153E-04,7.1399949804E-01,1.6347629040E+01},
			{1.0490219191E-04,5.5218331689E-01,1.5636862560E+01},
			{2.7213169735E-04,1.9217781072E+00,1.4985326620E+01},
			{2.8771359850E-05,2.2834929006E+00,1.4385913555E+01},
			{2.7055732417E-04,2.9936859119E+00,1.3832609188E+01},
			{1.1862895918E-04,-2.5581520809E+00,1.3320290329E+01},
			{7.9562328734E-05,-2.6941803972E+00,1.2844565674E+01},
			{1.8689729643E-04,-2.3185283681E+00,1.2401649617E+01},
			{8.0638670934E-05,-1.6240672909E+00,1.1988261296E+01},
			{1.8136793028E-04,-1.2375765810E+00,1.1601543190E+01},
			{1.6136708224E-04,-2.3978970641E+00,1.1238994965E+01},
			{1.2495378166E-04,6.0103896559E-01,1.0898419360E+01},
			{2.4145259184E-04,1.5398138470E+00,1.0577877614E+01},
			{2.2646069620E-04,2.0528858146E+00,1.0275652539E+01},
			{1.6499465983E-04,2.4674367246E+00,9.9902177466E+00},
			{1.2696608887E-04,-2.2342792246E+00,9.7202118616E+00},
			{1.8477887487E-04,-2.4500988264E+00,9.4644168126E+00},
			{2.1960431737E-04,-6.3960034023E-01,9.2217394584E+00},
			{1.0188446233E-04,-4.9474299296E-01,8.9911959720E+00},
			{2.1053922078E-04,-2.9518911738E+00,8.7718985092E+00},
			{1.4852096038E-03,2.2113244978E+00,8.5630437828E+00},
			{1.7162835158E-04,2.0271244979E+00,8.3639032297E+00},
			{1.1407557322E-04,1.5943235394E+00,8.1738145200E+00},
			{8.0313176267E-05,-2.9625652152E-01,7.9921741973E+00},
			{1.3482542279E-04,2.4574276126E+00,7.8184312800E+00},
			{1.3187051810E-04,-2.9775507141E+00,7.6520816783E+00},
			{1.3109424134E-04,-2.0719697657E+00,7.4926633100E+00},
			{1.3625561070E-04,-1.4991715456E+00,7.3397518139E+00},
			{6.3555475863E-05,1.5036934663E+00,7.1929567776E+00}};
	
	const double deltaZ_A_phi_T[51][3] = {
			{6.1031365843E-05,0.0,0.0},
			{5.8974717125E-04,1.7811473351E+00,3.5990974471E+02},
			{4.0854852247E-04,-2.6003551057E+00,1.7995487235E+02},                                                                                                 
			{3.4867701410E-03,1.5128086420E+00,1.1996991490E+02},                                                                                                  
			{7.8911821183E-04,1.9672588206E+00,8.9977436177E+01},                                                                                                  
			{6.9451672748E-04,-9.5592799286E-01,7.1981948942E+01},                                                                                                 
			{9.4651234414E-04,1.6532224066E+00,5.9984957452E+01},
			{4.0666916301E-04,2.6466373050E+00,5.1415677816E+01},
			{3.2267969308E-04,2.9303803378E+00,4.4988718089E+01},
			{1.1016303575E-03,-1.6497812365E+00,3.9989971634E+01},
			{1.0202613956E-03,-2.2047342894E+00,3.5990974471E+01},
			{2.3225751584E-04,-1.2424765126E+00,3.2719067701E+01},
			{2.5049337354E-04,2.2845075865E+00,2.9992478726E+01},
			{1.8276649393E-04,-2.6514459471E-01,2.7685364978E+01},
			{2.1991151848E-03,-2.6982071766E+00,2.5707838908E+01},
			{4.5762406162E-04,-3.9979919812E-01,2.3993982981E+01},
			{2.8537630491E-04,-4.3498146135E-01,2.2494359044E+01},
			{2.5601999128E-04,-9.5130402198E-01,2.1171161454E+01},
			{3.4131811378E-04,-2.6255431693E+00,1.9994985817E+01},
			{4.0785879696E-04,-2.9737380475E+00,1.8942618143E+01},
			{7.0419095424E-04,-2.9510473627E+00,1.7995487235E+01},
			{2.1776010725E-04,3.0084577024E+00,1.7138559272E+01},
			{3.6958604739E-05,-2.1961016080E+00,1.6359533850E+01},
			{3.6834571411E-04,4.9972426188E-02,1.5648249770E+01},
			{2.8645969469E-04,3.1083111914E-01,1.4996239363E+01},
			{3.0475897324E-04,9.4605913134E-02,1.4396389788E+01},
			{1.3821256068E-04,-1.2051938560E+00,1.3842682489E+01},
			{1.8788813573E-04,2.7052042092E+00,1.3329990545E+01},
			{1.5646459375E-04,-2.2659842917E+00,1.2853919454E+01},
			{4.7487733425E-04,-2.9503798986E+00,1.2410680852E+01},
			{1.7631710246E-04,-2.6475751136E+00,1.1996991490E+01},
			{8.5467498390E-05,-2.4013118071E+00,1.1609991765E+01},
			{3.3306771678E-04,-4.8102299618E-01,1.1247179522E+01},
			{3.2772469936E-04,-3.0997894318E-02,1.0906355900E+01},
			{3.2294873161E-04,2.0693484746E-01,1.0585580727E+01},
			{1.6459431559E-04,1.7563465081E+00,1.0283135563E+01},
			{2.0195800309E-04,-2.6410653936E+00,9.9974929086E+00},
			{3.8195320183E-04,-2.9463440179E+00,9.7272903976E+00},
			{3.3587917938E-04,-2.5961626323E+00,9.4713090713E+00},
			{1.7860306368E-04,-3.0561297896E+00,9.2284549926E+00},
			{3.1128186374E-04,2.3303732166E-01,8.9977436177E+00},
			{9.5328259266E-05,-3.1221250565E+00,8.7782864563E+00},
			{7.4697876787E-04,1.8867612798E+00,8.5692796359E+00},
			{2.0636941761E-04,-2.3087059333E+00,8.3699940630E+00},
			{1.4737174254E-04,-6.7120776279E-01,8.1797669252E+00},
			{2.8705226183E-05,-1.0659823844E+00,7.9979943269E+00},
			{9.9876318620E-05,9.8926442465E-01,7.8241248850E+00},
			{2.0844522051E-04,-2.9464022958E+00,7.6576541428E+00},
			{2.3603680020E-04,2.4934551141E+00,7.4981196815E+00},
			{4.7268543178E-05,2.2828282702E+00,7.3450968308E+00},
			{1.6503420594E-04,6.5779151054E-01,7.1981948942E+00}};
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	double MA=0,MZ=0;
	
	// Model4e
	cpeds_RT32_Model4e(AZ,ZD,&MA,&MZ);
	
	// Model5
	MA+=deltaA_A_phi_T[0][0];
	MZ+=deltaZ_A_phi_T[0][0];
	for (unsigned long i = 1; i < 51; i++) {
		MA+=deltaA_A_phi_T[i][0]*sin(-twoPI/deltaA_A_phi_T[i][2]*(AZ-180)+deltaA_A_phi_T[i][1]);
		MZ+=deltaZ_A_phi_T[i][0]*sin(-twoPI/deltaZ_A_phi_T[i][2]*(AZ-180)+deltaZ_A_phi_T[i][1]);
	}
	
	*dAZ=MA;
	*dZD=MZ;
}

/***************************************************************************************/
/***************************************************************************************/
/* ******************************************************************************************** */
void Model4g_DR2017sum(double AZ, double ZD, double* dAZ, double *dZD) {
	static const double p[9] = { 
			1.3055302387E-02, 6.5695042094E-04, -1.1498839756E-03, 2.5331659102E-02, 
			-5.3689330650E-02, -9.9878337964E-03, 8.1609365701E-03, -8.8716819295E-03,
			-2.0867022803E-03 
	}; // deg
	static const double q[7] = { 
			9.3408029523E-02, 3.3243649853E-03, -1.6482803900E-03, 3.6375218824E-02, 
			8.3740866915E+01, -3.5432521791E-03, 3.1894010297E-03 
	}; // deg

	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;

	double PI180=M_PI/180l;

	double MA=0,MZ=0;
	double A=AZ*PI180;
	double Z=ZD*PI180;
	double sinA=sin(A);
	double cosA=cos(A);
	double sin2A=sin(2.0*A);
	double cos2A=cos(2.0*A);
	double sinZ=sin(Z);
	double cosZ=cos(Z);
	double phi=q[4]*PI180;
	
	// Model4g
	MA+=p[0]+( p[1]*sinA-p[2]*cosA+p[3] )/tan(Z) + p[4]/sinZ;
	MA+=p[5]*sin2A + p[6]*cos2A + p[7]*sin(3.0*A)*cosZ + p[8]*cos(A/4)*sinZ;
	
	MZ+=q[0]+q[1]*cosA+q[2]*sinA+q[3]*sin(Z-phi)+q[5]*sin2A+q[6]*cos2A;
	
	*dAZ=MA;
	*dZD=MZ;	
}

/* ******************************************************************************************** */
void Model4gr_DR2017sum(double AZ, double ZD, double *dAZ, double *dZD) {
	static const double deltaA_A_phi_T[51][3] = {
	{-4.4193510922E-05,0.0,0.0},
	{2.0975624936E-04,2.3739467890E+00,3.5995884719E+02},
	{2.1751766083E-05,1.5495360203E+00,1.7997942359E+02},
	{1.0854141784E-03,3.8681293724E-01,1.1998628240E+02},
	{2.8992718338E-04,-5.6802627748E-01,8.9989711797E+01},
	{4.4884800373E-04,-2.3100821930E+00,7.1991769438E+01},
	{1.3836783637E-03,2.1806533044E+00,5.9993141198E+01},
	{2.2485302979E-04,-2.0194550302E+00,5.1422692456E+01},
	{2.9168988970E-04,-1.3217108418E+00,4.4994855899E+01},
	{4.9852469088E-04,-2.1598779336E+00,3.9995427466E+01},
	{1.0794846894E-03,-1.7852151662E+00,3.5995884719E+01},
	{2.5641483632E-04,7.0085045918E-01,3.2723531563E+01},
	{1.8949331692E-04,1.8596130245E+00,2.9996570599E+01},
	{2.1418246633E-04,-1.7356228089E+00,2.7689142092E+01},
	{3.7778697518E-03,-2.3768899200E+00,2.5711346228E+01},
	{2.4361396666E-04,-2.5090394999E+00,2.3997256479E+01},
	{5.8319735260E-05,5.3737954663E-01,2.2497427949E+01},
	{1.0606995216E-04,1.1037184595E+00,2.1174049835E+01},
	{2.2855964603E-04,9.3851112570E-01,1.9997713733E+01},
	{7.0998829811E-05,2.4093838268E+00,1.8945202484E+01},
	{6.6321461852E-05,-2.9284123918E+00,1.7997942359E+01},
	{7.2439121418E-05,9.6744771972E-01,1.7140897485E+01},
	{1.3421645214E-04,-1.2826209016E+00,1.6361765781E+01},
	{7.1669917803E-05,-1.2539991715E+00,1.5650384660E+01},
	{1.3419292877E-04,2.6903416019E+00,1.4998285300E+01},
	{7.2135984931E-05,1.7750133268E+00,1.4398353888E+01},
	{1.4133007721E-04,1.2725667229E+00,1.3844571046E+01},
	{1.3074342218E-04,-9.2603046603E-01,1.3331809155E+01},
	{1.8804772953E-04,1.0562022024E+00,1.2855673114E+01},
	{1.4432651780E-04,-2.1050951511E+00,1.2412374041E+01},
	{1.5301098036E-04,4.8135220440E-01,1.1998628240E+01},
	{6.4849844692E-05,-2.3111608408E+00,1.1611575716E+01},
	{9.8254637506E-05,-2.0511060374E+00,1.1248713975E+01},
	{1.4873687336E-04,2.2041035384E+00,1.0907843854E+01},
	{2.0465524855E-04,2.1426705020E+00,1.0587024917E+01},
	{5.2404953408E-05,1.6013987175E-01,1.0284538491E+01},
	{7.3192280446E-05,-5.3187502061E-01,9.9988568664E+00},
	{1.1610973398E-04,2.5834079595E-03,9.7286174916E+00},
	{1.7855434004E-04,-1.2104279212E+00,9.4726012418E+00},
	{1.3093863735E-04,-6.2825321514E-01,9.2297140305E+00},
	{7.8778521918E-05,1.2796068420E+00,8.9989711797E+00},
	{2.2924306996E-04,-2.4347101338E+00,8.7794840778E+00},
	{1.5743712910E-03,3.1154526452E+00,8.5704487426E+00},
	{1.4934451616E-04,2.8588683103E+00,8.3711359812E+00},
	{7.7648989107E-05,-5.3312799189E-01,8.1808828907E+00},
	{1.4690714922E-04,6.9852473964E-01,7.9990854931E+00},
	{1.1152735371E-04,-1.2362716023E+00,7.8251923302E+00},
	{1.4563583864E-04,7.6110026708E-01,7.6586988764E+00},
	{3.2630776628E-05,2.7724211399E+00,7.4991426498E+00},
	{5.1016035439E-05,-1.0547400732E+00,7.3460989222E+00},
	{7.7384865383E-05,2.1933886797E+00,7.1991769438E+00}};

	
	static const double deltaZ_A_phi_T[51][3] = {
	{-2.5743329454E-05,0.0,0.0},
	{1.0888419147E-03,1.9703846671E+00,3.5984000000E+02},
	{3.2779631257E-04,2.3228829176E-01,1.7992000000E+02},
	{3.3538375747E-03,1.0683871481E+00,1.1994666667E+02},
	{8.8579721532E-04,-9.2840596677E-01,8.9960000000E+01},
	{1.3097960608E-03,-1.0716042806E+00,7.1968000000E+01},
	{1.0329489375E-03,1.8140791721E+00,5.9973333333E+01},
	{7.0738376937E-04,2.5416556578E+00,5.1405714286E+01},
	{8.1513070404E-04,8.2274978759E-01,4.4980000000E+01},
	{7.2022419665E-04,-1.0600993711E+00,3.9982222222E+01},
	{1.3469141474E-03,-1.9495153480E+00,3.5984000000E+01},
	{4.2439036183E-04,-2.1877884209E+00,3.2712727273E+01},
	{4.1704499696E-04,-7.9244798144E-01,2.9986666667E+01},
	{1.0719289466E-04,-1.1349221839E+00,2.7680000000E+01},
	{2.7398256578E-03,-2.3547655158E+00,2.5702857143E+01},
	{5.2917906371E-04,2.3497715074E+00,2.3989333333E+01},
	{2.1367121957E-04,9.7712116989E-01,2.2490000000E+01},
	{2.2857817647E-04,3.0649068012E+00,2.1167058824E+01},
	{4.1003393592E-04,2.9838571439E+00,1.9991111111E+01},
	{7.9410225404E-05,2.3451489537E+00,1.8938947368E+01},
	{1.6587716128E-04,2.3072980750E+00,1.7992000000E+01},
	{1.5694209459E-04,-2.0978029399E+00,1.7135238095E+01},
	{2.4834484118E-04,1.9590465472E+00,1.6356363636E+01},
	{1.0659409485E-04,2.4764260337E+00,1.5645217391E+01},
	{2.8083694079E-04,5.5179511681E-01,1.4993333333E+01},
	{1.7226090363E-04,-2.5552958646E+00,1.4393600000E+01},
	{2.9180362474E-04,-1.5541185005E+00,1.3840000000E+01},
	{2.1454235714E-04,1.7581024072E+00,1.3327407407E+01},
	{3.2008877744E-04,-2.7356621638E+00,1.2851428571E+01},
	{1.9970436010E-04,2.8271290519E+00,1.2408275862E+01},
	{2.6278746174E-04,1.6332309892E+00,1.1994666667E+01},
	{2.7431926955E-04,-4.4190087375E-01,1.1607741935E+01},
	{2.1523781205E-04,-4.0257281499E-02,1.1245000000E+01},
	{3.5153204247E-04,-7.7181648138E-01,1.0904242424E+01},
	{1.1537338779E-04,1.5311899800E+00,1.0583529412E+01},
	{1.9568664703E-04,2.7246847996E+00,1.0281142857E+01},
	{2.0224506359E-05,-1.5975164302E+00,9.9955555556E+00},
	{1.3388223419E-04,1.6728302766E+00,9.7254054054E+00},
	{2.4921279610E-04,-4.6299939896E-01,9.4694736842E+00},
	{8.0020755040E-05,2.3434612356E+00,9.2266666667E+00},
	{2.7650412924E-04,-1.6456585032E+00,8.9960000000E+00},
	{3.1845435349E-04,-1.7015699433E-01,8.7765853659E+00},
	{1.0407394309E-03,3.0915556870E+00,8.5676190476E+00},
	{1.1406159625E-04,-1.4955281005E+00,8.3683720930E+00},
	{1.1850617820E-04,8.9402413206E-01,8.1781818182E+00},
	{2.4236626741E-04,-2.8235201264E+00,7.9964444444E+00},
	{1.3456229612E-04,-7.2969835810E-01,7.8226086957E+00},
	{2.1768170628E-04,-1.8423576771E+00,7.6561702128E+00},
	{2.8231451302E-04,1.3754856555E+00,7.4966666667E+00},
	{2.0760105789E-04,1.2682457793E+00,7.3436734694E+00},
	{5.2762227416E-05,-1.7212947719E+00,7.1968000000E+00}};

	double twoPI=M_PI*2.0;

	double MA=0,MZ=0;

	int i;
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	// Model4g
	Model4g_DR2017sum(AZ,ZD,&MA,&MZ);

	// Model4gr
	MA+=deltaA_A_phi_T[0][0];
	MZ+=deltaZ_A_phi_T[0][0];
	for (i = 1; i < 51; i++) {
		MA+=deltaA_A_phi_T[i][0]*sin(-twoPI/deltaA_A_phi_T[i][2]*(AZ-180)+deltaA_A_phi_T[i][1]);
		MZ+=deltaZ_A_phi_T[i][0]*sin(-twoPI/deltaZ_A_phi_T[i][2]*(AZ-180)+deltaZ_A_phi_T[i][1]);
	}
	
	*dAZ=MA;
	*dZD=MZ;
}

/* ******************************************************************************************** */
void Model6a_DR2017(double AZ, double ZD, double *dAZ, double *dZD, double Tmeteo) {
	static const double p[9] = { 
			1.3055302387E-02, 6.5695042094E-04, -1.1498839756E-03, 2.5331659102E-02, 
			-5.3689330650E-02, -9.9878337964E-03, 8.1609365701E-03, -8.8716819295E-03,
			-2.0867022803E-03 
	}; // deg
	static const double q[8] = { 
			7.8536724219E-02, 1.2478175768E-03, -1.4147289865E-03, 2.9390953980E-02, 
			7.5092710705E+01, -3.5510875355E-03, 3.2233855803E-03, 1.9587411237E+02
	}; // deg,deg,deg,deg,deg,deg,deg,um/degC
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;

	double PI180=M_PI/180l;

	double MA=0,MZ=0;
	double A=AZ*PI180;
	double Z=ZD*PI180;
	double sinA=sin(A);
	double cosA=cos(A);
	double sin2A=sin(2.0*A);
	double cos2A=cos(2.0*A);
	double sinZ=sin(Z);
	double cosZ=cos(Z);
	double phi=q[4]*PI180;
	double st_lt=q[7]*3.4e-6; // deg/degC, see 2017/1 tech rep. for details
	
	// Model6a (same as 4g)
	MA+=p[0]+( p[1]*sinA-p[2]*cosA+p[3] )/tan(Z) + p[4]/sinZ;
	MA+=p[5]*sin2A + p[6]*cos2A + p[7]*sin(3.0*A)*cosZ + p[8]*cos(A/4)*sinZ;
	
	MZ+=q[0]+q[1]*cosA+q[2]*sinA+q[3]*sin(Z-phi)+q[5]*sin2A+q[6]*cos2A + st_lt*Tmeteo;
	
	*dAZ=MA;
	*dZD=MZ;
}


/***************************************************************************************/
/*!
	\brief implements Model6r 
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	(see 2017/1 tech rep. for details)
	\date Jun 22, 2017, 1:38:30 PM
	\author Bartosz Lew
*/
void Model6ar_DR2017(double AZ, double ZD, double *dAZ, double *dZD, double Tmeteo) {
	static const double deltaA_A_phi_T[51][3] = {
	{-4.4193510922E-05,0.0,0.0},
	{2.0975624936E-04,2.3739467890E+00,3.5995884719E+02},
	{2.1751766083E-05,1.5495360203E+00,1.7997942359E+02},
	{1.0854141784E-03,3.8681293724E-01,1.1998628240E+02},
	{2.8992718338E-04,-5.6802627748E-01,8.9989711797E+01},
	{4.4884800373E-04,-2.3100821930E+00,7.1991769438E+01},
	{1.3836783637E-03,2.1806533044E+00,5.9993141198E+01},
	{2.2485302979E-04,-2.0194550302E+00,5.1422692456E+01},
	{2.9168988970E-04,-1.3217108418E+00,4.4994855899E+01},
	{4.9852469088E-04,-2.1598779336E+00,3.9995427466E+01},
	{1.0794846894E-03,-1.7852151662E+00,3.5995884719E+01},
	{2.5641483632E-04,7.0085045918E-01,3.2723531563E+01},
	{1.8949331692E-04,1.8596130245E+00,2.9996570599E+01},
	{2.1418246633E-04,-1.7356228089E+00,2.7689142092E+01},
	{3.7778697518E-03,-2.3768899200E+00,2.5711346228E+01},
	{2.4361396666E-04,-2.5090394999E+00,2.3997256479E+01},
	{5.8319735260E-05,5.3737954663E-01,2.2497427949E+01},
	{1.0606995216E-04,1.1037184595E+00,2.1174049835E+01},
	{2.2855964603E-04,9.3851112570E-01,1.9997713733E+01},
	{7.0998829811E-05,2.4093838268E+00,1.8945202484E+01},
	{6.6321461852E-05,-2.9284123918E+00,1.7997942359E+01},
	{7.2439121418E-05,9.6744771972E-01,1.7140897485E+01},
	{1.3421645214E-04,-1.2826209016E+00,1.6361765781E+01},
	{7.1669917803E-05,-1.2539991715E+00,1.5650384660E+01},
	{1.3419292877E-04,2.6903416019E+00,1.4998285300E+01},
	{7.2135984931E-05,1.7750133268E+00,1.4398353888E+01},
	{1.4133007721E-04,1.2725667229E+00,1.3844571046E+01},
	{1.3074342218E-04,-9.2603046603E-01,1.3331809155E+01},
	{1.8804772953E-04,1.0562022024E+00,1.2855673114E+01},
	{1.4432651780E-04,-2.1050951511E+00,1.2412374041E+01},
	{1.5301098036E-04,4.8135220440E-01,1.1998628240E+01},
	{6.4849844692E-05,-2.3111608408E+00,1.1611575716E+01},
	{9.8254637506E-05,-2.0511060374E+00,1.1248713975E+01},
	{1.4873687336E-04,2.2041035384E+00,1.0907843854E+01},
	{2.0465524855E-04,2.1426705020E+00,1.0587024917E+01},
	{5.2404953408E-05,1.6013987175E-01,1.0284538491E+01},
	{7.3192280446E-05,-5.3187502061E-01,9.9988568664E+00},
	{1.1610973398E-04,2.5834079595E-03,9.7286174916E+00},
	{1.7855434004E-04,-1.2104279212E+00,9.4726012418E+00},
	{1.3093863735E-04,-6.2825321514E-01,9.2297140305E+00},
	{7.8778521918E-05,1.2796068420E+00,8.9989711797E+00},
	{2.2924306996E-04,-2.4347101338E+00,8.7794840778E+00},
	{1.5743712910E-03,3.1154526452E+00,8.5704487426E+00},
	{1.4934451616E-04,2.8588683103E+00,8.3711359812E+00},
	{7.7648989107E-05,-5.3312799189E-01,8.1808828907E+00},
	{1.4690714922E-04,6.9852473964E-01,7.9990854931E+00},
	{1.1152735371E-04,-1.2362716023E+00,7.8251923302E+00},
	{1.4563583864E-04,7.6110026708E-01,7.6586988764E+00},
	{3.2630776628E-05,2.7724211399E+00,7.4991426498E+00},
	{5.1016035439E-05,-1.0547400732E+00,7.3460989222E+00},
	{7.7384865383E-05,2.1933886797E+00,7.1991769438E+00}};

	static const double deltaZ_A_phi_T[51][3] = {
	{-2.0909256301E-06,0.0,0.0},
	{8.0452407247E-04,1.9121875725E+00,3.5984000000E+02},
	{4.1463723223E-04,-8.4166067330E-01,1.7992000000E+02},
	{3.9997244504E-03,1.1526356413E+00,1.1994666667E+02},
	{7.4605803303E-04,-4.8484466393E-01,8.9960000000E+01},
	{1.2891286257E-03,-6.1409989891E-01,7.1968000000E+01},
	{8.3728542684E-04,1.5981762051E+00,5.9973333333E+01},
	{5.2804849596E-04,2.5874104020E+00,5.1405714286E+01},
	{1.9511950535E-04,4.1685576907E-01,4.4980000000E+01},
	{5.6870079464E-04,-9.7990230793E-01,3.9982222222E+01},
	{1.0772761912E-03,-1.7817927339E+00,3.5984000000E+01},
	{4.6695635229E-04,-8.9493339976E-01,3.2712727273E+01},
	{2.0772602215E-04,-2.3397482069E+00,2.9986666667E+01},
	{2.4280117040E-05,1.3472521231E+00,2.7680000000E+01},
	{2.2854243381E-03,-2.4399206125E+00,2.5702857143E+01},
	{7.1086610829E-05,1.6317507601E+00,2.3989333333E+01},
	{9.7730058661E-05,-6.7772194791E-02,2.2490000000E+01},
	{9.4113607159E-05,-2.3070838808E+00,2.1167058824E+01},
	{2.6327972273E-04,-2.3709450809E+00,1.9991111111E+01},
	{1.8356945408E-04,-1.7903389230E+00,1.8938947368E+01},
	{2.6183945905E-04,3.1167668772E+00,1.7992000000E+01},
	{1.0468907683E-04,-2.9771743726E+00,1.7135238095E+01},
	{1.9353413819E-04,1.8404871659E+00,1.6356363636E+01},
	{3.5042690636E-05,2.2022397122E+00,1.5645217391E+01},
	{1.8551387754E-04,6.6659725680E-01,1.4993333333E+01},
	{1.4481335021E-04,-2.2128585238E+00,1.4393600000E+01},
	{1.7942386200E-04,-1.5721928595E+00,1.3840000000E+01},
	{1.3391580940E-04,2.2503893843E+00,1.3327407407E+01},
	{2.0966209384E-04,2.9549250645E+00,1.2851428571E+01},
	{1.0442523919E-04,2.4342560951E+00,1.2408275862E+01},
	{1.9249136747E-04,2.1102181561E+00,1.1994666667E+01},
	{2.0712222310E-04,-1.1463794230E+00,1.1607741935E+01},
	{1.7265015359E-04,-3.5849644240E-01,1.1245000000E+01},
	{1.3363967362E-04,1.9824858508E-01,1.0904242424E+01},
	{1.2752796026E-04,1.7827899611E+00,1.0583529412E+01},
	{9.6964494004E-05,2.4757027425E+00,1.0281142857E+01},
	{6.1383996581E-05,-1.4179002054E+00,9.9955555556E+00},
	{2.0935022711E-04,-7.0095164664E-01,9.7254054054E+00},
	{1.6857588910E-04,-2.1253028536E-01,9.4694736842E+00},
	{4.5236991056E-06,-6.1641383706E-01,9.2266666667E+00},
	{2.2027483602E-04,-1.1096768003E+00,8.9960000000E+00},
	{1.2665595930E-04,4.1733550831E-01,8.7765853659E+00},
	{9.2160513055E-04,2.9919983381E+00,8.5676190476E+00},
	{1.6677608461E-04,-2.6791178297E+00,8.3683720930E+00},
	{1.1268019620E-04,2.0030724026E+00,8.1781818182E+00},
	{2.0899575575E-04,-2.9393082826E+00,7.9964444444E+00},
	{1.5170437568E-04,-2.2376710352E+00,7.8226086957E+00},
	{1.7000219780E-04,-1.8803770653E+00,7.6561702128E+00},
	{9.6047496464E-05,1.2973590132E+00,7.4966666667E+00},
	{8.1961123891E-05,1.6073261308E+00,7.3436734694E+00},
	{1.3990622468E-04,-8.4979239250E-01,7.1968000000E+00}};
	
	
	double twoPI=M_PI*2.0;

	double MA=0,MZ=0;

	int i;

	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	// Model6
	Model6a_DR2017(AZ,ZD,&MA,&MZ,Tmeteo);
	
	// Model6r
	MA+=deltaA_A_phi_T[0][0];
	MZ+=deltaZ_A_phi_T[0][0];
	for (i = 1; i < 51; i++) {
		MA+=deltaA_A_phi_T[i][0]*sin(-twoPI/deltaA_A_phi_T[i][2]*(AZ-180)+deltaA_A_phi_T[i][1]);
		MZ+=deltaZ_A_phi_T[i][0]*sin(-twoPI/deltaZ_A_phi_T[i][2]*(AZ-180)+deltaZ_A_phi_T[i][1]);
	}
	
	*dAZ=MA;
	*dZD=MZ;
}
/* ******************************************************************************************** */
void aging_linear_driftZD(double AZ, double ZD, double *dAZ, double *dZD,double JD) {
	// convert JD to years since 2015-01-01
	double JDref = 2457023.5; // JulianDay(2015,1,1);
	double t=(JD-JDref)/365.25;
	
	double A=19.15278693; // mdeg/yr
	double B=-32.06389286; // mdeg
	
	double dzd=A*t+B;

	*dZD+=dzd/1000; // convert to deg
	
}
/* ******************************************************************************************** */
void aging_linearT_driftZD(double AZ, double ZD, double *dAZ, double *dZD, double JD, double Tmeteo) {
	// convert JD to years since 2015-01-01
	double JDref = 2457023.5; // JulianDay(2015,1,1);
	double t=(JD-JDref)/365.25;
	
	double A=18.7101907; // mdeg/yr
	double B=-34.4520555; // mdeg
	double C=0.283225502; // mdeg/K
	
	double dzd=A*t+B+C*Tmeteo;

	*dZD+=dzd/1000; // convert to deg
}
