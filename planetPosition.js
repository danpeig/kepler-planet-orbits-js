/*
* Keplerian Elements for Approximate Positions of the Major Planets
* Based on the paper written by E. M. Standish, Solar System Dynamics Group JPL/Caltech
* https://ssd.jpl.nasa.gov/?planet_pos
* 
* Programmed by Daniel BP (www.danbp.org), 02/06/2019
* Updated on 23/11/2020
*/

function computeEclipticCoordinates(keplerianElementsTable, planetIndex, UNIXMilliTime){
	let deg2rad = Math.PI/180;
	let Teph = UNIXMilliTime/1000/86400 + 2440587.5; //Julian Ephemeris Date
	let T = (Teph-2451545)/36525; //Centuries past J2000.0
	
	//Calculations of the Keplerian elements
	let a0 = keplerianElementsTable[planetIndex][1]; //semi-major axis (initial)
	let at = keplerianElementsTable[planetIndex][7]; //semi-major axis (century variation)
	let a = a0+at*T; //semi-major axis
	let e0 = keplerianElementsTable[planetIndex][2]; //eccentricity (initial)
	let et = keplerianElementsTable[planetIndex][8]; //eccentricity (century variation)
	let e = e0+et*T; //eccentricity 
	let I0 = keplerianElementsTable[planetIndex][3]; //inclination (initial)
	let It = keplerianElementsTable[planetIndex][9]; //inclination (century variation)
	let I = I0 + It*T; //inclination
	let L0 = keplerianElementsTable[planetIndex][4]; //mean longitude (initial)
	let Lt = keplerianElementsTable[planetIndex][10]; //mean longitude (century variation)
	let L = L0+Lt*T; //mean longitude
	let W0 = keplerianElementsTable[planetIndex][5]; //longitude of perihelion (initial)
	let Wt = keplerianElementsTable[planetIndex][11]; //longitude of perihelion (century variation)
	let W = W0+Wt*T; //longitude of perihelion
	let O0 = keplerianElementsTable[planetIndex][6]; //longitude of ascending node (initial)
	let Ot = keplerianElementsTable[planetIndex][12]; //longitude of ascending node (century variation)
	let O = O0+Ot*T; //longitude of ascending node
	
	//Auxiliary parameters for the 3000 BC - 3000 AD period
	let b = keplerianElementsTable[planetIndex][13];
	let c = keplerianElementsTable[planetIndex][14];
	let s = keplerianElementsTable[planetIndex][15];
	let f = keplerianElementsTable[planetIndex][16];
	
	//Argument of perihelion w
	let w = W-O;
	
	//Mean anomaly M
	let M = L-W+b*Math.pow(T,2)+c*Math.cos(deg2rad*f*T)+s*Math.sin(deg2rad*f*T);
	if(M < -180 || M > 180) M = M % 360; //Modulus M to be between -180 and 180
	
	//Solve the Kepler Equation for the eccentric anomaly E
	let E = solveKeplerEquation(M,e); 
	
	//Heliocentric coordinates
	let x1 = a*(Math.cos(deg2rad*E)-e);
	let y1 = a*Math.sqrt(1-(e*e))*Math.sin(deg2rad*E);
	let z1 = 0;
	
	//Auxiliary trigonometric calculations
	let cosw = Math.cos(deg2rad*w);
	let sinw = Math.sin(deg2rad*w);
	let cosO = Math.cos(deg2rad*O); 
	let sinO = Math.sin(deg2rad*O);
	let cosI = Math.cos(deg2rad*I);
	let sinI = Math.sin(deg2rad*I);
	
	//J2000 Ecliptic Plane Coordinates
	let Xecl = (cosw*cosO - sinw*sinO*cosI)*x1 + (-sinw*cosO - cosw*sinO*cosI)*y1;
	let Yecl = (cosw*sinO + sinw*cosO*cosI)*x1 + (-sinw*sinO + cosw*cosO*cosI)*y1;
	let Zecl = (sinw*sinI)*x1 + (cosw*sinI)*y1;
	
	//Result array
	let EclipticCoordinates = [Xecl, Yecl, Zecl];
	
	return EclipticCoordinates;
}

function convertToICRF(EclipticCoordinates){
	let deg2rad = Math.PI/180;
	let E = 23.43928;
	let cosE = Math.cos(deg2rad*E);
	let sinE = Math.sin(deg2rad*E);
	let Xeq = EclipticCoordinates[0];
	let Yeq = cosE*EclipticCoordinates[1]-sinE*EclipticCoordinates[2];
	let Zeq = sinE*EclipticCoordinates[1]+cosE*EclipticCoordinates[2];
	
	//Result array
	let ICRFCoordinates = [Xeq, Yeq, Zeq];

	return ICRFCoordinates;
}

function solveKeplerEquation(M,e) {
	let e1 = 180/Math.PI*e;
	let deg2rad = Math.PI/180;
	let E0 = M+e1*Math.sin(M*deg2rad);
	let tolerance = 10e-6;
	let dE;
	let dM;
	let E1 = E0;
	let E2 = 0;
	let difE1E2 = 1;
	
	for (let i=0; (i<1000) && (difE1E2 > tolerance); i++){
		dM = M - (E1 - e1*Math.sin(deg2rad*E1));
		dE = dM/(1-e*Math.cos(deg2rad*E1));
		E2=E1+dE;
		difE1E2 = Math.abs(E1-E2);
	}
	return E2;
}


//Keplerian elements and their rates, with respect to the mean ecliptic
//and equinox of J2000, valid for the time-interval 1800 AD - 2050 AD.
//
//               a0              e0               I0                L0                 W0              O0                 at               et               It                 Lt               Wt               Ot             b           c           s         f
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
let keplerianElements2050 = [
["Mercury",   0.38709927,      0.20563593,      7.00497902,      252.25032350,     77.45779628,     48.33076593,        0.00000037 ,     0.00001906 ,     -0.00594749,    149472.67411175,      0.16047689,     -0.12534081,     0,          0,         0,         0],
["Venus",     0.72333566,      0.00677672,      3.39467605,      181.97909950,    131.60246718,     76.67984255,        0.00000390 ,     -0.00004107,     -0.00078890,    58517.81538729 ,      0.00268329,     -0.27769418,     0,          0,         0,         0],
["EM Bary",   1.00000261,      0.01671123,     -0.00001531,      100.46457166,    102.93768193,      0.0       ,   	    0.00000562 ,     -0.00004392,     -0.01294668,    35999.37244981 ,      0.32327364,     0.0        ,     0,          0,         0,         0],
["Mars",	  1.52371034,	   0.09339410,	    1.84969142, 	  -4.55343205, 	  -23.94362959, 	49.55953891,        0.00001847 , 	-0.00007882 ,	  -0.00813131,	  19140.30268499 ,	    0.44441088, 	-0.29257343, 	 0, 	     0, 	    0, 	   	   0],
["Jupiter",   5.20288700,      0.04838624,      1.30439695,       34.39644051,     14.72847983,    100.47390909,        -0.00011607,    -0.00013253 ,     -0.00183714,     3034.74612775 ,      0.21252668,     0.20469106 ,     0,          0,         0,         0],
["Saturn",    9.53667594,      0.05386179,      2.48599187,       49.95424423,     92.59887831,    113.66242448,        -0.00125060,    -0.00050991 ,      0.00193609,     1222.49362201 ,     -0.41897216,     -0.28867794,     0,          0,         0,         0],
["Uranus",   19.18916464,      0.04725744,      0.77263783,      313.23810451,    170.95427630,     74.01692503,        -0.00196176,    -0.00004397 ,     -0.00242939,     428.48202785  ,      0.40805281,      0.04240589,     0,          0,         0,         0],
["Neptune",  30.06992276,      0.00859048,      1.77004347,      -55.12002969,     44.96476227,    131.78422574,         0.00026291,     0.00005105 ,      0.00035372,      218.45945325 ,     -0.32241464,     -0.00508664,     0,          0,         0,         0],
["Pluto",    39.48211675,      0.24882730,     17.14001206,      238.92903833,    224.06891629,    110.30393684,        -0.00031596,     0.00005170 ,      0.00004818,      145.20780515 ,     -0.04062942,     -0.01183482,     0,          0,         0,         0]
];



//Keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000,
//valid for the time-interval 3000 BC -- 3000 AD.
//
//
//               a0              e0               I0                L0                 W0              O0                 at               et               It                 Lt               Wt               Ot                  b                c             s                f
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
let keplerianElements3000 = [
["Mercury",   0.38709843,      0.20563661,      7.00559432,      252.25166724,     77.45771895,     48.33961819,        0.00000000,      0.00002123,     -0.00590158,   149472.67486623,      0.15940013,     -0.12214182,     0,               0,              0,             0],
["Venus",     0.72332102,      0.00676399,      3.39777545,      181.97970850,    131.76755713,     76.67261496,       -0.00000026,     -0.00005107,      0.00043494,    58517.81560260,      0.05679648,     -0.27274174,     0,               0,              0,             0],
["EM Bary",   1.00000018,      0.01673163,     -0.00054346,      100.46691572,    102.93005885,     -5.11260389,       -0.00000003,     -0.00003661,     -0.01337178,    35999.37306329,      0.31795260,     -0.24123856,     0,               0,              0,             0],
["Mars",      1.52371243,      0.09336511,      1.85181869,       -4.56813164,    -23.91744784,     49.71320984,        0.00000097,      0.00009149,     -0.00724757,    19140.29934243,      0.45223625,     -0.26852431,     0,               0,              0,             0],
["Jupiter",   5.20248019,      0.04853590,      1.29861416,       34.33479152,     14.27495244,    100.29282654,       -0.00002864,      0.00018026,     -0.00322699,     3034.90371757,      0.18199196,      0.13024619,     -0.00012452,     0.06064060,     -0.35635438,   38.35125000],
["Saturn",    9.54149883,      0.05550825,      2.49424102,       50.07571329,     92.86136063,    113.63998702,       -0.00003065,     -0.00032044,      0.00451969,     1222.11494724,      0.54179478,     -0.25015002,     0.00025899,      -0.13434469,    0.87320147,    38.35125000],
["Uranus",   19.18797948,      0.04685740,      0.77298127,      314.20276625,    172.43404441,     73.96250215,       -0.00020455,     -0.00001550,     -0.00180155,      428.49512595,      0.09266985,      0.05739699,     0.00058331,      -0.97731848,    0.17689245,    7.67025000],
["Neptune",  30.06952752,      0.00895439,      1.77005520,      304.22289287,     46.68158724,    131.78635853,        0.00006447,      0.00000818,      0.00022400,      218.46515314,      0.01009938,     -0.00606302,     -0.00041348,     0.68346318,     -0.10162547,   7.67025000],
["Pluto",    39.48686035,      0.24885238,     17.14104260,      238.96535011,    224.09702598,    110.30167986,        0.00449751,      0.00006016,      0.00000501,      145.18042903,     -0.00968827,     -0.00809981,     -0.01262724,     0,              0,             0]
];
