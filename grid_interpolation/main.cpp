#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "include/armadillo"

//WGS84 constantes
const double a = 6378137;
const double b = 6.356752314245179e+06;
const double e2 = 0.006694379990141;

//GRS80 constantes
//const double a = 6378137;
//const double b = 6.356752314140356e+06;
//const double e2 = 0.006694380022901;

const double pi = 3.141592653589793;
const double torad = pi/180;

// namespaces to use the operators overloading
using namespace arma;
using namespace std;

void ecef2lla(double x, double y, double z, double* lat, double* lon, double* h)
{
	double ep,p,th,N,k,t;
	double e = sqrt(e2);

	ep = sqrt((pow(a,2)-pow(b,2))/pow(b,2));
	p = sqrt(pow(x,2)+pow(y,2));
	th = atan2(a*z,b*p);
	*lon = atan2(y,x);
	*lat = atan2((z+pow(ep,2)*b*pow(sin(th),3)),(p-e2*a*pow(cos(th),3)));
	N = a/sqrt(1-e2*pow(sin(*lat),2));
	*h = p/cos(*lat)-N;

	//return lon in range [0,2*pi)
	*lon = fmod(*lon,2*pi)/torad;
	*lat = *lat/torad;
	//DDG: extraído de la función matlab. No necesitamos esta corrección
	//correct for numerical instability in altitude near exact poles:
	//(after this correction, error is about 2 millimeters, which is about
	//the same as the numerical precision of the overall function)
	//k=abs(x)<1 & abs(y)<1;
	//alt(k) = abs(z(k))-b;
}

void lla2ecef(double lat, double lon, double alt, double* x, double* y, double* z)
{
	double N;	
	//intermediate calculation
	//(prime vertical radius of curvature)
	N = a/sqrt(1-e2*pow(sin(lat*torad),2));

	//results:
	*x = (N+alt)*cos(lat*torad)*cos(lon*torad);
	*y = (N+alt)*cos(lat*torad)*sin(lon*torad);
	*z = ((1-e2)*N + alt)*sin(lat*torad);
}

void lg2ct(double north,double east,double up,double lat,double lon,double* dx, double* dy, double* dz)
{
	mat R = zeros(3,3);
	mat C = zeros(3,1);

	R(0,0)=-sin(lat*torad)*cos(lon*torad);
	R(1,0)=-sin(lat*torad)*sin(lon*torad);
	R(2,0)=cos(lat*torad);
	R(0,1)=-sin(lon*torad);
	R(1,1)=cos(lon*torad);
	R(2,1)=0;
	R(0,2)=cos(lat*torad)*cos(lon*torad);
	R(1,2)=cos(lat*torad)*sin(lon*torad);
	R(2,2)=sin(lat*torad);

	C(0,0)=north;
	C(1,0)=east;
	C(2,0)=up;	

	mat L=R*C;

	*dx = L(0,0);
	*dy = L(1,0);
	*dz = L(2,0);
}

string todms(double deg)
{
	double d,id,rm,s,im;
	string dms="";
	ostringstream ss;

	d=abs(deg);
	id=floor(d);
	rm=(d-id)*60;
	im=floor(rm);
	s=(rm-im)*60;

	ss << id << " " << im << " " << s;
	dms = ss.str();
	return dms;
}

int isecef(double* result)
{
	double R = sqrt(pow(result[0],2) + pow(result[1],2) + pow(result[2],2));

	if (R>6e+06)
		return 1;
	else
		return 0;
}

void print_output(double* result)
{
	//detecto si lat lon h son en realidad XYZ
	int isxyz = isecef(result);

	// salida a stdout
	for (int i=0;i<=19;i++)
	{
		switch (i)
		{
			case 0:
			case 1:
				if (isxyz == 1)
					printf("%+13.4f ", result[i]);  // ECEF
				else
					printf("%+12.10f ", result[i]); // deg
				break;
			case 2:
			case 7:
				if (isxyz == 1)
					printf("%+13.4f ", result[i]);  // ECEF
				else
					printf("%+9.4f ", result[i]);  // alturas elipsiodales
				break;
			case 3:
			case 4:
				printf("%+9.4f ", result[i]);  // epochas
				break;
			case 5:
			case 6:
				if (isxyz == 1)
					printf("%+13.4f ", result[i]);  // ECEF
				else
					printf("%+12.10f ", result[i]); // deg
				break;
			default:
				printf("%+7.4f ", result[i]);  // todo el resto
				break;
		}			
	}	
	printf("\n");
	
}

void lle2utm(double lat,double lon,double lcm, double *NE)
{
	// utilizo WGS84
	double f=1-sqrt(1-e2);

	lat = lat*torad;
	lon = lon*torad;
	lcm = lcm*torad;

	double ko=0.9996;    // Scale factor
	double No=1e7;       // False northing (north)
	double Eo=0;         // False easting

	double lam=lon-lcm;
	  
	double RN=a/pow(1-e2*pow(sin(lat),2),0.5);
	double RM=a*(1-e2)/pow(1-e2*pow(sin(lat),2),1.5);
	double h2=e2*pow(cos(lat),2)/(1-e2);
	double t=tan(lat);
	double n=f/(2-f);

	//----- Helmert (1880) expansion & simplification of Bessel series (faster)
	double A0=1+pow(n,2)/4+pow(n,4)/64;
	double A2=3/2*(n-pow(n,3)/8);
	double A4=15/16*(pow(n,2)-pow(n,4)/4);
	double A6=35/48*pow(n,3);
	double A8=315/512*pow(n,4);
	double S=a/(1+n)*(A0*lat-A2*sin(2*lat)+A4*sin(4*lat)-A6*sin(6*lat)+A8*sin(8*lat));


	double E1=lam*cos(lat);
	double E2=pow(lam,3)*pow(cos(lat),3)/6*(1-pow(t,2)+h2);
	double E3=pow(lam,5)*pow(cos(lat),5)/120*(5-18*pow(t,2)+pow(t,4)+14*h2-58*pow(t,2)*h2+13*pow(h2,2)+4*pow(h2,3)-64*pow(t,2)*pow(h2,2)-24*pow(t,2)*pow(h2,3));
	double E4=pow(lam,7)*pow(cos(lat),7)/5040*(61-479*pow(t,2)+179*pow(t,4)-pow(t,6));
	double E=Eo+ko*RN*(E1+E2+E3+E4);

	double N1=S/RN;
	double N2=pow(lam,2)/2*sin(lat)*cos(lat);
	double N3=pow(lam,4)/24*sin(lat)*pow(cos(lat),3)*(5-pow(t,2)+9*h2+4*pow(h2,2));
	double N4=pow(lam,6)/720*sin(lat)*pow(cos(lat),5)*(61-58*pow(t,2)+pow(t,4)+ 270*h2-330*pow(t,2)*h2+445*pow(h2,2)+324*pow(h2,3)-680*pow(t,2)*pow(h2,2)+ 88*pow(h2,4)-600*pow(t,2)*pow(h2,3)-192*pow(t,2)*pow(h2,4));
	double N5=pow(lam,8)/40320*sin(lat)*pow(cos(lat),7)*(1385-311*pow(t,2)+543*pow(t,4)-pow(t,6));
	double N=No+ko*RN*(N1+N2+N3+N4+N5);

	NE[0] = N/1000 - 6200;
	NE[1] = E/1000;
}

void get_variable(double lat, double lon, string archivo, double *val)
{
	mat vel_ar;
	double R = 6371; // radio terrestre

	// cargar los archivos del modelo
	vel_ar.load(archivo,raw_ascii);

	mat vel_ar_rad = vel_ar*torad;

	// busco los puntos mas cercanos a los parametros ingresados
	mat pos = ones(vel_ar.n_rows,1);

	pos = join_horiz(pos*lat*torad, pos*lon*torad);

	mat dist = zeros(vel_ar.n_rows,1);
	
	// distancia sobre circulo mayor utilizando la formula estable de haverside
	mat d_phi = (vel_ar_rad.col(0) - pos.col(0))/2;
	mat d_lam = (vel_ar_rad.col(1) - pos.col(1))/2;

	dist = 2*asin(sqrt(square(sin(d_phi)) + cos(vel_ar_rad.col(0))%cos(pos.col(0))%square(sin(d_lam))))*R;

	// busco los indices ordenados
	uvec index = sort_index(dist);
	// obtengo un vector segun los indices ordenados por distancia
	mat inter = join_horiz(dist(index,zeros<uvec>(1,1)), vel_ar(index,linspace<uvec>(0,vel_ar.n_cols-1,vel_ar.n_cols)));	

	// me quedo con 4 elementos
	int points = 4;
	inter = inter.rows(0,points-1);

	// temporary vector to save the plane coordinates
	double tNE[2];
	double XY[2];

	mat NE = zeros(points,2);

	for(int i = 0;i<points;i++) 
	{
		lle2utm(inter(i,1),inter(i,2),lon,tNE);
		NE(i,0) = tNE[0];
		NE(i,1) = tNE[1];
	}

	// mmcc
	mat A = join_horiz(ones<mat>(points,1), join_horiz(NE.col(0),NE.col(1)));
	
	// convert the interpolation lat lon to plane coordinates
	lle2utm(lat,lon,lon,tNE);
	mat a = join_horiz(ones<mat>(1,1), join_horiz(ones<mat>(1,1)*tNE[0],ones<mat>(1,1)*tNE[1]));

	mat L1 = inter.col(3);
	mat L2 = inter.col(4);

	mat X1 = solve(A.t()*A, A.t()*L1);
	mat X2 = solve(A.t()*A, A.t()*L2);
	
	val[0] = as_scalar(a*X1);
	val[1] = as_scalar(a*X2);

	// debug: muestra la matriz de interpolación
	//cout << inter;	

	// interpolo por el método de inverso a la distancia (exponente = 1)
	//double i = 1;
	//val[0] = sum(pow(1/inter.col(0),i)%inter.col(3))/sum(pow(1/inter.col(0),i));
	//val[1] = sum(pow(1/inter.col(0),i)%inter.col(4))/sum(pow(1/inter.col(0),i));

	// debug
	//cout << "Latitud longitud ingresada: " << lat << " " << lon << endl;
	//cout << "Matriz de valores: " << endl << inter;
	//cout << "Valor interpolado N/E: " << val[0] << " " << val[1] << endl;

	return;
}

void apply_vel_ar(double lat, double lon, double h, double epoch, double target_epoch, double *result)
{
	double norte,este,norte_pos,este_pos,x,y,z;
	double vel_ar_lin[2];
	double vel_ar_log[2];
	double vel_ar_cos[2];

	int isxyz = 0;

	// epoca del terremoto de Maule
	double te = 2.0101589e+03;

	double d1 = epoch - te;
	double d2 = target_epoch - te;
	double d3 = target_epoch - epoch;

	//verifico si lat lon h son en realidad XYZ
	double R;
	R = sqrt(pow(lat,2) + pow(lon,2) + pow(h,2));

	if (R>6e+06)
	{
		//ECEF, convertir a lat lon h (guardo las coordenadas originales)
		x = lat;
		y = lon;
		z = h;
		
		ecef2lla(lat,lon,h,&lat,&lon,&h);
		
		clog << "Coordenadas ECEF detectadas > CONVERSION: " << lat << " " << lon << " " << h << endl;
		isxyz = 1;
	}	

	// determino la velocidad lineal
	get_variable(lat,lon,"vel-ar-lin.dat",vel_ar_lin);

	// determino los parametros logaritmicos
    if (lat > -30 || lat < -42)
    {
        vel_ar_log[0] = 0;
        vel_ar_log[1] = 0;
    }
    else
        get_variable(lat,lon,"vel-ar-log.dat",vel_ar_log);

	// determino el salto co-sismico
	get_variable(lat,lon,"vel-ar-cos.dat",vel_ar_cos);

	// velocidad lineal
	// entre épocas anteriores al sismo, sin pasar por el sismo solo aplicamos la componente lineal
	// esta componente está siempre presente así que se aplica siempre utilizando el signo de d3
	norte = vel_ar_lin[0]*d3;
	este  = vel_ar_lin[1]*d3;

	// hacia adelante, pasando por el sismo
	if (d1 < 0 && d2 > 0)
	{
		// cout << "hacia adelante " << d3 << endl;
		// cout << "log " << vel_ar_log[0] << " " << vel_ar_log[1] << endl;

		// estoy yendo "hacia adelante", aplico el salto co-sísmico positivo
		norte = norte + vel_ar_cos[0];
		este  =  este + vel_ar_cos[1];

		// agrego la componente logaritmica
		norte_pos = vel_ar_log[0]*log10(1+(d2/0.5));
		este_pos  = vel_ar_log[1]*log10(1+(d2/0.5));
		norte = norte + norte_pos;
		este  =  este + este_pos;
	}

	// hacia atras, pasando por el sismo
	if (d1 > 0 && d2 < 0)
	{
		// cout << "hacia atras" << endl;
		// estoy yendo hacia "atrás", aplico el salto co-sísmico negativo
		norte = norte - vel_ar_cos[0];
		este  =  este - vel_ar_cos[1];

		// agrego la componente logaritmica
		norte_pos = vel_ar_log[0]*log10(1+(d1/0.5));
		este_pos  = vel_ar_log[1]*log10(1+(d1/0.5));
		norte = norte - norte_pos;
		este  =  este - este_pos;
	}

	// entre épocas posteriores al sismo, sin pasar por el sismo (lineal + log)
	if (d1 > 0 && d2 > 0)
	{
		// agrego la componente logaritmica. La dirección es tenida en cuenta por d1 y d2.
		norte_pos = vel_ar_log[0]*log10(1+(d2/0.5)) - vel_ar_log[0]*log10(1+(d1/0.5));
		este_pos  = vel_ar_log[1]*log10(1+(d2/0.5)) - vel_ar_log[1]*log10(1+(d1/0.5));
		norte = norte + norte_pos;
		este  =  este + este_pos;
	}
	
	// convierto los valores de norte/este a lat lon
	// utilizo WGS84
	double v,r,dlat,dlon;

	v=a/sqrt(1-e2*pow(sin(lat*torad),2));
	r=v*(1-e2)/(1-e2*pow(sin(lat*torad),2));

	// OJO! Falta altura elipsoidal
	dlat=(norte/(r+h))/torad;
	dlon=(este/cos(lat*torad)/(v+h))/torad;

	if (isxyz == 1)
	{
		result[0] = x;
		result[1] = y;
		result[2] = z;
		// convierto la nueva latitud y longitud a ECEF
		lla2ecef(lat+dlat,lon+dlon,h,&x,&y,&z);
		result[5] = x;
		result[6] = y;
		result[7] = z;

		// Convierto las diferencias NEU a ECEF. Up = 0 porque no hay componente UP en esta versión de Vel-Ar
		lg2ct(norte,este,0,lat+dlat,lon+dlon,&x,&y,&z);
		result[8] = x;
		result[9] = y;
		result[10] = z;

		// lo mismo para el salto co-sísmico
		lg2ct(vel_ar_cos[0],vel_ar_cos[1],0,lat+dlat,lon+dlon,&x,&y,&z);
		result[11] = x;
		result[12] = y;
		result[13] = z;

		// componente logarítmica
		lg2ct(vel_ar_log[0],vel_ar_log[1],0,lat+dlat,lon+dlon,&x,&y,&z);
		result[14] = x;
		result[15] = y;
		result[16] = z;

		// componente lineal
		lg2ct(vel_ar_lin[0],vel_ar_lin[1],0,lat+dlat,lon+dlon,&x,&y,&z);
		result[17] = x;
		result[18] = y;
		result[19] = z;

		// DEBUG: 
		// clog << "Norte y este en ECEF > CONVERSION: " << x << " " << y << " " << z << endl;
		// clog << result[5]-result[0] << " " << result[6]-result[1] << " " << result[7]-result[2]  << endl;
	}
	else
	{
		result[0] = lat;
		result[1] = lon;
		result[2] = h;
		result[5] = lat+dlat;
		result[6] = lon+dlon;
		result[7] = h;
		result[8] = norte;
		result[9] = este;
		result[10] = 0;
		result[11] = vel_ar_cos[0];
		result[12] = vel_ar_cos[1];
		result[13] = 0;
		result[14] = vel_ar_log[0];
		result[15] = vel_ar_log[1];
		result[16] = 0;
		result[17] = vel_ar_lin[0];
		result[18] = vel_ar_lin[1];
		result[19] = 0;
	}
	// estos parametros no cambian si es ECEF o lla.
	result[3] = epoch;
	result[4] = target_epoch;

	//clog << " Aplicando Vel-Ar >> Inc. en norte: " << setw(7) << norte << " Inc. en este : " << setw(7) << este << endl;
	fprintf(stderr," Aplicando Vel-Ar >> Inc. en norte: %+6.3f Inc. en este : %+6.3f\n", norte, este);

	return;
}

void puntoxpunto()
{
	double lat,lon,h,epoch,target_epoch,norte,este,norte_pos,este_pos;
	double result[20];

	cout << "  Ingrese latitud (N = +; S = -) o coordenada ECEF X: ";
	cin >> lat;
	cout << " Ingrese longitud (E = +; W = -) o coordenada ECEF Y: ";
	cin >> lon;
        cout << "    Ingrese la altura elipsoidal o coordenada ECEF Z: ";
        cin >> h;
	cout << "      Ingrese la época de observación (ej: 2006.632): ";
	cin >> epoch;
	cout << "          Ingrese la época de destino (ej: 2012.125): ";
	cin >> target_epoch;

	apply_vel_ar(lat,lon,h,epoch,target_epoch,result);
	
	int isxyz = isecef(result);
	
	if (isxyz == 1)
	{
		cout << "-------------------------------------------" << endl;
		printf(" ECEF XYZ (%9.4f): %+13.4f %+13.4f %+13.4f\n",result[3],result[0],result[1],result[2]);

		printf(" Nueva X (%9.4f): %+13.4f\n",result[4],result[5]);
		printf(" Nueva Y (%9.4f): %+13.4f\n",result[4],result[6]);
		printf(" Nueva Z (%9.4f): %+13.4f\n",result[4],result[7]);
	}
	else
	{
		cout << "-------------------------------------------" << endl;
		printf("  Latitud ingresada (%9.4f): %+12.10f (%s)\n",result[3],result[0],todms(result[0]).c_str());
		printf(" Longitud ingresada (%9.4f): %+12.10f (%s)\n",result[3],result[1],todms(result[1]).c_str());
		printf("      Nueva latitud (%9.4f): %+12.10f (%s)\n",result[4],result[5],todms(result[5]).c_str());
		printf("     Nueva longitud (%9.4f): %+12.10f (%s)\n",result[4],result[6],todms(result[6]).c_str());
	}

	clog << endl << " Terminado." << endl;

	return;
}

void list_file(string entrada)
{
	mat ment;
	// string entrada="posgar_test.txt", salida="out.txt";
	double result[20];

	// leo el archivo de entrada
	ment.load(entrada,raw_ascii);

	for(int i=0;i<ment.n_rows;i++)
	{
		apply_vel_ar(ment(i,0),ment(i,1),ment(i,2),ment(i,3),ment(i,4),result);

	        print_output(result);
	}
	clog << endl << " Terminado." << endl;
	
	return;
}

void main_menu()
{
	char myChar  = {0};
	string input = "";
	string entrada = "";

	while (true) {
		system("clear");
		clog << "         Instituto Geográfico Nacional Argentino " << endl;
		clog << "               The University of Memphis         " << endl;
		clog << "         Modelo Velocidades Argentinas (Vel-Ar)  " << endl;
		clog << "        -----------------------------------------" << endl;
		clog << "        Seleccione una de las siguientes opciones" << endl;
		clog << "        1) Aplicar Vel-Ar a un solo punto"         << endl;
		clog << "        2) Aplicar Vel-Ar a un archivo de puntos"  << endl;
		clog << "        s) Salir del programa Vel-Ar" << endl;
		clog << endl;
		clog << "        Tambien puede invocar el modelo Vel-Ar utilizando la linea de comandos:" << endl;
		clog << "          vel-ar [archivo_de_puntos] > [archivo de salida]" << endl;
		clog << "          vel-ar [lat] [lon] [h] [epoca_inicio] [epoca_fin]" << endl;
		clog << "        El archivo de puntos debe tener el siguiente formato:" << endl;
		clog << "          lat lon altura epoca_inicio epoca_fin" << endl;
		clog << "        En todos los casos, el formato de salida es:" << endl;
		clog << "          lat lon altura t_inicio t_fin lat+offset lon+offset altura offset_norte offset_este b_norte b_este a_norte a_este p_norte p_este" << endl;

		clog << "~$ ";
		getline(cin, input);

		if ((input.length() == 1) && (input.compare("1") == 0 || input.compare("2") == 0 || input.compare("s") == 0))
		{
			myChar = input[0];
			char option = input[0];
			switch (option)
			{
				case '1':
					system("clear");
					puntoxpunto();
					return;
				case '2':
					system("clear");
					cout << " Ingrese el archivo de entrada: ";
					cin >> entrada;
					list_file(entrada);
					return;
				case 's':
					return;
			}
		}else{
			cout << "Selección inválida! Por favor, seleccione 1, 2 o (s)alir" << endl;
		}
	}

	return;

}

int main(int argc, char *argv[])
{
	string arg;
        double result[20];
	double argd[5];
        int i=0;
        
	//system("clear");
	clog << "         Instituto Geográfico Nacional Argentino " << endl;
	clog << "               The University of Memphis         " << endl;
	clog << "         Modelo Velocidades Argentinas (Vel-Ar)  " << endl;
	clog << "        -----------------------------------------" << endl;

	switch (argc)
	{
		case 2:
			if (string(argv[1]) == "-i")
                                main_menu();
			else
			{
				list_file(string(argv[1]));
			}	
			break;
		case 6:
			apply_vel_ar(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), result);

			print_output(result);

			break;
		case 1:
			// este es el caso sin argumentos. Pruebo leer el stdin
			// leo el stdin para ver si pasaron argumentos desde allí
			for (string line; getline(cin, line);) {
				stringstream stream(line);
                                i = 0;
				while( getline(stream, arg, ' ') )
                                {
                                        argd[i] = atof(arg.c_str());
					//cout << argd[i] << endl; DEBUG
                                        i++;
                                }
                                //verificar si las coordenadas son lat lon h o XYZ
                                apply_vel_ar(argd[0], argd[1], argd[2], argd[3], argd[4], result);

			        print_output(result);
			}
			break;
		default:
			clog << "El numero de argumentos no es valido. Opciones:" << endl;
			clog << "Para convertir un archivo de coordenadas de una epoca a otra:" << endl;
			clog << "  vel-ar [archivo_de_puntos] > [archivo de salida]" << endl;
			clog << "Para utilizar la interfaz interactiva:" << endl;
			clog << "  vel-ar -i" << endl;
			clog << "Para convertir una sola coordenada:" << endl;
			clog << "  vel-ar [lat] [lon] [h] [epoca_inicio] [epoca_fin]" << endl;
			clog << "Alternativamente es posible ingresar coordenadas ECEF (X,Y,Z) en vez de lat lon h" << endl;
			clog << " en cuyo caso todas las salidas también serán en ECEF" << endl;
	}
	
	return 0;
}

