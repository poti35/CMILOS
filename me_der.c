#include "defines.h"
#include "lib.h"
#include <string.h>
#include "fftw.h"
#include "milosUtils.h"
#include "convolution.h"


int funcionComponentFor(int n_pi,PRECISION iwlines,int numl,PRECISION *wex,float *nuxB,PRECISION *dfi,PRECISION *dshi,PRECISION LD,PRECISION A,int desp);

void Resetear_Valores_Intermedios(int nlambda);



/*
	E00	int eta0; // 0
	MF	int B;    
	VL	double vlos;
	LD	double dopp;
	A	double aa;
	GM	int gm; //5
	AZI	int az;
	B0	double S0;
	B1	double S1;
	MC	double mac; //9
		double alfa;		
*/

//extern PRECISION gp4_gp2_rhoq[NLAMBDA],gp5_gp2_rhou[NLAMBDA],gp6_gp2_rhov[NLAMBDA];
extern float *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
extern float *gp4_gp2_rhoq,*gp5_gp2_rhou,*gp6_gp2_rhov;
extern float * gp1,*gp2,*dt,*dti,*gp3,*gp4,*gp5,*gp6,*etai_2;
extern float *dgp1,*dgp2,*dgp3,*dgp4,*dgp5,*dgp6,*d_dt;
extern float * d_ei,*d_eq,*d_eu,*d_ev,*d_rq,*d_ru,*d_rv;
extern PRECISION *dfi,*dshi;
extern PRECISION CC,CC_2,sin_gm,azi_2,sinis,cosis,cosis_2,cosi,sina,cosa,sinda,cosda,sindi,cosdi,sinis_cosa,sinis_sina;
extern PRECISION *fi_p,*fi_b,*fi_r,*shi_p,*shi_b,*shi_r;
extern float *etain,*etaqn,*etaun,*etavn,*rhoqn,*rhoun,*rhovn;
extern float *etai,*etaq,*etau,*etav,*rhoq,*rhou,*rhov;
extern float *parcial1,*parcial2,*parcial3;
extern float *nubB,*nupB,*nurB;
PRECISION **uuGlobalInicial;
PRECISION **HGlobalInicial;
PRECISION **FGlobalInicial;
extern int FGlobal,HGlobal,uuGlobal;
//extern PRECISION *G, *GMAC; // VECTOR WITH GAUSSIAN CREATED FOR CONVOLUTION 
extern PRECISION *GMAC; // VECTOR WITH GAUSSIAN CREATED FOR CONVOLUTION 
extern float *G;
//extern VSLConvTaskPtr taskConv;

extern fftw_complex * inSpectraFwMAC, *inSpectraBwMAC, *outSpectraFwMAC, *outSpectraBwMAC;
extern fftw_complex * inFilterMAC, * inFilterMAC_DERIV, * outFilterMAC, * outFilterMAC_DERIV;
extern fftw_plan planForwardMAC, planBackwardMAC;
extern fftw_plan planFilterMAC, planFilterMAC_DERIV;

extern fftw_complex * fftw_G_PSF, * fftw_G_MAC_PSF, * fftw_G_MAC_DERIV_PSF;
extern fftw_complex * inPSF_MAC, * inMulMacPSF, * inPSF_MAC_DERIV, *inMulMacPSFDeriv, *outConvFilters, * outConvFiltersDeriv;
extern fftw_plan planForwardPSF_MAC, planForwardPSF_MAC_DERIV,planBackwardPSF_MAC, planBackwardPSF_MAC_DERIV;
extern fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
extern fftw_plan planForwardPSF, planBackwardPSF;

int me_der(Cuantic *cuantic,Init_Model *initModel,PRECISION * wlines,PRECISION *lambda,int nlambda,float *d_spectra,float *spectra, float * spectra_slight, PRECISION ah,PRECISION * slight,int calcSpectra, int filter)
{

	int nterms,numl;
	int lineas;

//	PRECISION *etain,*etaqn,*etaun,*etavn,*rhoqn,*rhoun,*rhovn,etaaux;
//	PRECISION *etai,*etaq,*etau,*etav,*rhoq,*rhou,*rhov;

	PRECISION E00,LD,A,B1,MC,ALF;
	int il,i,j,k;
	PRECISION E0;	
//	PRECISION ulos,*nubB,*nupB,*nurB;


    int odd,ishift,par;
    
	
	E00=initModel->eta0; 
	//MF=initModel->B;
	LD=initModel->dopp;
	A=initModel->aa;
	B1=-((initModel->S1)*ah);
	MC=initModel->mac;
	ALF=initModel->alfa;

	nterms=NTERMS; 
	numl=nlambda;
	lineas=wlines[0];

	//DEFINO UN VECTOR DE DERIVADAS
	//POR ORDEN SERAN param=[eta0,magnet,vlos,landadopp,aa,gamma,azi]	
	
	Resetear_Valores_Intermedios(nlambda);

	
/*	for(j=0;j<numl;j++){
		etai[j]=1.0;
	}*/

	//a radianes	
/*	CC=PI/180.0;
	AZI=AZI*CC;
	GM=GM*CC;
	sinis=sin(GM)*sin(GM);
	cosis=cos(GM)*cos(GM);
	cosi=cos(GM);
	sina=sin(2*AZI);
	cosa=cos(2*AZI);

	sinda=cos(2*AZI)*2*CC;
	cosda=-sin(2*AZI)*2*CC;
	sindi=cos(GM)*sin(GM)*2*CC;
	cosdi=-sin(GM)*CC;
*/

	//reserva de memoria para vectores auxiliares
//    u=calloc(nlambda,sizeof(double));
    
    for(il=0;il<lineas;il++) {
		//Line strength
	    E0=E00*cuantic[il].FO; //y sino se definio Fo que debe de pasar 0 o 1 ...??

		fi_p=fi_p+nlambda*il*sizeof(float);
		fi_b=fi_b+nlambda*il*sizeof(float);
		fi_r=fi_r+nlambda*il*sizeof(float);
		shi_p=shi_p+nlambda*il*sizeof(float);
		shi_b=shi_b+nlambda*il*sizeof(float);
		shi_r=shi_r+nlambda*il*sizeof(float);

		
//			Leer_Puntero_Calculos_Compartidos(3,&nubB,&nupB,&nurB);

			

		//central component					    					
		funcionComponentFor(cuantic[il].N_PI,wlines[il+1],numl,cuantic[il].WEP,nupB,dfi,dshi,LD,A,0);

		//blue component
		funcionComponentFor(cuantic[il].N_SIG,wlines[il+1],numl,cuantic[il].WEB,nubB,dfi,dshi,LD,A,1);

		//red component
		funcionComponentFor(cuantic[il].N_SIG,wlines[il+1],numl,cuantic[il].WER,nurB,dfi,dshi,LD,A,2);
/*			free(nubB);
		free(nurB);
		free(nupB);*/
		

		//dispersion profiles				
		PRECISION E0_2;
//		Leer_Puntero_Calculos_Compartidos(3,&parcial1,&parcial2,&parcial3);

//		Leer_Puntero_Calculos_Compartidos(7,&etain,&etaqn,&etaun,&etavn,&rhoqn,&rhoun,&rhovn);

		//derivadas respecto de eta0 
		/*for(i=0;i<numl;i++){
			d_ei[i]=d_ei[i]+etain[i]/E0;
			d_eq[i]=d_eq[i]+etaqn[i]/E0;
			d_eu[i]=d_eu[i]+etaun[i]/E0;
			d_ev[i]=d_ev[i]+etavn[i]/E0;
			d_rq[i]=d_rq[i]+rhoqn[i]/E0;
			d_ru[i]=d_ru[i]+rhoun[i]/E0;
			d_rv[i]=d_rv[i]+rhovn[i]/E0;
		}*/

		for(i=0;i<numl;i++){
			d_ei[i]=d_ei[i]+etain[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_eq[i]=d_eq[i]+etaqn[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_eu[i]=d_eu[i]+etaun[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_ev[i]=d_ev[i]+etavn[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_rq[i]=d_rq[i]+rhoqn[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_ru[i]=d_ru[i]+rhoun[i]/E0;
		}
		for(i=0;i<numl;i++){
			d_rv[i]=d_rv[i]+rhovn[i]/E0;
		}

		float cosi_2_E0;
		sinis_cosa=E0*sinis_cosa/2;
		sinis_sina=E0*sinis_sina/2;
		E0_2=E0/2.0;
		float sindi_cosa,sindi_sina,cosdi_E0_2,cosis_2_E0_2,sinis_E0_2;
		sindi_cosa=sindi*cosa;
		sindi_sina=sindi*sina;
		cosdi_E0_2=(E0_2)*cosdi;
		cosi_2_E0=E0*cosi/2.0;
		cosis_2_E0_2=cosis_2*E0_2;
		sinis_E0_2=sinis*E0_2;

		
		/*printf("-------------------> d_fi :  %d\n",il);	
		printf("\n\nd_fi_c=[");
		int In,kk1;
		for(kk1=0;kk1<3;kk1++)		
			for(In=0;In<4;In++)		
				for(i=0;i<nlambda;i++){
					printf("%.16e\n",dfi[i + (In)*numl + (numl*4*kk1)]);										  
				}	
		printf("];\n\n\n");*/
	
		
		for(j=1;j<5;j++){
			//derivadas de los perfiles de dispersion respecto de B,VL,LDOPP,A 
			for(i=0;i<numl;i++){
				float dfisum,aux;
				dfisum=	dfi[i + (j-1)*numl+ (numl*4)]+dfi[i + (j-1)*numl + (numl*4*2)];
				
				d_ei[j*numl+i] = d_ei[j*numl+i] + (dfi[i+ (j-1)*numl] * sinis_E0_2 + dfisum * cosis_2_E0_2);
				
				aux=dfi[(j-1)*numl+i]-dfisum/2;
				d_eq[j*numl+i]=d_eq[j*numl+i]+(aux)*sinis_cosa;
				d_eu[j*numl+i]=d_eu[j*numl+i]+(aux)*sinis_sina;
			}
			for(i=0;i<numl;i++){
				d_ev[j*numl+i]= d_ev[j*numl+i] +(dfi[(j-1)*numl+i+(numl*4*2)]-dfi[(j-1)*numl+i+(numl*4)])*cosi_2_E0;
			}
		}
		for(j=1;j<5;j++){
			for(i=0;i<numl;i++){
				float aux=dshi[(j-1)*numl+i]-(dshi[(j-1)*numl+i+(numl*4)]+dshi[(j-1)*numl+i+(numl*4*2)])/2.0;
				d_rq[j*numl+i]=d_rq[j*numl+i]+(aux)*sinis_cosa;
				d_ru[j*numl+i]=d_ru[j*numl+i]+(aux)*sinis_sina;
			}
			for(i=0;i<numl;i++){
				d_rv[j*numl+i]=d_rv[j*numl+i]+((dshi[(j-1)*numl+i+(numl*4*2)]-dshi[(j-1)*numl+i+(numl*4)]))*cosi_2_E0;
			}	
		}
		
		//derivadas de los perfiles de dispersion respecto de GAMMA
		float cosi_cosdi,sindi_E0_2;
		cosi_cosdi=cosi*cosdi*E0_2;
		sindi_E0_2=sindi*E0_2;
		for(i=0;i<numl;i++)
			d_ei[5*numl+i]=d_ei[5*numl+i]+fi_p[i]*sindi_E0_2+(parcial1[i])*cosi_cosdi;
		for(i=0;i<numl;i++)
			d_eq[5*numl+i]=d_eq[5*numl+i]+parcial2[i]*sindi_cosa;
		for(i=0;i<numl;i++)
			d_eu[5*numl+i]=d_eu[5*numl+i]+parcial2[i]*sindi_sina;
		for(i=0;i<numl;i++)
			d_ev[5*numl+i]=d_ev[5*numl+i]+(fi_r[i]-fi_b[i])*cosdi_E0_2;
		

		for(i=0;i<numl;i++)
			d_rq[5*numl+i]=d_rq[5*numl+i]+parcial3[i]*sindi_cosa;
		for(i=0;i<numl;i++)
			d_ru[5*numl+i]=d_ru[5*numl+i]+parcial3[i]*sindi_sina;
		for(i=0;i<numl;i++)
			d_rv[5*numl+i]=d_rv[5*numl+i]+(shi_r[i]-shi_b[i])*cosdi_E0_2;
		

		//derivadas de los perfiles de dispersion respecto de AZI
		float sinis_cosda,sinis_sinda;
		sinis_cosda=sinis*cosda;
		sinis_sinda=sinis*sinda;		
		for(i=0;i<numl;i++){				
			d_eq[6*numl+i]=d_eq[6*numl+i]+parcial2[i]*sinis_cosda;
		}

		for(i=0;i<numl;i++){
			d_eu[6*numl+i]=d_eu[6*numl+i]+parcial2[i]*sinis_sinda;
		}
		for(i=0;i<numl;i++){
			d_rq[6*numl+i]=d_rq[6*numl+i]+parcial3[i]*sinis_cosda;
		}
		for(i=0;i<numl;i++){
			d_ru[6*numl+i]=d_ru[6*numl+i]+parcial3[i]*sinis_sinda;
		}

/*		free(parcial1);
		free(parcial2);
		free(parcial3);
*/
	
/*		free(fi_p);
		free(fi_b);
		free(fi_r);
		free(shi_p);
		free(shi_b);
		free(shi_r);
*/
/*		free(etain);
		free(etaqn);
		free(etaun);
		free(etavn);
		free(rhoqn);
		free(rhoun);
		free(rhovn);
*/
	} //end for

            	
//	Leer_Puntero_Calculos_Compartidos(7,&etai,&etaq,&etau,&etav,&rhoq,&rhou,&rhov);

    //Los parametros de Stokes estan normalizados a la intensidad en el continuo (no)

    //bucle para derivadas de I,Q,U,V 
    //derivada de spectra respecto  E0,MF,VL,LD,A,gamma,azi

	//static PRECISION dtaux[nlambda];
	//PRECISION *dtaux,*etai_gp3,*ext1,*ext2,*ext3,*ext4;
	//PRECISION dtaux[nlambda], etai_gp3[nlambda],ext1[nlambda],ext2[nlambda],ext3[nlambda],ext4[nlambda];
	
	/*dtaux = calloc(nlambda,sizeof(PRECISION));
	etai_gp3 = calloc(nlambda,sizeof(PRECISION));
	ext1 = calloc(nlambda,sizeof(PRECISION));
	ext2 = calloc(nlambda,sizeof(PRECISION));
	ext3 = calloc(nlambda,sizeof(PRECISION));
	ext4 = calloc(nlambda,sizeof(PRECISION));*/

   for(i=0;i<numl;i++)
		dtaux[i]=(B1)/(dt[i]*dt[i]);

   for(i=0;i<numl;i++){
		etai_gp3[i]=etai[i]*gp3[i];
	}

   for(i=0;i<numl;i++){
		float aux=2*etai[i];
		ext1[i]=aux*etaq[i]+etav[i]*rhou[i]-etau[i]*rhov[i];
		ext2[i]=aux*etau[i]+etaq[i]*rhov[i]-etav[i]*rhoq[i];
		ext3[i]=aux*etav[i]+etau[i]*rhoq[i]-etaq[i]*rhou[i];
		ext4[i]=aux*gp1[i];
	}

	/*
	printf("-------------------> d_ei :  %d\n",il);	
	printf("\n\nd_ei_c=[");
	int In;
	for(In=0;In<NTERMS;In++)		
		for(i=0;i<nlambda;i++){
			printf("%.16e\n",d_ei[i+nlambda*In]);										  
		}	
	printf("];\n\n\n");
	*/
	
	
	//printf("-------------------> dpgx :  %d\n",il);
    for(il=0;il<7;il++){


    	for(i=0;i<numl;i++){
    		dgp1[i]=2.0*(etai[i]*d_ei[i+numl*il]-etaq[i]*d_eq[i+numl*il]-etau[i]*d_eu[i+numl*il]-etav[i]*d_ev[i+numl*il]  
				 +rhoq[i]*d_rq[i+numl*il]+rhou[i]*d_ru[i+numl*il]+rhov[i]*d_rv[i+numl*il]);
    	}
    
    	for(i=0;i<numl;i++){
    		dgp2[i]=rhoq[i]*d_eq[i+numl*il]+etaq[i]*d_rq[i+numl*il]+rhou[i]*d_eu[i+numl*il]+etau[i]*d_ru[i+numl*il]+
    		                    rhov[i]*d_ev[i+numl*il]+etav[i]*d_rv[i+numl*il];
    	}
	
    	for(i=0;i<numl;i++){
    		d_dt[i]=ext4[i]*d_ei[i+numl*il]+etai_2[i]*dgp1[i]-2*gp2[i]*dgp2[i];
    	}
		
    	for(i=0;i<numl;i++){
    		dgp3[i]=2.0*(etai[i]*d_ei[i+numl*il]+rhoq[i]*d_rq[i+numl*il]+rhou[i]*d_ru[i+numl*il]+rhov[i]*d_rv[i+numl*il]);
    	}

    	for(i=0;i<numl;i++){    		
    		d_spectra[i+numl*il]=-(((d_ei[i+numl*il]*gp3[i]+etai[i]*dgp3[i])*dt[i]-d_dt[i]*etai_gp3[i])*(dtaux[i]));
    	}
	
    	for(i=0;i<numl;i++){
    		dgp4[i]=d_ei[i+numl*il]*(ext1[i])+(etai_2[i])*d_eq[i+numl*il]+
    		etai[i]*(rhou[i]*d_ev[i+numl*il]+etav[i]*d_ru[i+numl*il]-rhov[i]*d_eu[i+numl*il]-etau[i]*d_rv[i+numl*il]);
    	}
		
    	for(i=0;i<numl;i++){
    		d_spectra[i+numl*il+numl*nterms]=((dgp4[i]+d_rq[i+numl*il]*gp2[i]+rhoq[i]*dgp2[i])*dt[i]-
    				d_dt[i]*(gp4_gp2_rhoq[i]))*(dtaux[i]);
    	}    
	
    	for(i=0;i<numl;i++){
    		dgp5[i]=d_ei[i+numl*il]*(ext2[i])+(etai_2[i])*d_eu[i+numl*il]+
    		etai[i]*(rhov[i]*d_eq[i+numl*il]+etaq[i]*d_rv[i+numl*il]-rhoq[i]*d_ev[i+numl*il]-etav[i]*d_rq[i+numl*il]);
    	}

    	for(i=0;i<numl;i++){
    		d_spectra[i+numl*il+(numl*nterms*2)]=((dgp5[i]+d_ru[i+numl*il]*gp2[i]+rhou[i]*dgp2[i])*dt[i]-
    				d_dt[i]*(gp5_gp2_rhou[i]))*(dtaux[i]);
    	}    

    	for(i=0;i<numl;i++){
    		dgp6[i]=d_ei[i+numl*il]*(ext3[i])+(etai_2[i])*d_ev[i+numl*il]+
    		etai[i]*(rhoq[i]*d_eu[i+numl*il]+etau[i]*d_rq[i+numl*il]-rhou[i]*d_eq[i+numl*il]-etaq[i]*d_ru[i+numl*il]);
    	}

		for(i=0;i<numl;i++){
    		d_spectra[i+numl*il+(numl*nterms*3)]=((dgp6[i]+d_rv[i+numl*il]*gp2[i]+rhov[i]*dgp2[i])*dt[i]-
    				d_dt[i]*(gp6_gp2_rhov[i]))*(dtaux[i]);
    	} 

		
		 // printf("-------------------> IL :  %d\n",il);
		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp1[i]);
			 
		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp2[i]);

		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp3[i]);
			 
		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp4[i]);

		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp5[i]);
			 
		 // for(i=0;i<numl;i++)
			 // printf("%.8e\n",dgp6[i]);
			 

		/*
		//printf("-------------------> dspectra :  %d\n",il);
		for(i=0;i<numl;i++)
			printf("%.8e\n",d_spectra[i+numl*il+(numl*nterms*0)]);			 
						 
		for(i=0;i<numl;i++)
			printf("%.8e\n",d_spectra[i+numl*il+(numl*nterms*1)]);			 
			 
		for(i=0;i<numl;i++)
			printf("%.8e\n",d_spectra[i+numl*il+(numl*nterms*2)]);			 

		for(i=0;i<numl;i++)
			printf("%.8e\n",d_spectra[i+numl*il+(numl*nterms*3)]);			 
		*/	 
    }
	//printf("------------------->  end dpgx :  %d\n",il);
	// printf("-------------------> end dspectra :  %d\n",il);

	
    //LA 7-8 RESPECTO B0 Y B1

	/*for(i=0;i<numl;i++){
		dti[i]=-(dti[i]*ah);		
		d_spectra[i+numl*8]=-dti[i]*etai_gp3[i];
		d_spectra[i+numl*8+(numl*nterms)]= dti[i]*(gp4_gp2_rhoq[i]);   		
		d_spectra[i+numl*8+(numl*nterms*2)]= dti[i]*(gp5_gp2_rhou[i]);
		d_spectra[i+numl*8+(numl*nterms*3)]= dti[i]*(gp6_gp2_rhov[i]);
		// S0 
		d_spectra[i+numl*7]=1;
    	d_spectra[i+numl*7+(numl*nterms)]=0;
    	d_spectra[i+numl*7+(numl*nterms*2)]=0;
    	d_spectra[i+numl*7+(numl*nterms*3)]=0;
		// azimuth stokes I &V
    	d_spectra[i+numl*6]=0;
    	d_spectra[i+numl*6+(numl*nterms*3)]=0;
	}*/
    for(i=0;i<numl;i++)
		dti[i]=-(dti[i]*ah);

    for(i=0;i<numl;i++){
    	d_spectra[i+numl*8]=-dti[i]*etai_gp3[i];
    }

    for(i=0;i<numl;i++){
    	d_spectra[i+numl*8+(numl*nterms)]= dti[i]*(gp4_gp2_rhoq[i]);   		
    }

    for(i=0;i<numl;i++){
    	d_spectra[i+numl*8+(numl*nterms*2)]= dti[i]*(gp5_gp2_rhou[i]);
    }

    for(i=0;i<numl;i++){
    	d_spectra[i+numl*8+(numl*nterms*3)]= dti[i]*(gp6_gp2_rhov[i]);
    }

	//S0
   for(i=0;i<numl;i++){
    	d_spectra[i+numl*7]=1;
    	d_spectra[i+numl*7+(numl*nterms)]=0;
    	d_spectra[i+numl*7+(numl*nterms*2)]=0;
    	d_spectra[i+numl*7+(numl*nterms*3)]=0;
	}

	//azimuth stokes I &V
    for(i=0;i<numl;i++){
    	d_spectra[i+numl*6]=0;
    	d_spectra[i+numl*6+(numl*nterms*3)]=0;				  
	}
	
	
/*        
    free(etai);
    free(etaq);
    free(etau);
    free(etav);

    free(rhoq);
    free(rhou);
    free(rhov);
*/
	/*free(etai_gp3);
	free(ext1);
	free(ext2);
	free(ext3);
	free(ext4);
	free(dtaux);*/
    //MACROTURBULENCIA
                
	 int macApplied = 0;
    if(MC > 0.0001){
		 
		macApplied = 1;
		odd=(numl%2);		
		int startShift = numl/2;
		if(odd) startShift+=1;		
		

//		printf("dentro MC \n");
    	//la 9 respecto MC
    	//convolucion del espectro original
    	//PRECISION *g1 =fgauss(MC,lambda,numl,wlines[1],0);  // Gauss Function
		fgauss(MC,lambda,numl,wlines[1],0);  // Gauss Function
		// CALCULATE DERIV OF G1 
		PRECISION ild = (wlines[1] * MC) / 2.99792458e5; //Sigma
		PRECISION centro = lambda[(int)numl / 2];		  //center of the axis
    	//PRECISION * g2=fgauss(MC,lambda,numl,wlines[1],1);  // Derivative of the Gauss F. with respect macro
		// CALCULATE FFT OF G1 

		if(filter){// if there is PSF filter convolve both gaussian and use the result as the signal to convolve
			for(i=0;i<numl;i++){ // copy gmac to
				inPSF_MAC[i] = (GMAC[i]) + 0 * _Complex_I;
				inPSF_MAC_DERIV[i] = (GMAC[i] / MC * ((((lambda[i] - centro) / ild) * ((lambda[i] - centro) / ild)) - 1.0)) + 0 * _Complex_I;
			}
			fftw_execute(planForwardPSF_MAC);
			fftw_execute(planForwardPSF_MAC_DERIV);
			/*for(i=0;i<numl;i++){ // multiply both fft gaussians
				outConvFilters[i] = fftw_G_PSF[i] * (fftw_G_MAC_PSF[i]/numl);
				outConvFiltersDeriv[i] = fftw_G_PSF[i] * (fftw_G_MAC_DERIV_PSF[i]/numl);
			}			
			for(i=0,ishift=numl/2;i<numl/2;i++,ishift++){
				outFilterMAC[ishift]= outConvFilters[i]*numl;
				outFilterMAC_DERIV[ishift]= outConvFiltersDeriv[i]*numl;
			}
			for(i=numl/2,ishift=0;i<numl;i++,ishift++){
				outFilterMAC[ishift]= outConvFilters[i]*numl;
				outFilterMAC_DERIV[ishift]= outConvFiltersDeriv[i]*numl;
			}*/
			for(i=0;i<numl;i++){ // multiply both fft gaussians
				inMulMacPSF[i] = fftw_G_PSF[i] * (fftw_G_MAC_PSF[i]/numl);
				inMulMacPSFDeriv[i] = fftw_G_PSF[i] * (fftw_G_MAC_DERIV_PSF[i]/numl);
			}
			fftw_execute(planBackwardPSF_MAC);
			fftw_execute(planBackwardPSF_MAC_DERIV);
			for(i=0,ishift=startShift;i<numl/2;i++,ishift++){
				inFilterMAC[ishift]= outConvFilters[i]*numl;
				inFilterMAC_DERIV[ishift]= outConvFiltersDeriv[i]*numl;
			}
			for(i=(numl/2),ishift=0;i<numl;i++,ishift++){
				inFilterMAC[ishift]= outConvFilters[i]*numl;
				inFilterMAC_DERIV[ishift]= outConvFiltersDeriv[i]*numl;
			}
		}
		else{
			for(i=0;i<numl;i++){
				inFilterMAC[i] = GMAC[i] + 0 * _Complex_I;
				//inFilterMAC_DERIV[i] = g2[i]  + 0 * _Complex_I;
				inFilterMAC_DERIV[i] = (GMAC[i] / MC * ((((lambda[i] - centro) / ild) * ((lambda[i] - centro) / ild)) - 1.0)) + 0 * _Complex_I;
			}
		}
		fftw_execute(planFilterMAC);
		fftw_execute(planFilterMAC_DERIV);

		
    	for(il=0;il<4;il++){
    		//if(fabs(mean(spectra+numl*il,numl)) >= 1.e-25){
				for(i=0;i<numl;i++){
					inSpectraFwMAC[i] = spectra[numl*il+i] + 0 * _Complex_I;
				} 
				fftw_execute(planForwardMAC);
				for(i=0;i<numl;i++){
					//outSpectraFwMAC[i] = outSpectraFwMAC[i]/numl;
					inSpectraBwMAC[i]=(outSpectraFwMAC[i]/numl)*(outFilterMAC_DERIV[i]/numl);
    			}
				fftw_execute(planBackwardMAC);
				//shift: -numl/2
				for(i=0,ishift=startShift;i<numl/2;i++,ishift++){
					d_spectra[ishift+9*numl+numl*nterms*il]=creal(outSpectraBwMAC[i])*numl;
				}
				for(i=(numl/2),ishift=0;i<numl;i++,ishift++){
					d_spectra[ishift+9*numl+numl*nterms*il]=creal(outSpectraBwMAC[i])*numl;
				}  					
    			if(calcSpectra){
					for(i=0;i<numl;i++){
						inSpectraBwMAC[i]=(outSpectraFwMAC[i]/numl)*(outFilterMAC[i]/numl);
					}
					fftw_execute(planBackwardMAC);
					//shift: -numl/2
					for(i=0,ishift=numl/2;i<numl/2;i++,ishift++){
						spectra[ishift+il*numl]=creal(outSpectraBwMAC[i])*numl;
					}
					for(i=numl/2,ishift=0;i<numl;i++,ishift++){
						spectra[ishift+il*numl]=creal(outSpectraBwMAC[i])*numl;
					} 
				}  		
    		//}
    	}
		
		
	   for(par=0;par<4;par++){
    		//Va hasta 8 porque la macro no la convoluciono
    		//seria directamente dmacro=I*dg    		
	    	for(il=0;il<9;il++){
				if(il!=7){
	    			//if(fabs(mean(d_spectra+numl*il+numl*nterms*par,numl)) >= 1.e-25){
						for(i=0;i<numl;i++){
							inSpectraFwMAC[i] = d_spectra[(numl*il+numl*nterms*par)+i] + 0 * _Complex_I;
						} 
						fftw_execute(planForwardMAC);
						for(i=0;i<numl;i++){
							inSpectraBwMAC[i]=(outSpectraFwMAC[i]/numl)*(outFilterMAC[i]/numl);
						}
						fftw_execute(planBackwardMAC);  			

						//shift 
						for(i=0,ishift=startShift;i<numl/2;i++,ishift++){
							d_spectra[ishift+il*numl+numl*nterms*par]=creal(outSpectraBwMAC[i])*numl;
						}
						for(i=(numl/2),ishift=0;i<numl;i++,ishift++){
							d_spectra[ishift+il*numl+numl*nterms*par]=creal(outSpectraBwMAC[i])*numl;
						}  
					//}
	    		}
	    	}
    	}
		//free(g1);
		//free(g2);

   
    }//end if(MC > 0.0001)


	// stray light factor 

	if(slight!=NULL){

		/*printf("\nSLIGHT A APLICAR \n");
		for(i=0;i<numl;i++){
			printf("%lf\t%lf\t%lf\t%lf\n",slight[i],slight[i+numl],slight[i+(numl*2)],slight[i+(numl*3)]);
		}
		printf("\n");*/

		// Response Functions 
	   for(par=0;par<NPARMS;par++){
	    	for(il=0;il<NTERMS;il++){
				for(i=0;i<numl;i++){
					d_spectra[numl*il+numl*nterms*par+i]=d_spectra[numl*il+numl*nterms*par+i]*ALF;

					if(il==10){ //Magnetic filling factor Response function
						d_spectra[numl*il+numl*nterms*par+i]=spectra_slight[numl*par+i]-slight[numl*par+i];
					}
				}
	    	}

		
    	}
		/*for(i=0;i<numl*NPARMS;i++){
			spectra[i] = spectra[i]*ALF+slight[i]*(1.0-ALF);
		}*/
	}
	if(!macApplied && filter){
		int h;
		int odd=(numl%2);
		int startShift = (numl)/2;
		if(odd) startShift+=1;

		for (j = 0; j < NPARMS; j++)
		{
			for (i = 0; i < NTERMS; i++)
			{
				if (i != 7)	{																														 //no convolucionamos S0
					// copy to inSpectra
					for(k=0;k<(numl);k++){
						inSpectraFwPSF[k] = d_spectra[(numl * i + numl * NTERMS * j) + k] + 0 * _Complex_I;
					}
					fftw_execute(planForwardPSF);
					for(h=0;h<numl;h++){
						inSpectraBwPSF[h] = (outSpectraFwPSF[h]/numl) * fftw_G_PSF[h];
					}
					fftw_execute(planBackwardPSF);   			
					//shift 
					for(h=0,ishift=startShift;h<numl/2;h++,ishift++){
						d_spectra[ishift+ numl * i + numl * NTERMS * j]=creal(outSpectraBwPSF[h])*numl;
					}
					for(h=(numl/2),ishift=0;h<numl;h++,ishift++){
						d_spectra[ishift+numl * i + numl * NTERMS * j]=creal(outSpectraBwPSF[h])*numl;
					}
				}
			}
		}

		//response_functions_convolution(&nlambda);
	}
	ResetPointerShareCalculation();
	
	return 1;
	
}


/*
 * 
 */
int funcionComponentFor(int n_pi,PRECISION iwlines,int numl,PRECISION *wex,float *nuxB,PRECISION *dfi,PRECISION *dshi,PRECISION LD,PRECISION A,int desp)
{
	PRECISION *uu;
	int i,j;
	PRECISION dH_u[numl],dF_u[numl],auxCte[numl];
	
	PRECISION *H,*F;
	
	/*PRECISION *dH_u,*dF_u,*auxCte;
	dH_u = malloc(numl * sizeof(PRECISION));
	dF_u = malloc(numl * sizeof(PRECISION));
	auxCte = malloc(numl * sizeof(PRECISION));*/
	

	for(j=0;j<numl;j++){
		auxCte[j]=(-iwlines)/(VLIGHT*LD);
	}

	

//	H=NULL;
//	F=NULL;
//	printf("MEDER\n");

	//component
	for(i=0;i<n_pi;i++){

/*		uu=uuGlobal+NLAMBDA*sizeof(PRECISION)*i;
		H=HGlobal+NLAMBDA*sizeof(PRECISION)*i;
		F=FGlobal+NLAMBDA*sizeof(PRECISION)*i;*/
		uu=uuGlobalInicial[uuGlobal+i];
		F=FGlobalInicial[HGlobal+i];
		H=HGlobalInicial[FGlobal+i];

/*	printf("uu %d\n",uu);
	printf("H %d\n",H);
	printf("F %d\n",F);
	*/			
//		Leer_Puntero_Calculos_Compartidos(2,&H,&F);
//		printf("aaa mi3 here nto for uu %f\n",H[0]);
		for(j=0;j<numl;j++){
			dH_u[j]=((4*A*F[j])-(2*uu[j]*H[j]))*wex[i];
		}

		for(j=0;j<numl;j++){
			dF_u[j]=(RR-A*H[j]-2*uu[j]*F[j])*wex[i];//
		}

		/*
		printf("DFU %d\n",desp);
		for(j=0;j<numl;j++)
			printf("%.16e\n",dF_u[j]);
		*/
		
		for(j=0;j<numl;j++){
			uu[j]=-uu[j]/LD;
		}

		//numl*4*3
		// a   b c
		//col a,fil b,desplazamiento c
		//b*a+numl+(numl*4*c)
		//dfi
		for(j=0;j<numl;j++){
			//(*,0,desp)=>0*numl+j+(numl*4*0)
			dfi[j+(numl*4*desp)]=dfi[j+(numl*4*desp)]+dH_u[j]*(-nuxB[i]);
		}
		
		for(j=0;j<numl;j++){
			//(*,1,desp)=>1*numl+j+(numl*4*0)
			dfi[numl+j+(numl*4*desp)]=dfi[numl+j+(numl*4*desp)]+dH_u[j]*auxCte[j];
		}

		for(j=0;j<numl;j++){
			//(*,2,desp)=>1*numl+j+(numl*4*0)
			dfi[2*numl+j+(numl*4*desp)]=dfi[2*numl+j+(numl*4*desp)]+(dH_u[j]*uu[j]); 											
		}								

		for(j=0;j<numl;j++){
			//(*,3,desp)=>1*numl+j+(numl*4*0)
			dfi[3*numl+j+(numl*4*desp)]=dfi[3*numl+j+(numl*4*desp)]+(-2*dF_u[j]);//dH_a[j];						
		}
		
		//dshi
		for(j=0;j<numl;j++){
			//(*,0,desp)=>0*numl+j+(numl*4*0)
			dshi[j+(numl*4*desp)]=dshi[j+(numl*4*desp)]+(dF_u[j])*(-nuxB[i]);
		}
						
		for(j=0;j<numl;j++){
			//(*,1,desp)=>1*numl+j+(numl*4*0)
			dshi[numl+j+(numl*4*desp)]=dshi[numl+j+(numl*4*desp)]+dF_u[j]*auxCte[j];
		}

		for(j=0;j<numl;j++){
			//(*,2,desp)=>1*numl+j+(numl*4*0)
			dshi[2*numl+j+(numl*4*desp)]=dshi[2*numl+j+(numl*4*desp)]+(dF_u[j]*uu[j]); 											
		}								
		
		for(j=0;j<numl;j++){
			//(*,3,desp)=>1*numl+j+(numl*4*0)
			dshi[3*numl+j+(numl*4*desp)]=dshi[3*numl+j+(numl*4*desp)]+(dH_u[j]/2);						
		}									

//		free(uu);
//		free(H);
//		free(F);

	}


	uuGlobal=uuGlobal+n_pi;
	HGlobal=HGlobal+n_pi;
	FGlobal=FGlobal+n_pi;

	return 1;	
}


void Resetear_Valores_Intermedios(int nlambda){
	

	
	memset(d_ei , 0, (nlambda*7)*sizeof(float));
	memset(d_eq , 0, (nlambda*7)*sizeof(float));
	memset(d_ev , 0, (nlambda*7)*sizeof(float));
	memset(d_eu , 0, (nlambda*7)*sizeof(float));
	memset(d_rq , 0, (nlambda*7)*sizeof(float));
	memset(d_ru , 0, (nlambda*7)*sizeof(float));
	memset(d_rv , 0, (nlambda*7)*sizeof(float));
	memset(dfi , 0, (nlambda*4*3)*sizeof(PRECISION));
	memset(dshi , 0, (nlambda*4*3)*sizeof(PRECISION));

	/*
	int i=0;
	for(i=0;i<nlambda*7;i++){
		d_ei[i]=0;
	}

	for(i=0;i<nlambda*7;i++){
		d_eq[i]=0;
	}

	for(i=0;i<nlambda*7;i++){
		d_ev[i]=0;
	}

	for(i=0;i<nlambda*7;i++){
		d_eu[i]=0;
	}

	for(i=0;i<nlambda*7;i++){
		d_rq[i]=0;
	}
	for(i=0;i<nlambda*7;i++){
		d_ru[i]=0;
	}
	for(i=0;i<nlambda*7;i++){
		d_rv[i]=0;
	}
	for(i=0;i<nlambda*4*3;i++)
		dfi[i]=0;

	for(i=0;i<nlambda*4*3;i++)
		dshi[i]=0;

	*/


}
