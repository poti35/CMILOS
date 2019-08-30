#include "mpiTools.h"
#include "defines.h"

void buildMpiType( MPI_Datatype * mpiInitModel){ 
   	const int nitemsStructInitModel = 11;
	int blocklenghtInitModel [11] = {1,1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype typesInitModel [11] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
		
	MPI_Aint offsetsInitModel [11];
	offsetsInitModel[0] = offsetof(Init_Model, eta0);
	offsetsInitModel[1] = offsetof(Init_Model, B);
	offsetsInitModel[2] = offsetof(Init_Model, vlos);
	offsetsInitModel[3] = offsetof(Init_Model, dopp);
	offsetsInitModel[4] = offsetof(Init_Model, aa);
	offsetsInitModel[5] = offsetof(Init_Model, gm);
	offsetsInitModel[6] = offsetof(Init_Model, az);
	offsetsInitModel[7] = offsetof(Init_Model, S0);
	offsetsInitModel[8] = offsetof(Init_Model, S1);
	offsetsInitModel[9] = offsetof(Init_Model, mac);
	offsetsInitModel[10] = offsetof(Init_Model, alfa);
	MPI_Type_create_struct(nitemsStructInitModel, blocklenghtInitModel, offsetsInitModel, typesInitModel, mpiInitModel);
	MPI_Type_commit(mpiInitModel);

	/* CREATE A TYPE FOR STRUCT VPIXEL  */

	/* const int nitemsStructVPixel = 3;
	int blocklenghtVPixel [3] = {nlambda,nlambda*NPARMS,1};
	MPI_Datatype typesVPixel [3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
	
	MPI_Aint offsetsVPixels [3];
	offsetsVPixels[0] = offsetof(vpixels, vLambda);
	offsetsVPixels[1] = offsetof(vpixels, spectro);
	offsetsVPixels[2] = offsetof(vpixels, nLambda);
	MPI_Type_create_struct(nitemsStructVPixel, blocklenghtVPixel, offsetsVPixels, typesVPixel, mpiVPixels);
	MPI_Type_commit(mpiVPixels);
	*/
	/* CREATE  A TYPE FOR STRUCT CUANTIC*/ 

	/* const int nitemsStructCuantic = 12;
	int blocklenghtCuantic [12] = {1,1, N_SIG, N_PI, N_SIG, N_SIG, N_PI, N_SIG,1,1,1,1};
	MPI_Datatype typesCuantic [12] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	
	MPI_Aint offsetsCuantic [12];
	offsetsCuantic[0] = offsetof(Cuantic, N_PI);
	offsetsCuantic[1] = offsetof(Cuantic, N_SIG);
	offsetsCuantic[2] = offsetof(Cuantic, NUB);
	offsetsCuantic[2] = offsetof(Cuantic, NUP);
	offsetsCuantic[2] = offsetof(Cuantic, NUR);
	offsetsCuantic[2] = offsetof(Cuantic, WEB);
	offsetsCuantic[2] = offsetof(Cuantic, WEP);
	offsetsCuantic[2] = offsetof(Cuantic, WER);
	offsetsCuantic[2] = offsetof(Cuantic, GL);
	offsetsCuantic[2] = offsetof(Cuantic, GU);
	offsetsCuantic[2] = offsetof(Cuantic, GEFF);
	offsetsCuantic[2] = offsetof(Cuantic, FO);
	MPI_Type_create_struct(nitemsStructCuantic, blocklenghtCuantic, offsetsCuantic, typesCuantic, mpiCuantic);
	MPI_Type_commit(mpiCuantic);
	*/

}


