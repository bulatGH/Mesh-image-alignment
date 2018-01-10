// Cpp_Project2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "AllTypes.h"
#include "BasicClasses.h"
#include "Filtering.h"
#include "Object_3D.h"
#include "Procrustes.h"
#include "Segmenter.h"
#include <math.h> 


float newScale = 0.5f;

void DentalCasts_Test()
{
	std::srand(std::time(0));
	bool downCheck = false;
	std::string tImageFile  = "G:\\MY PROJECTS\\_x_Teeth and Jaws\\Teeth\\Training Images\\Original\\Rouseau Elva Max\\Rouseau, Elva - Max Final Plan.pr";
	std::string tMeshFile   = "G:\\MY PROJECTS\\_x_Teeth and Jaws\\Teeth\\Training Images\\Original\\Rouseau Elva Max\\Rouseau_Cast Model_Up.stl";
	std::string tMeshResult = "G:\\MY PROJECTS\\_x_Teeth and Jaws\\Teeth\\Training Images\\Original\\Rouseau Elva Max\\Rouseau_Cast Model_Up_Result.stl";
	
	DWORD start = GetTickCount();
	
	// read mesh with rescaling it
	shared_ptr<Mesh> mainObject = Mesh::ReadSTLMeshFlip(tMeshFile.c_str(), 1.0f / newScale);
	
	std::vector<bool> tselectedVertices;
	tselectedVertices.resize(mainObject->Vertices.size());
	for (size_t i = 0; i < mainObject->Vertices.size(); i++)
	{
		tselectedVertices.push_back(true);
	}

	shared_ptr<sImages3D_S> mainImage1 = ReadStraightforwardly::ReadOneXLFormat(tImageFile);

	// I rescale image to isotropic with pixel size of "newScale". Do not forget that we also need to slightly modify the origin after rescaling
	float tnewOriginX = mainImage1->Origin[0] + (newScale / 2 - mainImage1->Scale[0] / 2);
	float tnewOriginY = mainImage1->Origin[1] + (newScale / 2 - mainImage1->Scale[1] / 2);
	float tnewOriginZ = mainImage1->Origin[2] + (newScale / 2 - mainImage1->Scale[2] / 2);
	shared_ptr<Images3D_S> mainImageNorm = ImageRescaling::RescaleISO(mainImage1, newScale)->CreateNormalizedImages(-1000, 2000);//->CreateNormalizedImages(2000);

	//I create a gradient image to optimize mesh adaptaion
	shared_ptr<GradientImage_3D> mainGradImage = GradientFiltering::GradientImage(mainImageNorm);
	
	//It is improtant to remove origin information from the mesh representation
	mainObject->MoveObject(-tnewOriginY / newScale, -tnewOriginX / newScale, -tnewOriginZ / newScale);
	

	shared_ptr<RigidRegistration> segmenter(new RigidRegistration(mainObject));
	double tchange = 100000; //this variable indicates how much the mesh has moved at each step.
	int tsteps = 0; 
	int blockStop = 100; // in case the mesh did not come to a stable position after "blockStop" steps, we optimization is terminated
	while (tchange > 0.25)
    {
        tchange = segmenter->DetectOptimalPointsSelectedVertices(mainGradImage, tselectedVertices);
		
		tsteps++;
		if(tsteps > blockStop)
		{
			break;
		}
	}
	

	DWORD elapsed = GetTickCount() - start;
	printf ("It took (%f seconds) to register mesh.\n",((float)elapsed)/CLOCKS_PER_SEC);

	//return origin value back to the mesh definition
	mainObject->MoveObject(tnewOriginY / newScale, tnewOriginX / newScale, tnewOriginZ / newScale);
    mainObject->ScaleBrutally(newScale, newScale, newScale);
	mainObject->SaveAsSTL(tMeshFile.c_str(), tMeshResult.c_str());
	printf ("finshed.\n");
	std::getchar();
	return;
}

int _tmain(int argc, _TCHAR* argv[])
{
	clock_t time1 = clock();
	DentalCasts_Test();
	clock_t time2 = clock();
	double tval = double(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << tval << '\n';

	std::cin >> tval;
	return 0;
}

