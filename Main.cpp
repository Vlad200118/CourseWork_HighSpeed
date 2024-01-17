#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "GasModel.h"
#include "Variables.h"
#include "InputOutputData.h"
#include "Mesh.h"
#include "Functions.h"
#include "BoundaryConditions.h"
#include "Initialization.h"
#include "MonitorPoints.h"
#include "FluxScheme.h"
#include "TVDSchemes.h"

int main()
{
	InputData inputData("input.txt");
	OutputData outputData;

	inputData.Write(std::cout);

	
	GasModel gas(inputData.Cp(), inputData.gamma(), inputData.molMass());



	std::vector<ConservativeVar> w;
	std::vector<ConservativeVar> wNext;
	std::vector<vectorFlux> flux;
	vectorFlux  freeFlux;


	ResizeVectors(inputData.NX() + 2, w, wNext);
	ResizeVectors(inputData.NX() + 1, flux);

	//Calculate MESH 
	Mesh mesh;
	mesh.BuildMesh(inputData.NX(), inputData.Length());

	//MONITORING POINTS
	MonitorPoints monPoints(mesh);


	BoundaryConditions boundCond (inputData.NX(), gas, inputData.BCFileName());


	FluxScheme* fluxScheme;

	if (inputData.FluxScheme() == 1)
	{	
		fluxScheme = new AUSM_Standard(gas);
	}
	else if (inputData.FluxScheme() == 2)
	{
		fluxScheme = new AUSM_Plus(gas);
	}
	else 
	{
		fluxScheme = new Godunov(gas);
	}

	TVDScheme* tvdScheme;
	if (inputData.TVDScheme() == 0)
	{
		tvdScheme = new TVDScheme(false);
	}
	else
	{
		tvdScheme = new TVDScheme(true, (TVDScheme::TVDFunctions)inputData.TVDScheme());
	}

	

	//INITIALIZATION
	Initialization init;
	init.initialize(mesh, w, gas, inputData.initialField(InputData::left), inputData.initialField(InputData::right), inputData.diaphragmPos());

	outputData.WriteData("initialField.txt", mesh, w, gas, naturalVariables);
	
	
	double dX = mesh.GetMinStep();
	double deltaTime = inputData.CFL() * dX / GetVelocityForCFL(inputData, gas);
	size_t NT = size_t(inputData.TimeEnd() / deltaTime);
	std::cout << "delta Time = " << deltaTime << std::endl;
	std::cout << "NT = " << NT << std::endl;
	std::cout << "ResultTime = " << NT * deltaTime << std::endl;



	//SOLVE
	std::cout << "Run solution ..." << std::endl;
	double time = 0.0;
	size_t iTime = 1;
	while (time < inputData.TimeEnd())
	//for (size_t iTime = 1; iTime <= NT; iTime++)
	{
		if ((inputData.TimeEnd() - time) < deltaTime)
		{
			deltaTime = inputData.TimeEnd() - time;
		}
		time += deltaTime;
		boundCond.EvaluateBoundVelocity(time);
		mesh.RebuildMesh(boundCond ,deltaTime);

		for (size_t i = 0; i < flux.size(); i++)
		{
			NaturalVar leftVar, rightVar;
			if (tvdScheme->isOn() and (i > 0 and i < flux.size()-1))
			{
				leftVar = tvdScheme->ReconstructVar(w[i - 1].convertToNaturalVar(gas),
					w[i].convertToNaturalVar(gas),
					w[i + 1].convertToNaturalVar(gas), TVDScheme::TVDSide::left);

				rightVar = tvdScheme->ReconstructVar(w[i].convertToNaturalVar(gas),
					w[i+1].convertToNaturalVar(gas),
					w[i+2].convertToNaturalVar(gas), TVDScheme::TVDSide::right);
			}
			else
			{
				leftVar = w[i].convertToNaturalVar(gas);
				rightVar = w[i + 1].convertToNaturalVar(gas);
			}
			flux[i] = fluxScheme->CalculateFlux(leftVar, rightVar, mesh.FaceVelocity(i));

			flux[i] = flux[i] - (w[i] + w[i+1]) * 0.5 * mesh.FaceVelocity(i);
		}


		for (size_t i = 1; i < wNext.size() - 1 ; i++)
		{
			freeFlux.rho_U = 0.0;
			freeFlux.rho_U2_plus_p = w[i].convertToNaturalVar(gas).p;
			freeFlux.rho_U_H = 0.0;
			
			wNext[i] = w[i] * (mesh.OldCellVolume(i) / mesh.CellVolume(i)) - (flux[i]*mesh.S(i) - flux[i - 1]*mesh.S(i-1) - freeFlux*(mesh.S(i)- mesh.S(i-1)))
				* (deltaTime / mesh.CellVolume(i));
			
		}

		boundCond.Evaluate(wNext, time);

		if (iTime % 10 == 0)
			monPoints.UpdateMonitors(time, wNext, gas);

		std::copy(wNext.begin(), wNext.end(), w.begin());

		if (iTime % 1000 == 0)
			std::cout << (double)iTime / NT * 100 << "%" << std::endl;

		iTime++;
	}

	std::cout << "Time = " << time << " s." << std::endl;
	std::cout << "End solution" << std::endl;


	std::cout << "Write Result: ";

	try
	{
		//std::ofstream meshOut("D1_9.txt");
		//mesh.WriteMesh(meshOut);
		outputData.WriteData(inputData.ResFileName(), mesh, w, gas, naturalVariables);
		std::cout << "OK" << std::endl;;
	}
	catch (const std::exception&)
	{
		std::cout << "ERROR" << std::endl;
	}


	return 0;
}