#pragma once
#include <vector>
#include <fstream>
#include "Mesh.h"
#include "Variables.h"
#include "GasModel.h"
#include <string>


class MonitorPoints
{
private:
	struct index
	{
		size_t left;
		size_t right;
	};
	std::vector<index> indexes;
	std::vector<double> positions;
	std::vector<std::ofstream> fstreams;
	std::vector<std::string> fileNames;
	size_t size;
	Mesh* mesh;
public:
	MonitorPoints(Mesh& mesh): mesh(&mesh)
	{
		double x;
		std::string fileName;
		std::ifstream monFin("monitorPoints.mon");
		while (!monFin.eof())
		{
			monFin >> x >> fileName;
			positions.push_back(x);
			fileNames.push_back(fileName);
		}

		size = positions.size();
		indexes.resize(size);
		fstreams.resize(size);

		CalculateIndexes();
		OpenWriteFiles();
	}

	~MonitorPoints()
	{
		CloseWriteFiles();
	}



	void UpdateMonitors(const double  flow_time, const std::vector <ConservativeVar> & var, const GasModel& gas)
	{
		CalculateIndexes();
		for (size_t i = 0; i < size; i++)
		{
			NaturalVar varL = var[indexes[i].left].convertToNaturalVar(gas);
			NaturalVar varR = var[indexes[i].right].convertToNaturalVar(gas);

			double leftLen = abs(mesh->CellCenter(indexes[i].left) - positions[i]);
			double rightLen = abs(mesh->CellCenter(indexes[i].right) - positions[i]);


			double u = LinearInterpolation(varL.u, varR.u, leftLen, rightLen);
			double p = LinearInterpolation(varL.p, varR.p, leftLen, rightLen);
			double rho = LinearInterpolation(varL.rho, varR.rho, leftLen, rightLen);
			double e = LinearInterpolation(varL.InternalEnergy(gas), varR.InternalEnergy(gas), leftLen, rightLen);

			fstreams[i] << flow_time << "\t" << rho << "\t" << u << "\t" << p << "\t" << e << std::endl;
		}
	
	}


private:
	void CalculateIndexes()
	{
		for (size_t i = 0; i < size; i++)
		{
			if (positions[i] >= mesh->Node(0) && positions[i] <= mesh->Node(mesh->GetNX()))
			{
				for (size_t j = 0; j < mesh->GetNX() + 1; j++)
				{
					if (positions[i] >= mesh->CellCenter(j) && positions[i] <= mesh->CellCenter(j+1))
					{
						indexes[i].left = j;
						indexes[i].right = j + 1;
						break;
					}
				}
			}
		}
		
	}

	void OpenWriteFiles()
	{
		for (size_t i = 0; i < size; i++)
		{
			fstreams[i].open(fileNames[i]);
			fstreams[i] << "\'VARIABLES = \"flow-time\",\"rho\",\"u\",\"p\",\"e\"\'" << std::endl;
		}
	}


	void CloseWriteFiles()
	{
		for (size_t i = 0; i < size; i++)
		{
			fstreams[i].close();
		}
	}

	/*double Average(const double l, const double r) const 
	{
		return (l + r) / 2.0;
	}*/

	double LinearInterpolation(const double varL, const double varR, 
		const double leftLen, const double rightLen) const
	{
		return (varL * rightLen + varR * leftLen) / (leftLen + rightLen);
	}



};