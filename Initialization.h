#pragma once
#include <vector>
#include "Variables.h"
#include "GasModel.h"
#include "Functions.h"
#include "Mesh.h"
class Initialization
{
public:
	void initialize(const Mesh& mesh, std::vector <ConservativeVar>& w,
		const GasModel& gasModel,
		const NaturalVar& left, const NaturalVar& right,
		const double diaphragmPos)
	{
		for (size_t i = 0; i < w.size(); i++)
		{
			NaturalVar temp;
			if (mesh.CellCenter(i) < diaphragmPos)
			{
				temp.rho = left.rho;
				temp.p = left.p;
				temp.u = left.u;
			}
			else
			{
				temp.rho = right.rho;
				temp.p = right.p;
				temp.u = right.u;
			}
			w[i] = temp.convertToConsVar(gasModel);
		}
	}


};
