#pragma once
#include <fstream>
#include <string>
#include "Variables.h"
#include "BoundaryConditions.h"
#include "Mesh.h"

class InputData
{
private:
	size_t _NX;
	double _length, _CFL;
	double _gamma, _Cp, _molMass;
	double _timeEnd;
	NaturalVar _left, _right;
	double _diaphragmPos;
	int  _leftBoundType, _rightBoundType;

	std::string resFileName;
	std::string bcFileName;

	int fluxScheme = 0;
	int tvdScheme = 0;

public:

	InputData(const std::string& path)
	{
		ReadInputFile(path);
	}


	void ReadInputFile(const std::string path)
	{
		std::string temp;
		std::ifstream fin(path);
		//MESH
		fin >> temp;
		fin >> temp >> _length;
		fin >> temp >> _NX;
		fin >> temp >> bcFileName;

		//SOLVER_SETTINGS
		fin >> temp;
		fin >> temp >> _timeEnd;
		fin >> temp >> _CFL;
		fin >> temp >> fluxScheme;
		fin >> temp >> tvdScheme;

		//GAS_MODEL
		fin >> temp;
		fin >> temp >> _gamma;
		fin >> temp >> _Cp;
		fin >> temp >> _molMass;

		//INITIAL_FIELD
		fin >> temp;
		fin >> temp >> _diaphragmPos;
		fin >> temp >> _left.p;
		fin >> temp >> _left.rho;
		fin >> temp >> _left.u;
		fin >> temp >> _right.p;
		fin >> temp >> _right.rho;
		fin >> temp >> _right.u;

		//RESULT_FILE_NAME
		fin >> temp;
		fin >> temp >> resFileName;



	}


	void Write(std::ostream& out)
	{
		out << "Input Data :\n";
		out << "length = " << _length << "\n";
		out << "NX = " << _NX << "\n";
		out << "CFL = " << _CFL << "\n";
		out << "End Time = " << _timeEnd << "\n\n";

	}
	double TimeEnd() const
	{
		return _timeEnd;
	}
	size_t NX() const
	{
		return _NX;
	}
	double CFL() const
	{
		return _CFL;
	}
	double Length() const
	{
		return _length;
	}
	double gamma() const
	{
		return _gamma;
	}
	double Cp()
	{
		return _Cp;
	}
	double molMass()
	{
		return _molMass;
	}
	double diaphragmPos() const
	{
		return _diaphragmPos;
	}
	int FluxScheme() const
	{
		return fluxScheme;
	}
	int TVDScheme() const
	{
		return tvdScheme;
	}
	std::string ResFileName() const
	{
		return resFileName;
	}
	std::string BCFileName() const
	{
		return bcFileName;
	}



	enum Type { left, right };
	NaturalVar initialField(const Type type) const
	{
		switch (type)
		{
		case left:
			return _left;
			break;
		case right:
			return _right;
			break;

		default:
			break;
		}
	}

	BoundType boundType(const BoundSide type) const
	{
		switch (type)
		{
		case BoundSide::left:
			return BoundType(_leftBoundType);
			break;
		case BoundSide::right:
			return BoundType(_rightBoundType);
			break;

		default:
			break;
		}
	}
};


enum typeOutputData
{
	naturalVariables,
	ConservativeVariables
};
class OutputData
{


public:
	void WriteData(const std::string& path, const Mesh& mesh,
		const std::vector<ConservativeVar>& w,
		const GasModel& gas, typeOutputData type)
	{
		std::ofstream fout(path);

		switch (type)
		{
		case naturalVariables:
			fout << "\'VARIABLES = \"position\",\"density\",\"velocity\",\"pressure\", \"energy\"\'" << std::endl;
			break;
		case ConservativeVariables:
			fout << "\'VARIABLES = \"position\",\"density\",\"density*velocity\",\"density*energy\"\'" << std::endl;
			break;
		default:
			break;
		}

		for (size_t i = 1; i < mesh.GetNX(); i++)
		{
			switch (type)
			{
			case naturalVariables:
			{
				NaturalVar natVar = w[i].convertToNaturalVar(gas);
				fout << mesh.CellCenter(i) << "\t" << natVar << "\t" << natVar.InternalEnergy(gas) << "\n";
				break;
			}
			case ConservativeVariables:
				fout << mesh.CellCenter(i) << "\t" << w[i] << "\n";
				break;
			default:
				break;
			}
		}
	}

};
