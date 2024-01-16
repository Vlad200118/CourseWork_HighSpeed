#pragma once
#include <vector>
#include "Variables.h"
#include <map>
#include <fstream>
//#include "Functions.h"

class BoundaryConditions;

enum BoundType { wall, open, plunger, symmetry };
enum BoundSide { left = 1, right = -1 };

class BoundaryType
{
private:
	size_t borderInd;
	size_t overseasInd;
	BoundSide side;
	GasModel* gas;
	double boundaryVelocity = 0.0;
	friend class Wall;
	friend class Open;
	friend class Plunger;
	friend class Symmetry;

public:
	BoundaryType(const size_t& meshSize, BoundSide side, GasModel& gas)
	{
		this->gas = &gas;
		this->side = side;
		if (side == BoundSide::left)
		{
			borderInd = 1;
			overseasInd = 0;
		}
		else
		{
			borderInd = meshSize;
			overseasInd = meshSize + 1;
		}

	}
	virtual void Evaluate(std::vector<ConservativeVar>& vector, const double time) const = 0;

	double BoundaryVelocity()
	{
		return boundaryVelocity;
	}

	virtual void EvaluateBoundVelocity(const double time) = 0;
};


class Wall : public BoundaryType
{
public:
	Wall(const size_t& meshSize, BoundSide side, GasModel& gas, const double wallVelocity = 0.0) : BoundaryType(meshSize, side, gas)
	{
		boundaryVelocity = wallVelocity;
	}

	void Evaluate(std::vector<ConservativeVar>& vector, const double time = 0.0) const
	{
		NaturalVar borderVar = vector[borderInd].convertToNaturalVar(*gas);
		NaturalVar overseasVar;

		overseasVar.p = borderVar.p;
		overseasVar.rho = borderVar.rho;
		overseasVar.u = 2.0 * boundaryVelocity - borderVar.u;
		vector[overseasInd] = overseasVar.convertToConsVar(*gas);
	}

	void SetWallVelocity(const double& newValue)
	{
		boundaryVelocity = newValue;
	}

	void EvaluateBoundVelocity(const double time) override
	{
		//std::cout << "WALL = " << boundaryVelocity;
	}

};

class Open : public BoundaryType
{
private:
	double backflowPressure;
	double backflowTemperature;
public:
	Open(const size_t& meshSize, BoundSide side, GasModel& gas, double BFP, double BFT) :
		BoundaryType(meshSize, side, gas)
	{
		backflowPressure = BFP;
		backflowTemperature = BFT;
		BoundaryType::boundaryVelocity = 0.0;
	}
	void Evaluate(std::vector<ConservativeVar>& vector, const double time = 0.0) const
	{
		NaturalVar bordeeVar = vector[borderInd].convertToNaturalVar(*gas);
		NaturalVar overseasVar;
		double speedSound = gas->SpeedSound(bordeeVar.p, bordeeVar.rho);


		if (bordeeVar.u * side <= 0.0)
		{
			//Outlet
			overseasVar.rho = bordeeVar.rho;
			overseasVar.u = bordeeVar.u;
			if (abs(bordeeVar.u) < speedSound)
			{
				overseasVar.p = backflowPressure;
			}
			else
			{
				overseasVar.p = bordeeVar.p;
			}
		}
		else
		{
			//Inlet
			double mach = bordeeVar.u / speedSound;
			double backflowDensity = gas->Rho(backflowPressure, backflowTemperature);

			double temp = (1.0 + mach * mach * (gas->gamma() - 1.0) / 2.0);
			double powerDensity = -1.0 / (gas->gamma() - 1.0);
			double powerPressure = powerDensity * gas->gamma();

			overseasVar.rho = backflowDensity * std::pow(temp, powerDensity);
			overseasVar.p = backflowPressure * std::pow(temp, powerPressure);

			if (abs(bordeeVar.u) < speedSound)
			{
				overseasVar.u = bordeeVar.u;
			}
			else
			{
				overseasVar.u = mach * gas->SpeedSound(backflowPressure, backflowDensity);
			}
		}

		vector[overseasInd] = overseasVar.convertToConsVar(*gas);
	}


	void EvaluateBoundVelocity(const double time) override {}

};


class Plunger : public Wall
{
private:
	double frequency = 0.0;
	double amplitude = 0.0;

public:
	Plunger(const size_t& meshSize, BoundSide side, GasModel& gas, const double freq, const double ampl) :
		Wall(meshSize, side, gas)
	{
		frequency = freq;
		amplitude = ampl;
	}


	void Evaluate(std::vector<ConservativeVar>& vector, const double time)
	{
		Wall::Evaluate(vector);
	}

	void EvaluateBoundVelocity(const double time) override
	{
		boundaryVelocity = frequency * amplitude * sin(frequency * time);
	}


};


class Symmetry : public BoundaryType
{
public:
	Symmetry(const size_t& meshSize, BoundSide side, GasModel& gas) : BoundaryType(meshSize, side, gas)
	{
		boundaryVelocity = 0.0;
	}

	void Evaluate(std::vector<ConservativeVar>& vector, const double time = 0.0) const
	{
		NaturalVar borderVar = vector[borderInd].convertToNaturalVar(*gas);
		NaturalVar overseasVar;

		overseasVar.p = borderVar.p;
		overseasVar.rho = borderVar.rho;
		overseasVar.u = borderVar.u;
		vector[overseasInd] = overseasVar.convertToConsVar(*gas);
	}

	void EvaluateBoundVelocity(const double time) override
	{
		//std::cout << "SYMMETRY = " << boundaryVelocity;
	}


};

class BoundaryConditions
{
private:
	std::map<std::string, BoundaryType*> bounds;

public:

	BoundaryConditions(const size_t meshSize, GasModel& gas, const std::string& path)
	{
		ReadBCFileName(meshSize, gas, path);
	}

	void Evaluate(std::vector<ConservativeVar>& vector, const double time)
	{
		for (auto& boundary : bounds)
		{
			boundary.second->Evaluate(vector, time);
		}
	}

	auto GetBoundary(const std::string name) 
	{
		return bounds[name];
	}

	void EvaluateBoundVelocity(const double time)
	{
		for (auto& boundary : bounds)
		{
			boundary.second->EvaluateBoundVelocity(time);
		}
	}

	~BoundaryConditions()
	{
		for (auto& boundary : bounds)
		{
			delete boundary.second;
		}
	}

private:
	void ReadBCFileName(const size_t meshSize, GasModel& gas, const std::string& path)
	{
		std::string boundSide;
		std::string line;
		int boundType;
		BoundSide side;

		std::ifstream fin(path);
		for (size_t i = 0; i < 2; i++)
		{
			fin >> boundSide >> boundType;
			side = boundSide == "left" ? BoundSide::left : BoundSide::right;
			if (boundType == BoundType::wall)
			{
				double wallVelocity;
				fin >> wallVelocity;
				bounds[boundSide] = new Wall(meshSize, side, gas, wallVelocity);
			}
			else if (boundType == BoundType::open)
			{
				double backflowPressure, backflowTemperature;
				fin >> backflowPressure >> backflowTemperature;
				bounds[boundSide] = new Open(meshSize, side,gas, backflowPressure, backflowTemperature);

			}
			else if (boundType == BoundType::plunger)
			{
				double freq, ampl;
				fin >> freq >> ampl;
				bounds[boundSide] = new Plunger(meshSize, side, gas, freq, ampl);
			}
			else if (boundType == BoundType::symmetry)
			{
				bounds[boundSide] = new Symmetry(meshSize, side, gas);
			}
		}

	}


};
