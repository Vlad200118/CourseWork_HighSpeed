#pragma once
#include "Variables.h"
#include "GasModel.h"
#include "InputOutputData.h"
#include <iostream>
#include <fstream>
#include <vector>


void WriteData(const std::string& path, std::vector<double>& x, std::vector<ConservativeVar>& w)
{
	std::ofstream fout(path);

	for (size_t i = 1; i < x.size() - 1; i++)
	{
		fout << x[i] << "\t" << w[i] << "\n";
	}
}

inline double GetVelocityForCFL(const InputData& input, const GasModel& gas)
{
	double maxSpeedSound = std::max(gas.SpeedSound(input.initialField(InputData::left).p, input.initialField(InputData::left).rho),
		gas.SpeedSound(input.initialField(InputData::right).p, input.initialField(InputData::right).rho));
	double velocity = 10.0;
	while ((maxSpeedSound / velocity) > 1.0)
	{
		velocity += 10.0;
	}
	return velocity;
	//return maxSpeedSound;
}


template<class T>
void ResizeVectors(const size_t& size, std::vector<T>& vector)
{
	vector.resize(size);
};
template<class T, class ... Args>
void ResizeVectors(const size_t& size, std::vector<T>& vector, Args&... args)
{
	vector.resize(size);
	ResizeVectors(size, args...);
};