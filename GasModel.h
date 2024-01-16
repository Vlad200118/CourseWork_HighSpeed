#pragma once
#include <math.h>
class GasModel
{
private:
	double _Cp;
	double _Cv;
	double _gamma;
	double _molMass;
	double R_M;
	const double R = 8.31446262;
public:
	GasModel() = delete;
	GasModel(const double Cp, const double gamma, const double molMass) :
		_Cp(Cp), _gamma(gamma), _molMass(molMass)
	{
		R_M = R / _molMass;
		_Cv = _Cp - R_M;
	}
	double T(const double pressure, const double density) const
	{
		return pressure / density / R_M;
	}
	double P(const double temperature, const double density) const
	{
		return temperature * density * R_M;
	}
	double Rho(const double pressure, const double temperature) const
	{
		return pressure / temperature / R_M;
	}

	double gamma() const
	{
		return _gamma;
	}
	double Cp() const
	{
		return _Cp;
	}
	double Cv() const
	{
		return _Cv;
	}
	double SpeedSound(const double pressure, const double density) const
	{
		return sqrt(_gamma * pressure / density);
	}
	double MolWieght() const
	{
		return _molMass;
	}
};
