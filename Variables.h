#pragma once
#include "GasModel.h"
#include <iostream>

class ConservativeVar;
struct LimiterArg;

struct Limiter
{
	double limitRho;
	double limitP;
	double limitU;

	static Limiter GetMin(const Limiter& var1, const Limiter& var2)
	{

		return Limiter{ var1.limitRho < var2.limitRho ? var1.limitRho : var2.limitRho ,
								var1.limitP < var2.limitP ? var1.limitP : var2.limitP,
									var1.limitU < var2.limitU ? var1.limitU : var2.limitU };

	}

	static Limiter GetMax(const Limiter& var1, const Limiter& var2)
	{

		return Limiter{ var1.limitRho > var2.limitRho ? var1.limitRho : var2.limitRho ,
								var1.limitP > var2.limitP ? var1.limitP : var2.limitP,
									var1.limitU > var2.limitU ? var1.limitU : var2.limitU };

	}

};

struct NaturalVar
{
	double p;
	double rho;
	double u;

	NaturalVar() { p = 0.0; rho = 0.0; u = 0.0; };
	NaturalVar(const double rho, const double p, const double u)
		:rho(rho), p(p), u(u) {};

	NaturalVar operator- (const NaturalVar& other) const
	{
		return NaturalVar(rho - other.rho, p - other.p, u - other.u);
	}

	NaturalVar operator+ (const NaturalVar& other) const
	{
		return NaturalVar(rho + other.rho, p + other.p, u + other.u);
	}

	NaturalVar operator* (const double val) const
	{
		return NaturalVar(rho * val, p * val, u * val);
	}

	NaturalVar operator/ (const NaturalVar& other) const
	{
		return NaturalVar(rho / other.rho, p / other.p, u / other.u);
	}


	ConservativeVar convertToConsVar(const GasModel& gas);

	NaturalVar operator* (const Limiter& other) const
	{
		return NaturalVar(rho * other.limitRho, p * other.limitP, u * other.limitU);
	}

	double InternalEnergy(const GasModel& gas) const
	{
		return gas.T(p, rho) * gas.Cv();
	}

};


struct LimiterArg
{
	double argRho;
	double argP;
	double argU;

	void CalculateArg(const NaturalVar& left, const NaturalVar& center, const NaturalVar& right, bool reverse = false)
	{
		NaturalVar lower;
		NaturalVar temp;
		if (reverse)
		{
			lower = right - center;
			temp = (center - left)/ lower ;
		}
		else
		{
			lower = center - left;
			temp = (right - center) / lower;
		}
		*this = { lower.rho == 0 || temp.rho <= 0 ? 0 : temp.rho,
			      lower.p == 0   || temp.p <= 0 ? 0 : temp.p,
				  lower.u == 0   || temp.u <= 0 ? 0 : temp.u };


	}

	LimiterArg operator+ (const double& other) const
	{
		return LimiterArg{ argRho + other, argP + other, argU + other };
	}

	LimiterArg operator* (const double& other) const
	{
		return LimiterArg{ argRho * other, argP * other, argU * other };
	}

	LimiterArg operator / (const LimiterArg& other) const
	{
		return LimiterArg{ argRho / other.argRho, argP  / other.argP, argU / other.argU };
	}

	LimiterArg operator * (const LimiterArg& other) const
	{
		return LimiterArg{ argRho * other.argRho, argP * other.argP, argU * other.argU};
	}

	LimiterArg operator + (const LimiterArg& other) const
	{
		return LimiterArg{ argRho + other.argRho, argP + other.argP, argU + other.argU };
	}

	operator Limiter () const
	{
		return Limiter{ this->argRho, this->argP, this->argU };
	}

};




class vectorFlux
{

private:
	friend class ConservativeVar;
	friend std::ostream& operator<< (std::ostream&, const vectorFlux&);
public:
	double rho_U{ 0 };
	double rho_U2_plus_p{ 0 };
	double rho_U_H{ 0 };
public:
	vectorFlux() {};
	vectorFlux(const double rho_U, const double rho_U2_plus_p, const double rho_U_H) :
		rho_U(rho_U), rho_U2_plus_p(rho_U2_plus_p), rho_U_H(rho_U_H) {}

	vectorFlux operator* (double value)
	{
		return vectorFlux(rho_U * value, rho_U2_plus_p * value, rho_U_H * value);
	}
	vectorFlux operator/ (double value)
	{
		return vectorFlux(rho_U / value, rho_U2_plus_p / value, rho_U_H / value);
	}

	vectorFlux operator- (const vectorFlux& other)
	{
		return vectorFlux(rho_U - other.rho_U, rho_U2_plus_p - other.rho_U2_plus_p, rho_U_H - other.rho_U_H);
	}

	vectorFlux operator- (const ConservativeVar& other);


	vectorFlux operator+ (const vectorFlux& other)
	{
		return vectorFlux(rho_U + other.rho_U, rho_U2_plus_p + other.rho_U2_plus_p, rho_U_H + other.rho_U_H);
	}

	vectorFlux& operator= (const vectorFlux& other)
	{
		if (&other != this)
		{
			this->rho_U = other.rho_U;
			this->rho_U2_plus_p = other.rho_U2_plus_p;
			this->rho_U_H = other.rho_U_H;
		}
		return *this;
	}
	vectorFlux(const vectorFlux& other)
	{
		this->rho_U = other.rho_U;
		this->rho_U2_plus_p = other.rho_U2_plus_p;
		this->rho_U_H = other.rho_U_H;
	}

};

class ConservativeVar
{
private:
	double rho{0};
	double rho_U{ 0 };
	double rho_E{ 0 };
	friend std::ostream& operator<< (std::ostream&, const ConservativeVar&);
public:
	ConservativeVar() {};
	ConservativeVar(double rho, double rho_U, double rho_E) : rho(rho), rho_U(rho_U), rho_E(rho_E) {}


	ConservativeVar(const ConservativeVar& other)
	{
		rho = other.rho;
		rho_U = other.rho_U;
		rho_E = other.rho_E;
	}

	ConservativeVar& operator= (const ConservativeVar& other)
	{
		if (&other != this)
		{
			this->rho = other.rho;
			this->rho_U = other.rho_U;
			this->rho_E = other.rho_E;
		}
		return *this;
	}

	ConservativeVar operator+ (const ConservativeVar& other) const
	{
		return ConservativeVar{ rho + other.rho, rho_U + other.rho_U, rho_E + other.rho_E };
	}

	ConservativeVar operator+ (const vectorFlux& flux) const
	{
		return ConservativeVar{ rho + flux.rho_U, rho_U + flux.rho_U2_plus_p, rho_E + flux.rho_U_H };
	}
	ConservativeVar operator- (const vectorFlux& flux) const
	{
		return ConservativeVar{ rho - flux.rho_U, rho_U - flux.rho_U2_plus_p, rho_E - flux.rho_U_H };
	}

	ConservativeVar operator* (const double& value) const
	{
		return ConservativeVar{ rho * value, rho_U * value, rho_E * value };
	}

	void ReverseVelocity()
	{
		rho_U = -rho_U;
	}

	NaturalVar convertToNaturalVar(const GasModel& gas) const
	{
		NaturalVar result;
		result.rho = rho;
		result.u = rho_U / rho;
		double T = (rho_E / rho - result.u * result.u / 2.0) / gas.Cv();
		result.p = gas.P(T, result.rho);
		return result;
	}

	double Rho() const
	{
		return rho;
	}

	double Rho_U() const
	{
		return rho_U;
	}

	double Rho_E() const
	{
		return rho_E;
	}

};


vectorFlux vectorFlux::operator-(const ConservativeVar& other)
{
	return vectorFlux(rho_U - other.Rho(), rho_U2_plus_p - other.Rho_U(), rho_U_H - other.Rho_E());
}

ConservativeVar NaturalVar::convertToConsVar(const GasModel& gas)
{
	double E = gas.Cv() * gas.T(p, rho) + u * u / 2.0;
	return ConservativeVar(rho, rho * u, rho * E);
}

std::ostream& operator<< (std::ostream& out, const vectorFlux& vector)
{
	out << vector.rho_U << "\t" << vector.rho_U2_plus_p << "\t" << vector.rho_U_H;
	return out;
}

std::ostream& operator<< (std::ostream& out, const ConservativeVar& vector)
{
	out << vector.rho << "\t" << vector.rho_U << "\t" << vector.rho_E;
	return out;
}

std::ostream& operator<< (std::ostream& out, const NaturalVar& vector)
{
	out << vector.rho << "\t" << vector.u << "\t" << vector.p;
	return out;
}