#pragma once
#include "Variables.h"
#include "GasModel.h"

enum splitSign { minus = -1, plus = 1 };

class FluxScheme
{
public:
	GasModel* gas;

	FluxScheme(GasModel& gas) :gas(&gas) {};
	virtual vectorFlux CalculateFlux(const NaturalVar& varLeft, const NaturalVar& varRight, double faceVelocity=0.0) = 0;
	~FluxScheme()
	{
		gas = nullptr;
	}

};


class AUSM_Standard : public FluxScheme
{
private:
	enum splitSign { minus = -1, plus = 1 };
	double machLeft{0}, machRight{ 0 }, machFace{0};
	double soundSpeedLeft{0}, soundSpeedRight{ 0 };
public:
	AUSM_Standard(GasModel& gas) : FluxScheme(gas) {}

	vectorFlux CalculateFlux(const NaturalVar& varLeft, const NaturalVar& varRight, double faceVelocity=0.0) override
	{
		CalculateSoundSpeed(varLeft, varRight);
		CalculateMachNumbers(varLeft, varRight);
		CalculateMachFace();
		
		vectorFlux fluxPressure = CalculatePressureFluxFace(varLeft, varRight);


		vectorFlux fluxContinuity = CalculateContinuityFluxFace(varLeft, varRight);


		return fluxContinuity + fluxPressure;
	}

private:


	double CalculateMachSplit(const double& mach, splitSign sign) const
	{

		if (abs(mach) <= 1.0)
		{
			return sign * (mach + sign) * (mach + sign) / 4.0;
		}
		else
		{
			return (mach + sign * abs(mach)) / 2.0;
		}

	}

	double CalculatePressureSplit(const double mach, splitSign sign) const
	{

		if (abs(mach) >= 1.0)
		{
			return (mach + sign * abs(mach)) / mach * 0.5;
		}
		else
		{
			return (mach + sign) * (mach + sign) * (2.0 - sign * mach) * 0.25;
		}
	}

	vectorFlux CalculateContinuityFlux(const NaturalVar& var) const 
	{
		double H = gas->Cp() * gas->T(var.p, var.rho) + var.u * var.u / 2.0;
		return vectorFlux(var.rho, var.rho * var.u, var.rho * H);
	}

	vectorFlux CalculateContinuityFluxFace(const NaturalVar& varLeft, const NaturalVar& varRight) const
	{
		if (machFace >= 0.0)
		{
			return CalculateContinuityFlux(varLeft) * machFace * soundSpeedLeft;
		}
		else
		{
			return CalculateContinuityFlux(varRight) * machFace * soundSpeedRight;
		}
	}

	vectorFlux CalculatePressureFluxFace(const NaturalVar& varLeft, const NaturalVar& varRight) const 
	{
		double facePressure = CalculatePressureSplit(machLeft, plus) * varLeft.p +
			CalculatePressureSplit(machRight, minus) * varRight.p;
		return vectorFlux(0.0, facePressure, 0.0);
	}

	void CalculateMachNumbers(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		machLeft = varLeft.u / soundSpeedLeft;
		machRight = varRight.u / soundSpeedRight;
	}

	void CalculateSoundSpeed(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		soundSpeedLeft = gas->SpeedSound(varLeft.p, varLeft.rho);
		soundSpeedRight = gas->SpeedSound(varRight.p, varRight.rho);
	}

	void CalculateMachFace()
	{
		machFace = CalculateMachSplit(machLeft, plus) + CalculateMachSplit(machRight, minus);
	}
};

class AUSM_Plus : public FluxScheme
{
private:
	double soundSpeedFace = 0.0;
	const double alpha = 3.0 / 16.0;
	const double beta = 1.0 / 8.0;
	double machLeft{ 0 }, machRight{ 0 }, machFace{0};
public:
	AUSM_Plus(GasModel& gas) : FluxScheme(gas) {}

	vectorFlux CalculateFlux(const NaturalVar& varLeft, const NaturalVar& varRight, double faceVelocity = 0.0) override
	{
		CalculateSoundSpeed(varLeft, varRight);
		CalculateMachNumbers(varLeft, varRight);
		CalculateMachFace();
		

		vectorFlux fluxContinuity = CalculateContinuityFluxFace(varLeft, varRight);
		vectorFlux fluxPressure = CalculatePressureFluxFace(varLeft, varRight);

		return fluxContinuity + fluxPressure;
	}

private:

	double CalculateMachSplit(const double& mach, const splitSign sign)
	{

		if (abs(mach) >= 1.0)
		{
			return (mach + sign * abs(mach)) * 0.5; 
		}
		else
		{
			return sign * (mach + sign) * (mach + sign) * 0.25 + sign * beta * (mach*mach - 1) * (mach * mach - 1);
		}

	}

	double CalculatePressureSplit(const double& mach, const splitSign sign)
	{

		if (abs(mach) >= 1.0)
		{
			return  (mach + sign * abs(mach)) / mach * 0.5;
		}
		else
		{
			return  (mach + sign) * (mach + sign) * (2.0 - sign * mach) / 4.0 + 
				alpha*sign * mach * (mach*mach - 1.0) * (mach * mach - 1.0);
		}
	}

	double CalculateMachFaceSplit(const splitSign sign)
	{
		return (machFace + sign * abs(machFace))*0.5;
	}

	void CalculateMachFace()
	{
		machFace =  CalculateMachSplit(machLeft, plus) + CalculateMachSplit(machRight, minus);
	}

	vectorFlux CalculateContinuityFlux(const NaturalVar& var)
	{
		double H = gas->Cp() * gas->T(var.p, var.rho) + var.u * var.u / 2.0;
		return vectorFlux(var.rho , var.rho * var.u, var.rho * H);
	}

	vectorFlux CalculateContinuityFluxFace(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		return (CalculateContinuityFlux(varLeft) * CalculateMachFaceSplit(plus) +
			CalculateContinuityFlux(varRight) * CalculateMachFaceSplit(minus))* soundSpeedFace;
	}

	vectorFlux CalculatePressureFluxFace(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		double facePressure = CalculatePressureSplit(machLeft, plus) * varLeft.p +
			CalculatePressureSplit(machRight, minus) * varRight.p;
		return vectorFlux(0.0, facePressure, 0.0);
	}

	void CalculateSoundSpeed(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		//double speedStar = sqrt(abs(varLeft.u * varRight.u));
		//double gamma = gas->gamma();
		//GasModel* gas = this->gas;
		////double speedStar =sqrt( 2.0 * (gamma - 1) / (gamma + 1) * (gas->Cp() * gas->T(varLeft.p, varLeft.rho) + varLeft.u * varLeft.u / 2.0));
		////double sp2 = sqrt(2.0 * (gamma - 1) / (gamma + 1) * (gas->Cp() * gas->T(varRight.p, varRight.rho) + varRight.u * varRight.u / 2.0));
		//
		//
		//auto speedStar = [&gamma, &gas](const NaturalVar& var)
		//{
		//	return sqrt(2.0 * (gamma - 1) / (gamma + 1) * (gas->Cp() * gas->T(var.p, var.rho) + var.u * var.u / 2.0));
		//};

		//auto speedTilda = [&speedStar](const NaturalVar& var)
		//{
		//	double s = speedStar(var);
		//	return s * s / std::max(s, abs(var.u));
		//};

		////soundSpeedFace = std::min(speedTilda(varLeft), speedTilda(varRight));
	
		////soundSpeedFace = sqrt (gas->SpeedSound(varLeft.p, varLeft.rho) * gas->SpeedSound(varRight.p, varRight.rho));
		soundSpeedFace = (gas->SpeedSound(varLeft.p, varLeft.rho) + gas->SpeedSound(varRight.p, varRight.rho))/2.0;

	}

	void CalculateMachNumbers(const NaturalVar& varLeft, const NaturalVar& varRight)
	{
		machLeft = varLeft.u / soundSpeedFace;
		machRight = varRight.u / soundSpeedFace;
	}
};

class Godunov : public FluxScheme
{
private:
	double PBIG{0}, UBIG{ 0 }, UDOT{ 0 };
	double RBIG1{ 0 }, RBIG2{ 0 };
	double DL1{ 0 }, DL2{ 0 };
	double DP1{ 0 }, DP2{ 0 };

	NaturalVar varFace;

public:

	Godunov(GasModel& gas) : FluxScheme(gas) {}

	vectorFlux CalculateFlux(const NaturalVar& varLeft, const NaturalVar& varRight, double faceVelocity = 0.0) override
	{
		UDOT = faceVelocity;
		RASPAD(varLeft, varRight);
		POTOK(varLeft, varRight);

		double flux1 = varFace.rho * varFace.u;
		double flux2 = varFace.rho * varFace.u * varFace.u + varFace.p;
		double H = gas->Cp() * gas->T(varFace.p, varFace.rho) + (varFace.u * varFace.u) / 2.0;
		double flux3 = varFace.rho * varFace.u * H;

		return vectorFlux(flux1, flux2, flux3);
	}
private:

	void POTOK(const NaturalVar& left, const NaturalVar& right)
	{
		double AK1, AK2;
		double CZJM, CZJP, CZT;
		AK1 = gas->gamma();
		AK2 = gas->gamma();
		varFace.u = UBIG;
		varFace.p = PBIG;

		if (UBIG - UDOT > 0.0)
		{
			if ((DL1 - DL2) == 0.0)
			{
				if (DL1 - UDOT >= 0.0)
				{
					varFace.rho = left.rho;
					varFace.u = left.u;
					varFace.p = left.p;
				}
				else
				{
					varFace.rho = RBIG1;
				}
			}
			else
			{
				if (DL2 - UDOT <= 0.0)
					varFace.rho = RBIG1;
				else
				{
					if (DL1 - UDOT >= 0.0)
					{
						varFace.rho = left.rho;
						varFace.u = left.u;
						varFace.p = left.p;
					}
					else
					{
						CZJM = sqrt(AK1 * left.p / left.rho);
						CZT = ((AK1 - 1.) * (left.u - UDOT) + 2. * CZJM) / (AK1 + 1.);
						varFace.u = CZT + UDOT;
						varFace.p = left.p * pow(CZT / CZJM, 2. * AK1 / (AK1 - 1.));
						varFace.rho = AK1 * varFace.p / (CZT * CZT);
					}
				}
			}
		}
		else
		{
			if (DP1 - DP2 == 0.0)
			{
				if (DP1 - UDOT <= 0.0)
				{
					varFace.rho = right.rho;
					varFace.u = right.u;
					varFace.p = right.p;
				}
				else
					varFace.rho = RBIG2;
			}
			else
			{
				if (DP1 - UDOT >= 0.0)
					varFace.rho = RBIG2;
				else
				{
					if (DP2 - UDOT <= 0.0)
					{
						varFace.rho = right.rho;
						varFace.u = right.u;
						varFace.p = right.p;
					}
					else
					{
						CZJP = sqrt(AK2 * right.p / right.rho);
						CZT = (-(AK2 - 1.) * (right.u - UDOT) + 2. * CZJP) / (AK2 + 1.);
						varFace.u = -CZT + UDOT;
						varFace.p = right.p * pow(CZT / CZJP, 2. * AK2 / (AK2 - 1.));
						varFace.rho = AK2 * varFace.p / (CZT * CZT);
					}
				}
			}
		}
		return;
	}


	double FGOD(const double P, const double R, const double AK)
	{
		double PI = PBIG / P;
		double ST = 0.5 * (AK - 1.0) / AK;

		if (PI >= 1.0)
			return sqrt(P / R) * (PI - 1.0) / sqrt(0.5 * (AK + 1.0) * PI + 0.5 * (AK - 1.0));
		else
			return 2.0 / (AK - 1.0) * sqrt(AK * P / R) * (pow(PI, ST) - 1.0);

	}

	double DIVFGOD(const double P, const double R, const double AK)
	{
		double PI = PBIG / P;
		double ST = 0.5 * (AK - 1.0) / AK;

		if (PI >= 1.0)

			return ((AK + 1.0) * PI + (3.0 * AK - 1.0))
			/ (4.0 * AK * R * sqrt(AK * P / R) *
				sqrt(pow(0.5 * (AK + 1.0) / AK * PI + 0.5 * (AK - 1.0) / AK, 3)));
		else

			return 1.0 / AK / PBIG * sqrt(AK * P / R) * pow(PI, ST);
	}

	void RASPAD(const NaturalVar& left, const NaturalVar& right)
	{
		bool FLAG;
		double AK1, AK2, C1, C2, ST1, ST2, A1, A2, C2Z, C1Z, PI, EPSP, EPS;
		double PBIGN;

		EPSP = 1.0E-8;
		EPS = 1.0E-8;

		AK1 = gas->gamma();
		AK2 = gas->gamma();
		C1 = sqrt(AK1 * left.p / left.rho);
		C2 = sqrt(AK2 * right.p / right.rho);

		if ((left.u - right.u + 2.0 * (C1 / (AK1 - 1.0) + C2 / (AK2 - 1.0))) < 0.0)
		{
			PBIG = 1.0e-05;
			RBIG1 = 0.0;
			RBIG2 = 0.0;
			return;
		}

		FLAG = false;
		PBIG = (left.p * right.rho * C2 + right.p * left.rho * C1 + (left.u - right.u) * left.rho * C1 * right.rho * C2) /
			(left.rho * C1 + right.rho * C2);

		if (PBIG < 0.0)
			PBIG = 1.0E-05;

	CalculatePBIG:
		PBIGN = PBIG - (FGOD(left.p, left.rho, AK1) + FGOD(right.p, right.rho, AK2) - (left.u - right.u)) /
			(DIVFGOD(left.p, left.rho, AK1) + DIVFGOD(right.p, right.rho, AK2));

		if (PBIGN < 0.0)
		{
			std::cout << PBIG << "\t" << FLAG << "\n";
			PBIG = 1.0E-05;
			RBIG1 = 0.0;
			RBIG2 = 0.0;
			if (FLAG) return;
			FLAG = true;
			goto CalculatePBIG;
		}
		if (abs(PBIGN / PBIG - 1.0) > EPS)
		{
			PBIG = PBIGN;
			goto CalculatePBIG;
		}

		ST1 = 0.5 * (AK1 - 1.0) / AK1;
		ST2 = 0.5 * (AK2 - 1.0) / AK2;
		if (PBIG >= (left.p - EPS))
			A1 = sqrt(left.rho * (0.5 * (AK1 + 1.0) * PBIG + 0.5 * (AK1 - 1.0) * left.p));
		else
		{
			PI = PBIG / left.p;
			A1 = 0.5 * (AK1 - 1.0) / AK1 * left.rho * C1 * (1.0 - PI) / (1.0 - pow(PI, ST1));
		}


		if (PBIG >= (right.p - EPS))
			A2 = sqrt(right.rho * (0.5 * (AK2 + 1.0) * PBIG + 0.5 * (AK2 - 1.0) * right.p));
		else
		{
			PI = PBIG / right.p;
			A2 = 0.5 * (AK2 - 1.0) / AK2 * right.rho * C2 * (1.0 - PI) / (1.0 - pow(PI, ST2));
		}

		UBIG = (A1 * left.u + A2 * right.u + left.p - right.p) / (A1 + A2);

		if ((left.p - EPS) < PBIG)
		{
			RBIG1 = left.rho * A1 / (A1 - left.rho * (left.u - UBIG));
			DL1 = left.u - A1 / left.rho;
			DL2 = DL1;
		}
		else
		{
			C1Z = C1 + 0.5 * (AK1 - 1.0) * (left.u - UBIG);
			RBIG1 = AK1 * PBIG / (C1Z * C1Z);
			DL1 = left.u - C1;
			DL2 = UBIG - C1Z;
		}


		if ((right.p - EPS) < PBIG)
		{
			RBIG2 = right.rho * A2 / (A2 + right.rho * (right.u - UBIG));
			DP1 = right.u + A2 / right.rho;
			DP2 = DP1;
		}
		else
		{
			C2Z = C2 - 0.5 * (AK2 - 1.0) * (right.u - UBIG);
			RBIG2 = AK2 * PBIG / (C2Z * C2Z);
			DP2 = right.u + C2;
			DP1 = UBIG + C2Z;
		}
		return;
	}
};