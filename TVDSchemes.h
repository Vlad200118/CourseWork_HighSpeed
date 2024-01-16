#pragma once
#include "Variables.h"
#include "GasModel.h"

class TVDScheme
{
public:
	enum TVDFunctions
	{
		Minmod = 1,
		Superbee = 2,
		vanLeer = 3,
		vanAlbado = 4
	};
	enum TVDSide
	{
		left = 0,
		right = 1
	};



private:
	LimiterArg r{0,0,0};
	TVDFunctions tvdFunc = vanAlbado;
	bool flagIsOn = false;

public:

	TVDScheme(const bool isOn, TVDFunctions tvdFunc = vanAlbado)
	{
		flagIsOn = isOn;
		this->tvdFunc = tvdFunc;
	}
	NaturalVar ReconstructVar(const NaturalVar& left, const NaturalVar& center , const NaturalVar& right, TVDSide side )
	{

		Limiter limter;

		if (side == TVDSide::left)
		{
			r.CalculateArg(left, center, right, false);
			limter = TVDLimiter();

			return (center + (right - left) * 0.25 * limter);
			//return (center + (center - left) * 0.5 * limter);
		}
		else
		{
			r.CalculateArg(left, center, right, true);
			//r.CalculateArg(left, center, right, false);
			limter = TVDLimiter();

			return (center - (right - left) * 0.25 * limter);
			//return (center - (right - center) * 0.5 * limter);
			//return (center - (center - left) * 0.5 * limter);

		}

	}


	Limiter TVDLimiter()
	{
		LimiterArg E{ 1.0,1.0,1.0 };
		switch (tvdFunc)
		{
		case TVDScheme::Minmod:
			return Limiter::GetMin(E * 2.0 / (r + 1.0), r / (r + 1) * 2.0);
			//return Limiter::GetMin(E, r);
			break;
		case TVDScheme::Superbee:
			return Limiter::GetMax(Limiter::GetMin(r/(r+1)*4.0, E * 2.0 /(r+1)), Limiter::GetMin(r / (r+1)*2.0, E *4.0 /(r+1)));
			//return Limiter::GetMax(Limiter::GetMin(r*2.0, E), Limiter::GetMin(r , E*2.0));
			break;
		case TVDScheme::vanLeer:
			return  r / (r + 1) / (r + 1) * 4.0;
			//return r * 2.0 / (r + E);
			break;
		default:
			return  r / (r * r + 1) * 2.0 ;
			//return (r * r + r) / (r * r + E);
 			break;

		}
	
	}

	bool isOn() const 
	{
		return flagIsOn;
	}


};