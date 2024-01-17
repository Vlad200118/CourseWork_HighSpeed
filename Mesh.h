#pragma once
#include <vector>
#include "BoundaryConditions.h"

class Mesh
{
private:
	std::vector <double> cellCenter;
	std::vector <double> surface;
	std::vector <double> cellVolume;
	std::vector <double> oldCellVolume;

	std::vector <double> oldPosFace;
	std::vector <double> newPosFace;

	std::vector <double> velocityFace;

	double leftBoundCoord{ 0 }, rightBoundCoord{0};
	size_t _NX{ 0 };
public:
	void BuildMesh(const size_t NX, const double length)
	{
		_NX = NX;
		cellCenter.resize(NX+2);

		newPosFace.resize(NX+1);
		oldPosFace.resize(NX + 1);

		surface.resize(NX+1);
		velocityFace.resize(NX + 1);
		cellVolume.resize(NX+2);
		oldCellVolume.resize(NX + 2);

		leftBoundCoord = 0.0;
		rightBoundCoord = length;
		CalculateCoordinates();
		CalculateSquare();
		CalculateVolume();
		std::copy(cellVolume.begin(), cellVolume.end(), oldCellVolume.begin());
		CalculateVelocityFace();

	}

	~Mesh()
	{
		cellCenter.clear();
		surface.clear();
		cellVolume.clear();

		oldPosFace.clear();
		newPosFace.clear();

		velocityFace.clear();
	}

	double S(int index) const
	{
		return surface[index];
	}
	double CellVolume(int index) const
	{
		return cellVolume[index];
	}
	double OldCellVolume(int index) const
	{
		return oldCellVolume[index];
	}
	double Node(int index) const
	{
		return newPosFace[index];
	}
	double CellCenter(int index) const
	{
		return cellCenter[index];
	}

	size_t GetNX() const 
	{
		return _NX;
	}
	double FaceVelocity(int index)const
	{
		return velocityFace[index];
	}


	double GetMinStep() const
	{
		double minStep = 1000000;
		for (size_t i = 2; i < cellCenter.size() - 1; i++)
		{
			minStep = std::min(minStep, abs(cellCenter[i] - cellCenter[i - 1]));
		}
		return minStep;
	}


	void RebuildMesh(BoundaryConditions& boundCond, const double& stepTime)
	{
		leftBoundCoord += boundCond.GetBoundary("left")->BoundaryVelocity()* stepTime;
		rightBoundCoord += boundCond.GetBoundary("right")->BoundaryVelocity()* stepTime;

		CalculateCoordinates();
		CalculateSquare();
		CalculateVolume();
		CalculateVelocityFace(stepTime);
	
	}

	void WriteMesh(std::ostream& out)
	{
		out << "x\tS\n";
		for (size_t i = 0; i < surface.size(); i++)
		{
			out << newPosFace[i] << "\t" << surface[i] << std::endl;
		}
	}

private:

	void CalculateSquare()
	{
		//Нужен закон

		double D = 1.9;

		double x1{ 0.4 * rightBoundCoord }, x2{ 0.6 * rightBoundCoord };

		for (size_t i = 0; i < surface.size(); i++)
		{
			surface[i] =  0.0 * newPosFace[i] + 1.0;

			/*if (newPosFace[i] < x1)
			{
				surface[i] = (1 - D) / (x1 - leftBoundCoord) * newPosFace[i] + D;
			}
			else if (newPosFace[i] >= x1 && newPosFace[i] <= x2)
			{
				surface[i] = 1.0;
			}
			else
			{
				surface[i] = (D - 1) / (leftBoundCoord + rightBoundCoord - x2) / (leftBoundCoord + rightBoundCoord - x2) * (newPosFace[i] - x2) * (newPosFace[i] - x2) + 1.0;
			}*/
			
		}



	}



	void CalculateVelocityFace(const double& stepTime = 0.0)
	{
		if (stepTime != 0.0)
		{
			for (size_t i = 0; i < velocityFace.size(); i++)
			{
				velocityFace[i] = (newPosFace[i] - oldPosFace[i]) / stepTime;
			}
		}
	}


	void CalculateCoordinates()
	{

		std::copy(newPosFace.begin(), newPosFace.end(), oldPosFace.begin());

		double dX = (rightBoundCoord - leftBoundCoord) / _NX;

		for (size_t i = 0; i < cellCenter.size(); i++)
		{
			if (i == 0)
				cellCenter[i] = leftBoundCoord - dX / 2.0;
			else
				cellCenter[i] = cellCenter[i - 1] + dX;
		}

		for (size_t i = 0; i < newPosFace.size(); i++)
		{
			if (i == 0)
				newPosFace[i] = leftBoundCoord;
			else
				newPosFace[i] = newPosFace[i - 1] + dX;
		}


	}

	void CalculateVolume()
	{
		std::copy(cellVolume.begin(), cellVolume.end(), oldCellVolume.begin());
		for (size_t i = 1; i < cellVolume.size()-1; i++)
		{
			cellVolume[i] = abs(newPosFace[i] - newPosFace[i-1]) * (surface[i] + surface[i - 1]) * 0.5;
		}
		cellVolume[0] = cellVolume[1];
		cellVolume[cellVolume.size()-1] = cellVolume[cellVolume.size() - 2];
	}

};
