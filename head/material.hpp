#pragma once
#include "Eigen/Dense"
#include "interaction.hpp"
#define M_PIf 3.14159265358979323846f
class BSDF
{
public:
	BSDF()
	{
		isSpecular = false;
	}

	// Evaluate the BSDF
	// The information in @Interaction contains ray's direction, normal
	// and other information that you might need
	virtual Eigen::Vector3f eval(Interaction& _interact) = 0;

	// Sample a direction based on the BSDF
	// The sampled direction is stored in @Interaction
	// The PDF of this direction is returned
	virtual float sample(Interaction& _interact) = 0;

	// Mark if the BSDF is specular
	bool isSpecular;
	float clamp(float x) { return x < 0.0f ? 0.0f : x > 1.0f ? 1.0f : x; }
	Eigen::Vector3f rotate(Eigen::Vector3f oriVec, Eigen::Vector3f oriNormal, Eigen::Vector3f newNormal)
	{
		Eigen::Matrix3f I, tmp, rotation;
		Eigen::Vector3f a = oriNormal.normalized();
		Eigen::Vector3f b = newNormal.normalized();
		if (a == -b)
			return -oriVec;
		if (a == b)
			return oriVec;
		Eigen::Vector3f v = a.cross(b);
		float s = v.norm();
		float c = a.dot(b);
		tmp << 0.0f, -v.z(), v.y(), v.z(), 0.0f, -v.x(), -v.y(), v.x(), 0.0f;
		I << 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f;
		rotation = I + tmp + tmp * tmp / (1.0f + c);
		return rotation * oriVec;
	}
};


class IdealDiffuse : public BSDF
{
public:

	Eigen::Vector3f eval(Interaction& _interact)
	{
		// TODO
		Eigen::Vector3f result;
		Eigen::Vector3f N = _interact.normal.normalized();
		Eigen::Vector3f V = _interact.outputDir.normalized();
		Eigen::Vector3f L = _interact.inputDir.normalized();
		float NdotL = N.dot(L);
		float NdotV = N.dot(V);
		if (NdotL <= 0.0f || NdotV <= 0.0f)
		{
			result.x() = 0.0f;
			result.y() = 0.0f;
			result.z() = 0.0f;
			return result;
		}
		result = clamp(NdotL) * _interact.surfaceColor / M_PIf;
		return result;
	};

	float sample(Interaction& _interact)
	{
		// TODO
		float rand1 = (float)std::rand() / (float)RAND_MAX;
		float r = std::sqrtf(rand1);
		float rand2 = (float)std::rand() / (float)RAND_MAX;
		float phi = 2.0f * M_PIf * (rand2);
		float x = r * std::cosf(phi);
		float y = r * std::sinf(phi);
		float z = std::sqrtf(std::fmaxf(0.0f, 1.0f - x * x - y * y));
		Eigen::Vector3f outputDir(x, y, z);
		outputDir = rotate(outputDir, { 0.0f,0.0f,1.0f }, _interact.normal);
		_interact.outputDir = outputDir.normalized();
		return std::fabsf(_interact.inputDir.dot(_interact.normal)) / M_PIf;
	};
};

class IdealSpecular : public BSDF
{
public:

	IdealSpecular() {
		isSpecular = true;
	}

	Eigen::Vector3f eval(Interaction& _interact)
	{
		// TODO
		float epslon = 1e-4;
		Eigen::Vector3f result;
		Eigen::Vector3f N = _interact.normal.normalized();
		Eigen::Vector3f V = _interact.outputDir.normalized();
		Eigen::Vector3f L = _interact.inputDir.normalized();
		float NdotL = N.dot(L);
		Eigen::Vector3f tmp = -L + 2.0f * NdotL * N;
		if (std::fabsf(tmp.x() - V.x()) >= epslon && std::fabsf(tmp.y() - V.y()) >= epslon && std::fabsf(tmp.y() - V.y()) >= epslon)
		{
			result.x() = 0.0f;
			result.y() = 0.0f;
			result.z() = 0.0f;
			return result;
		}
		result = clamp(NdotL) * _interact.surfaceColor;
		return result;
	};

	float sample(Interaction& _interact)
	{
		// TODO
		Eigen::Vector3f L = _interact.inputDir.normalized();
		Eigen::Vector3f N = _interact.normal.normalized();
		float LdotN = L.dot(N);
		Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
		_interact.outputDir = outputDir.normalized();
		return 1.0f;
	};
};