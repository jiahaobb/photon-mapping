#pragma once
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include "ray.hpp"
#include "interaction.hpp"
#define M_PIf 3.14159265358979323846f

class Light
{
public:
	Light(Eigen::Vector3f pos, Eigen::Vector3f color)
		: m_Pos(std::move(pos)),
		m_Color(std::move(color))
	{
	}
	virtual ~Light() = default;

	// Sample a point on the light's surface, if light is not delta light
	// Set PDF and the surface position
	// Return color of light
	virtual Eigen::Vector3f SampleSurfacePos(Eigen::Vector3f& sampled_lightPos, float& pdf)=0;
	virtual Eigen::Vector3f SampleLightDir(float& pdf) = 0;
	// Determine if light is hit, if light is not delta light
	virtual bool isHit(Ray* ray, Interaction* interaction) = 0;
	
	Eigen::Vector3f m_Pos;
	Eigen::Vector3f m_Color;
};

class AreaLight:public Light
{
public:
	AreaLight(Eigen::Vector3f pos, Eigen::Vector3f color) : Light(pos, color)
	{
	}
	
	Eigen::Vector3f SampleSurfacePos(Eigen::Vector3f& sampled_lightPos, float& pdf) override {
		// TODO
		float rand1 = (float)std::rand() / (float)RAND_MAX;
		float r = std::sqrtf(rand1);
		float rand2 = (float)std::rand() / (float)RAND_MAX;
		float phi = 2.0f * M_PIf * (rand2);
		float x = r * std::cosf(phi);
		float z = r * std::sinf(phi);
		Eigen::Vector3f localPos(x, 0.0f, z);
		sampled_lightPos = localPos + m_Pos;
		pdf = 1.0f / M_PIf;
		return m_Color;
	}

	Eigen::Vector3f SampleLightDir(float& pdf) override
	{
		float rand1 = (float)std::rand() / (float)RAND_MAX;
		float r = std::sqrtf(rand1);
		float rand2 = (float)std::rand() / (float)RAND_MAX;
		float phi = 2.0f * M_PIf * (rand2);
		float x = r * std::cosf(phi);
		float z = r * std::sinf(phi);
		float y = -sqrtf(1.0f - x * x - z * z);
		pdf = 1.0f / M_PIf;
		return { x,y,z };
	}

	bool isHit(Ray* ray, Interaction* interaction) override {
		// TODO
		if (ray->m_Dir.y() == 0)
			return false;
		float t = (m_Pos.y() - ray->m_Ori.y()) / ray->m_Dir.y();
		if (t > ray->m_fMax || t < ray->m_fMin)
			return false;
		if (interaction->isInteraction == true && interaction->entryDist < t)
			return false;
		float x = ray->m_Ori.x() + t * ray->m_Dir.x();
		float z = ray->m_Ori.z() + t * ray->m_Dir.z();
		if ((x - m_Pos.x()) * (x - m_Pos.x()) + (z - m_Pos.z()) * (z - m_Pos.z()) > 1)
			return false;
		return true;
	}
};