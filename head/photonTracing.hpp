#pragma once
#include "Eigen/Dense"
#include <vector>
// #include <omp.h>
#include <cstdlib>
#include <ctime>
#include "material.hpp"
#include "light.hpp"
#include "scene.hpp"
#include "kdTree.hpp"
//return number of photons
int globalPhotonTracing(Scene* scene, Map &photonMap, int n)
{
	Light* light = scene->lights[0];
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		bool firstHit = true;
		float lightPosPDF, lightDirPDF;
		Eigen::Vector3f lightPos, lightColor, lightDir;
		lightColor = light->SampleSurfacePos(lightPos, lightPosPDF);
		lightDir = light->SampleLightDir(lightDirPDF);
		Ray currRay(lightPos, lightDir);
		Interaction surfaceInteraction;
		while (1)
		{
			bool intersection = scene->intersection(&currRay, surfaceInteraction);
			if (intersection == false)
				break;
			if (((BSDF*)surfaceInteraction.material)->isSpecular == true)
			{
				firstHit = false;
				surfaceInteraction.inputDir = -currRay.m_Dir;
				//((BSDF*)surfaceInteraction.material)->sample(surfaceInteraction);
				Eigen::Vector3f L = surfaceInteraction.inputDir.normalized();
				Eigen::Vector3f N = surfaceInteraction.normal.normalized();
				float LdotN = L.dot(N);
				if (LdotN >= 0)
				{
					float n1 = 1.0f;
					float n2 = 1.5f;
					float cos1 = LdotN / (L.norm() * N.norm());
					float sin1 = sqrtf(1 - cos1 * cos1);
					float sin2 = n1 * sin1 / n2;
					float angle1 = acosf(cos1);
					float angle2 = asinf(sin2);
					float reflectRatio = 0.5f * ((sinf(angle1 - angle2) * sinf(angle1 - angle2)) / (sinf(angle1 + angle2) * sinf(angle1 + angle2)) + (tanf(angle1 - angle2) * tanf(angle1 - angle2)) / (tanf(angle1 + angle2) * tanf(angle1 + angle2)));
					float rand = (float)std::rand() / (float)RAND_MAX;
					if (rand < reflectRatio) //reflect
					{
						Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
					else //refract
					{
						Eigen::Vector3f outputDir = (n1 * LdotN / n2 - sqrtf(1 - n1 * n1 * (1 - LdotN * LdotN) / (n2 * n2))) * N - n1 * L / n2;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
				}
				else
				{
					float n1 = 1.5f;
					float n2 = 1.0f;
					float temp = 1 - n1 * n1 * (1 - LdotN * LdotN) / (n2 * n2);
					if (temp < 0) //total internal reflection
					{
						Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
					else
					{
						LdotN = -LdotN;
						N = -N;
						float cos1 = LdotN / (L.norm() * N.norm());
						float sin1 = sqrtf(1 - cos1 * cos1);
						float sin2 = n1 * sin1 / n2;
						float angle1 = acosf(cos1);
						float angle2 = asinf(sin2);
						float reflectRatio = 0.5f * ((sinf(angle1 - angle2) * sinf(angle1 - angle2)) / (sinf(angle1 + angle2) * sinf(angle1 + angle2)) + (tanf(angle1 - angle2) * tanf(angle1 - angle2)) / (tanf(angle1 + angle2) * tanf(angle1 + angle2)));
						float rand = (float)std::rand() / (float)RAND_MAX;
						if (rand < reflectRatio) //reflect
						{
							Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
							surfaceInteraction.outputDir = outputDir.normalized();
						}
						else //refract
						{
							Eigen::Vector3f outputDir = (n1 * LdotN / n2 - sqrtf(temp)) * N - n1 * L / n2;
							surfaceInteraction.outputDir = outputDir.normalized();
						}
					}
				}
				currRay.m_Ori = surfaceInteraction.entryPoint;
				currRay.m_Dir = surfaceInteraction.outputDir;
			}
			else
			{
				surfaceInteraction.inputDir = -currRay.m_Dir;
				((BSDF*)surfaceInteraction.material)->sample(surfaceInteraction);
				currRay.m_Ori = surfaceInteraction.entryPoint;
				currRay.m_Dir = surfaceInteraction.outputDir;
				if (firstHit)
					firstHit = false;
				else 
				{
					photonMap.store(surfaceInteraction.entryPoint, surfaceInteraction.inputDir);
					++count;
				}
				float rand = (float)std::rand() / (float)RAND_MAX;
				if (rand > 0.95f)
					break;
			}
		}
	}
	return count;
}
//return number of photons
int causticsPhotonTracing(Scene* scene, Map& photonMap, Eigen::Vector3f objPos, int n)
{
	Light* light = scene->lights[0];
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		float lightPosPDF;
		Eigen::Vector3f lightPos, lightColor, lightDir;
		while (1)
		{
			lightColor = light->SampleSurfacePos(lightPos, lightPosPDF);
			Eigen::Vector3f relativePos = objPos - lightPos;
			float rand1 = 2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f;
			lightDir.x() = relativePos.x() + rand1;
			float rand2 = 2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f;
			lightDir.y() = relativePos.y() + rand2;
			float rand3 = 2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f;
			lightDir.z() = relativePos.z() + rand3;
			Ray currRay(lightPos, lightDir);
			Interaction surfaceInteraction;
			bool intersection = scene->intersection(&currRay, surfaceInteraction);
			if (intersection == false)
				continue;
			if (((BSDF*)surfaceInteraction.material)->isSpecular == true)
				break;
		}
		Ray currRay(lightPos, lightDir);
		Interaction surfaceInteraction;
		while (1)
		{
			bool intersection = scene->intersection(&currRay, surfaceInteraction);
			if (intersection == false)
				break;
			if (((BSDF*)surfaceInteraction.material)->isSpecular == true)
			{
				surfaceInteraction.inputDir = -currRay.m_Dir;
				//((BSDF*)surfaceInteraction.material)->sample(surfaceInteraction);
				Eigen::Vector3f L = surfaceInteraction.inputDir.normalized();
				Eigen::Vector3f N = surfaceInteraction.normal.normalized();
				float LdotN = L.dot(N);
				if (LdotN >= 0)
				{
					float n1 = 1.0f;
					float n2 = 1.5f;
					float cos1 = LdotN / (L.norm() * N.norm());
					float sin1 = sqrtf(1 - cos1 * cos1);
					float sin2 = n1 * sin1 / n2;
					float angle1 = acosf(cos1);
					float angle2 = asinf(sin2);
					float reflectRatio = 0.5f * ((sinf(angle1 - angle2) * sinf(angle1 - angle2)) / (sinf(angle1 + angle2) * sinf(angle1 + angle2)) + (tanf(angle1 - angle2) * tanf(angle1 - angle2)) / (tanf(angle1 + angle2) * tanf(angle1 + angle2)));
					float rand = (float)std::rand() / (float)RAND_MAX;
					if (rand < reflectRatio) //reflect
					{
						Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
					else //refract
					{
						Eigen::Vector3f outputDir = (n1 * LdotN / n2 - sqrtf(1 - n1 * n1 * (1 - LdotN * LdotN) / (n2 * n2))) * N - n1 * L / n2;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
				}
				else
				{
					float n1 = 1.5f;
					float n2 = 1.0f;
					float temp = 1 - n1 * n1 * (1 - LdotN * LdotN) / (n2 * n2);
					if (temp < 0) //total internal reflection
					{
						Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
						surfaceInteraction.outputDir = outputDir.normalized();
					}
					else
					{
						LdotN = -LdotN;
						N = -N;
						float cos1 = LdotN / (L.norm() * N.norm());
						float sin1 = sqrtf(1 - cos1 * cos1);
						float sin2 = n1 * sin1 / n2;
						float angle1 = acosf(cos1);
						float angle2 = asinf(sin2);
						float reflectRatio = 0.5f * ((sinf(angle1 - angle2) * sinf(angle1 - angle2)) / (sinf(angle1 + angle2) * sinf(angle1 + angle2)) + (tanf(angle1 - angle2) * tanf(angle1 - angle2)) / (tanf(angle1 + angle2) * tanf(angle1 + angle2)));
						float rand = (float)std::rand() / (float)RAND_MAX;
						if (rand < reflectRatio) //reflect
						{
							Eigen::Vector3f outputDir = -L + 2.0f * LdotN * N;
							surfaceInteraction.outputDir = outputDir.normalized();
						}
						else //refract
						{
							Eigen::Vector3f outputDir = (n1 * LdotN / n2 - sqrtf(temp)) * N - n1 * L / n2;
							surfaceInteraction.outputDir = outputDir.normalized();
						}
					}
				}
				currRay.m_Ori = surfaceInteraction.entryPoint;
				currRay.m_Dir = surfaceInteraction.outputDir;
			}
			else
			{
				surfaceInteraction.inputDir = -currRay.m_Dir;
				((BSDF*)surfaceInteraction.material)->sample(surfaceInteraction);
				currRay.m_Ori = surfaceInteraction.entryPoint;
				currRay.m_Dir = surfaceInteraction.outputDir;
				photonMap.store(surfaceInteraction.entryPoint, surfaceInteraction.inputDir);
				++count;
				break;
			}
		}
	}
	return count;
}