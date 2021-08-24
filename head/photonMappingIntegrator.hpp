#pragma once
#include "photonTracing.hpp"
#include "integrator.hpp"
#include "material.hpp"
#include "kdTree.hpp"
#include <cmath>
// #include <omp.h>
#define PHOTON_NUM 1000000
class PhotonMappingIntegrator : public Integrator
{
public:
	PhotonMappingIntegrator(Scene* scene, Camera* camera)
		: Integrator(scene, camera)
	{
	}

	// main render loop
	void render(Map &global,Map &caustic) override
	{
		
		for (int dx = 0; dx < camera->m_Film.m_Res.x(); dx++)
		{
			for (int dy = 0; dy < camera->m_Film.m_Res.y(); dy++)
			{
				Eigen::Vector3f L(0.0f, 0.0f, 0.0f);
				Ray ray = camera->generateRay(dx, dy);
				Interaction surfaceInteraction;
				int samples = 200;
				for (int i = 0; i < samples; ++i)
				{
					L += 3.0f * radiance(&surfaceInteraction, &ray);	//direct light
				}
				L = L / samples;

				Ray specular_Ray = camera->generateRay(dx, dy);		//specular light
				Interaction specular_SurfaceInteraction;
				bool specular_interaction = scene->intersection(&specular_Ray, specular_SurfaceInteraction);
				if (specular_interaction) {
					if (((BSDF*)specular_SurfaceInteraction.material)->isSpecular == true) {
						float materialPDF, lightPDF;
						Eigen::Vector3f out_Dir, surface_Color, pos,surface_Norm, materialBRDF;
						Eigen::Vector3f color(0, 0, 0);
						Eigen::Vector3f beta(1, 1, 1);
						surface_Norm = specular_SurfaceInteraction.normal;
						specular_SurfaceInteraction.inputDir = -specular_Ray.m_Dir.normalized();
						out_Dir = -specular_SurfaceInteraction.inputDir + 2.0f * specular_SurfaceInteraction.inputDir.dot(surface_Norm) * surface_Norm;
						specular_Ray.m_Ori = specular_SurfaceInteraction.entryPoint;
						specular_Ray.m_Dir = out_Dir;
						for (int i = 0; i < 5; i++) {
							if (scene->lights[0]->isHit(&specular_Ray, &specular_SurfaceInteraction))
							{
								color += beta.cwiseProduct(scene->lights[0]->m_Color);
							}
							specular_interaction = scene->intersection(&specular_Ray, specular_SurfaceInteraction);
							if(specular_interaction){
								color += beta.cwiseProduct(scene->lights[0]->m_Color).cwiseProduct(specular_SurfaceInteraction.surfaceColor);
								specular_SurfaceInteraction.inputDir = -specular_Ray.m_Dir.normalized();
								materialPDF = ((BSDF*)specular_SurfaceInteraction.material)->sample(specular_SurfaceInteraction);
								materialBRDF = ((BSDF*)specular_SurfaceInteraction.material)->eval(specular_SurfaceInteraction);
								if (materialPDF == 0.0f || (materialBRDF.x() == 0.0f && materialBRDF.y() == 0.0f && materialBRDF.z() == 0.0f))
									break;
								specular_Ray.m_Ori = specular_SurfaceInteraction.entryPoint;
								specular_Ray.m_Dir = specular_SurfaceInteraction.outputDir;
								beta = beta.cwiseProduct(materialBRDF) * std::fabsf(specular_SurfaceInteraction.outputDir.dot(specular_SurfaceInteraction.normal)) / materialPDF;
							}
							else
								break; 
						}
						L += 0.01f * color;
					}
				}

				Ray ray_photon = camera->generateRay(dx, dy);
				Interaction surfaceInteraction_photon;
				bool intersection_photon = scene->intersection(&ray_photon,surfaceInteraction_photon);
				Eigen::Vector3f irradiancePhoton(0, 0, 0);
				Eigen::Vector3f surfaceNormPhoton = surfaceInteraction_photon.normal.normalized();
				Eigen::Vector3f surfaceColorPhoton = surfaceInteraction_photon.surfaceColor;
				surfaceInteraction_photon.inputDir = -ray_photon.m_Dir;
				((BSDF*)surfaceInteraction_photon.material)->sample(surfaceInteraction_photon);
				Eigen::Vector3f BRDF = ((BSDF*)surfaceInteraction_photon.material)->eval(surfaceInteraction_photon);
				int photon_num = global.stored_photons;
				if (intersection_photon) {
					if (((BSDF*)surfaceInteraction_photon.material)->isSpecular != true) {
						Nearest_photons np(1000,surfaceInteraction_photon.entryPoint,3);
						global.locate_photons(&np, 1);
						Photon** photons = np.get_photons();
						for (int i = 1;i <= np.curr_num;i++) {
							if (photons[i]->dir.dot(surfaceNormPhoton) > 0) {
								Eigen::Vector3f input_dir = photons[i]->dir.normalized();
								irradiancePhoton += surfaceColorPhoton * (input_dir.dot(surfaceNormPhoton.normalized()));
							}
							irradiancePhoton *= (1 / M_PIf) * (1 / np.dist[0])/photon_num;
							L += 20.0f * irradiancePhoton.cwiseProduct(BRDF);
						}
					}
				}

				Ray ray_caustic = camera->generateRay(dx, dy);
				Interaction surfaceInteraction_caustic;
				bool intersection_caustic = scene->intersection(&ray_caustic, surfaceInteraction_caustic);
				Eigen::Vector3f irradianceCaustics(0, 0, 0);
				Eigen::Vector3f surfaceNormCaustics = surfaceInteraction_caustic.normal.normalized();
				Eigen::Vector3f	surfaceColorCaustics = surfaceInteraction_caustic.surfaceColor;
				surfaceInteraction_caustic.inputDir = -ray_caustic.m_Dir;
				((BSDF*)surfaceInteraction_caustic.material)->sample(surfaceInteraction_caustic);
				Eigen::Vector3f BRDF_caustic = ((BSDF*)surfaceInteraction_caustic.material)->eval(surfaceInteraction_caustic);
				int photon_num1 = caustic.stored_photons;
				if (intersection_caustic) {
					if (((BSDF*)surfaceInteraction_caustic.material)->isSpecular != true) {
						Nearest_photons np1(1000, surfaceInteraction_caustic.entryPoint, 2);
						caustic.locate_photons(&np1, 1);
						Photon** photons_Caustic = np1.get_photons();
						for (int i = 1;i <= np1.curr_num;i++) {
							if (photons_Caustic[i]->dir.dot(surfaceNormCaustics) > 0) {
								Eigen::Vector3f input_dir = photons_Caustic[i]->dir.normalized();
								irradianceCaustics += surfaceColorCaustics * (input_dir.dot(surfaceNormCaustics.normalized()));
							}
						}
						irradianceCaustics *= (1 / M_PIf) * (1 / np1.dist[0])/ photon_num1;
						L += 20.0f * irradianceCaustics.cwiseProduct(BRDF_caustic);
					}
				}

				camera->setPixel(dx, dy, L);
			}
		}
		
	}

	// radiance of a specific point
	Eigen::Vector3f radiance(Interaction* interaction, Ray* ray) override
	{
		//// TODO: Calculate color here
		Eigen::Vector3f L(0.0f, 0.0f, 0.0f);
		//Eigen::Vector3f beta(1.0f, 1.0f, 1.0f);
		bool specularBounce = false;
		Ray currRay = *ray;
		Interaction surfaceInteraction;
		bool intersection = scene->intersection(&currRay, surfaceInteraction);
		if (scene->lights[0]->isHit(&currRay, &surfaceInteraction) == true)
			L += scene->lights[0]->m_Color;
		if (intersection) {
			Eigen::Vector3f lightPos, lightColor, materialBRDF;
			float lightPDF, materialPDF;
			lightColor = scene->lights[0]->SampleSurfacePos(lightPos, lightPDF);
			Eigen::Vector3f lightDir = lightPos - surfaceInteraction.entryPoint;
			Ray shadowRay(surfaceInteraction.entryPoint, lightDir, 1e-3f, lightDir.norm());
			if (!scene->intersection(&shadowRay))
				L += (lightColor.cwiseProduct(surfaceInteraction.surfaceColor)) / lightPDF;
			return 0.1f * L;
		}
	}
};