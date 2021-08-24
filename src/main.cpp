
#include <iostream>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "parallelogram.hpp"
#include "scene.hpp"
#include "camera.hpp"
#include "photonMappingIntegrator.hpp"
#include "triangleMesh.hpp"
#include "photonTracing.hpp"
inline float clamp(float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline unsigned char toInt(float x) { unsigned char c(pow(clamp(x), 1 / 2.2) * 255 + .5);return c; }

int main()
{
	/*
	 * 1. Camera Setting
	 */
	Eigen::Vector3f cameraPosition(0, 0, 10);
	Eigen::Vector3f cameraLookAt(0, 0, 0);
	Eigen::Vector3f cameraUp(0, 1, 0);
	float verticalFov = 45;
	Eigen::Vector2i filmRes(500, 500);

	Camera camera(cameraPosition, cameraLookAt, cameraUp, verticalFov, filmRes);
	

	/*
	 * 2. Basic geometry setting
	 * Setting the cornell box.
	 */
	Eigen::Vector3f p_p0(-10, -10, -10);
	Eigen::Vector3f p_s0(1, 0, 0);
	Eigen::Vector3f p_s1(0, 1, 0);
	Eigen::Vector3f p_normal(0, 0, 1);
	Eigen::Vector3f p_color(0.7, 0.7, 0.7);
	Parallelogram backWall(p_p0, p_s0 * 20, p_s1 * 20, p_normal, p_color);

	p_p0 = Eigen::Vector3f(-6, -7, -11);
	p_s0 = Eigen::Vector3f(0, 0, 1);
	p_s1 = Eigen::Vector3f(0, 1, 0);
	p_normal = Eigen::Vector3f(1, 0, 0);
	p_color = Eigen::Vector3f(1, 0, 0);
	Parallelogram leftWall(p_p0, p_s0 * 20, p_s1 * 20, p_normal, p_color);

	p_p0 = Eigen::Vector3f(6, -7, -11);
	p_s0 = Eigen::Vector3f(0, 0, 1);
	p_s1 = Eigen::Vector3f(0, 1, 0);
	p_normal = Eigen::Vector3f(-1, 0, 0);
	p_color = Eigen::Vector3f(0, 1, 0);
	Parallelogram rightWall(p_p0, p_s0 * 20, p_s1 * 20, p_normal, p_color);

	p_p0 = Eigen::Vector3f(-7, -6, -11);
	p_s0 = Eigen::Vector3f(1, 0, 0);
	p_s1 = Eigen::Vector3f(0, 0, 1);
	p_normal = Eigen::Vector3f(0, 1, 0);
	p_color = Eigen::Vector3f(0.7, 0.7, 0.7);
	Parallelogram floor(p_p0, p_s0 * 20, p_s1 * 20, p_normal, p_color);

	p_p0 = Eigen::Vector3f(-7, 6, -11);
	p_s0 = Eigen::Vector3f(1, 0, 0);
	p_s1 = Eigen::Vector3f(0, 0, 1);
	p_normal = Eigen::Vector3f(0, -1, 0);
	p_color = Eigen::Vector3f(0.7, 0.7, 0.7);
	Parallelogram ceiling(p_p0, p_s0 * 20, p_s1 * 20, p_normal, p_color);

	/*
	 * 3. Triangle mesh setting
	 */
	std::string obj_1_filePos("../resources/p.obj");
	Eigen::Vector3f obj_1_color(1,1,1);
	Eigen::Affine3f obj_1_transform;
	obj_1_transform = Eigen::Translation3f(3, -2, -8) * Eigen::Scaling(0.5f);
	TriangleMesh mesh_1(obj_1_color, obj_1_filePos);
	mesh_1.applyTransformation(obj_1_transform);
	mesh_1.buildUniformGrid();


	/*
	 * 4. Light setting
	 */
	//Light light(Eigen::Vector3f(0, 5.99, 5), Eigen::Vector3f(1, 1, 1));
	AreaLight light(Eigen::Vector3f(0.0f, 5.8f, -5.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f));
	

	/*
	 * 5. Material setting
	 */
	BSDF* diffuseMat = new IdealDiffuse();
	BSDF* specularMat = new IdealSpecular();
	backWall.material = diffuseMat;
	floor.material = diffuseMat;
	leftWall.material = diffuseMat;
	rightWall.material = diffuseMat;
	ceiling.material = diffuseMat;
	mesh_1.material = specularMat;


	/*
	 * 5. Scene integration
	 */
	Scene scene;
	scene.addShape(&backWall);
	scene.addShape(&leftWall);
	scene.addShape(&rightWall);
	scene.addShape(&floor);
	scene.addShape(&ceiling);
	scene.addShape(&mesh_1);
	scene.addLight(&light);
	Map globalPhoton(100000, { 1.0f,1.0f,1.0f }), causticsPhoton(100000, { 1.0f,1.0f,1.0f });
	globalPhotonTracing(&scene, globalPhoton,10000);
	causticsPhotonTracing(&scene, causticsPhoton, { 3,-2,8 }, 10000);
	globalPhoton.balance();
	causticsPhoton.balance();
	/*std::cout << "size of pos " << photon << std::endl;
	for (int i = 0; i < pos.size(); ++i)
		std::cout << pos[i].x() << " " << pos[i].y() << " " << pos[i].z() << std::endl;*/
	/*
	 * 6. Select and execute integrator
	 */
	PhotonMappingIntegrator integrator(&scene, &camera);
	integrator.render(globalPhoton, causticsPhoton);

	/*
	 * 7. Output image to file
	 */
	std::string outputPath = "./output.png";
	std::vector<unsigned char> outputData;
	outputData.reserve(int(filmRes.x() * filmRes.y() * 3));
	for (Eigen::Vector3f v : camera.m_Film.pixelSamples)
	{
		for (int i = 0; i < 3; i++)
		{
			outputData.push_back(toInt(v[i]));
		}
	}
	stbi_flip_vertically_on_write(true);
	stbi_write_png(outputPath.c_str(), filmRes.x(), filmRes.y(), 3, outputData.data(), 0);

	return 0;
}