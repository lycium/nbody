// n-body simulator
// thomas ludwig 2015
// thomas.ludwig@gmail.com

#define USE_KAHAN 1

#include <GL/glut.h>
#include <unistd.h>
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN

#include <random>

#include <algorithm>
#include <thread>
#include <mutex>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <stdint.h>

#include <vector>

#if USE_KAHAN
#include "kahan.h"
#endif

volatile bool quit = false;
volatile bool reset = false;

const int image_scale = 2;
const int image_div = 3;
const int xres = 1920 * image_scale / image_div;
const int yres = 1080 * image_scale / image_div;

double start_time;


template<typename real>
class vec2
{
public:
	real x, y;

	inline vec2<real>() { }
	inline vec2<real>(const vec2<real> & v) : x(v.x), y(v.y) { }
	inline vec2<real>(real v) : x(v), y(v) { }
	inline vec2<real>(real x_, real y_) : x(x_), y(y_) { }

	inline vec2<real> operator+(const vec2<real> & rhs) const { return vec2<real>(x + rhs.x, y + rhs.y); }
	inline vec2<real> operator-(const vec2<real> & rhs) const { return vec2<real>(x - rhs.x, y - rhs.y); }
	inline vec2<real> operator*(const real rhs) const { return vec2<real>(x * rhs, y * rhs); }

	inline const vec2<real> & operator=(const vec2<real> & rhs) { x = rhs.x; y = rhs.y; return *this; }

	inline const vec2<real> & operator+=(const vec2<real> & rhs) { x += rhs.x; y += rhs.y; return *this; }
	inline const vec2<real> & operator-=(const vec2<real> & rhs) { x -= rhs.x; y -= rhs.y; return *this; }
	inline const vec2<real> & operator*=(const real rhs)  { x *= rhs; y *= rhs; return *this; }

	inline real length2() const { return x * x + y * y; }
	inline real length() const { return std::sqrt(length2()); }

	inline void normalise(const real len) { const real s = len / length(); x *= s; y *= s; }
};

template<typename real>
inline static double dot(const vec2<real> & lhs, const vec2<real> & rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }

typedef vec2<double> vec2d;
typedef vec2<float> vec2f;


template<typename real>
class vec3
{
public:
	real x, y, z;

	inline vec3<real>() {x=0;y=0;z=0; }
	inline vec3<real>(const vec3<real> & v) : x(v.x), y(v.y), z(v.z) { }
	inline vec3<real>(real v) : x(v), y(v), z(v) { }
	inline vec3<real>(real x_, real y_, real z_) : x(x_), y(y_), z(z_) { }

	inline vec3<real> operator+(const vec3<real> & rhs) const { return vec3<real>(x + rhs.x, y + rhs.y, z + rhs.z); }
	inline vec3<real> operator-(const vec3<real> & rhs) const { return vec3<real>(x - rhs.x, y - rhs.y, z - rhs.z); }
	inline vec3<real> operator*(const real rhs)  const { return vec3<real>(x * rhs, y * rhs, z * rhs); }

	inline const vec3<real> & operator=(const vec3<real> & rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

	inline const vec3<real> & operator+=(const vec3<real> & rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	inline const vec3<real> & operator-=(const vec3<real> & rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	inline const vec3<real> & operator*=(const real rhs)  { x *= rhs; y *= rhs; z *= rhs; return *this; }

	inline real length2() const { return x * x + y * y + z * z; }
	inline real length() const { return std::sqrt(length2()); }

	inline void normalise(const real len) { const real s = len / length(); x *= s; y *= s; z *= s; }
};

template<typename real>
inline static double dot(const vec3<real> & lhs, const vec3<real> & rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }


#if USE_KAHAN
typedef vec3<KahanAdder<double> > vec3d; // Use Kahan summation to minimise roundoff error
#else
typedef vec3<double> vec3d;
#endif

typedef vec3<float> vec3f; // Used for path storage and simple stuff that doesn't need tons of precision

struct timespec ts;

class nbody
{
public:
	nbody(uint32_t seed = 17468) { init(seed); }

	void init(uint32_t seed)
	{
		std::lock_guard<std::mutex> lock(mutex);

		t_current = 0;
		t_step = 0.0000025;
		num_steps = 0;

		std::mt19937 mt(seed);

		vec3f init_points[num_particles];
		int restarts = 0;
		for (int z = 0; z < num_particles; ++z)
		{
			vec3f p;
			size_t tries = 0;
			const size_t max_tries = 1 << 19;
			do
			{
				p = vec3f(mt(), mt(), mt()) * (2.0f / 4294967296) - 1;
				if (p.length2() > 1) continue;

				bool far_enough = true;
				const float t_min = 0.96f;
				const float t_min2 = t_min * t_min;
				for (int i = 0; i < z; ++i) if ((p - init_points[i]).length2() < t_min2) { far_enough = false; break; }
				if (!far_enough) continue;

				break;
			} while (++tries < max_tries);

			if (tries >= max_tries)
			{
				++restarts;
				z = 0;
				continue;
			}

			init_points[z] = p;
		}

		std::cout << "init with seed " << seed << " took " << restarts << " restarts" << std::endl;

		for (int z = 0; z < num_particles; ++z)
		{
			pos[z] = vec3d(init_points[z].x, init_points[z].y, init_points[z].z);

			const double v0 = mt() * (2.0 / 4294967296) - 1;
			const double v1 = mt() * (2.0 / 4294967296) - 1;
			const double v2 = mt() * (2.0 / 4294967296) - 1;
			const vec3d d(v0, v1, v2);
			vel[z] = d * (2 / (d.length() + 1.0)); //d * (1.0 / d.length());

			mass[z] = 1.0;

			particle_paths[z].resize(0);
		}

		for (int z = 0; z < num_particles; ++z)
			pos[z] *= 0.2;
	}

	void predictionStep(
		const vec3d * const __restrict v_in,
		const vec3d * const __restrict a_in,
		const double dt,
		vec3d * const __restrict v_out,
		vec3d * const __restrict a_out)
	{
		for (int z = 0; z < num_particles; ++z)
		{
			p_tmp[z] = pos[z] + v_in[z] * dt;
			v_out[z] = vel[z] + a_in[z] * dt;
		}
		
		for (int i = 0; i < num_particles; ++i)
		{
			const vec3d & p_i = p_tmp[i];
			vec3d F_i ;
			for (int j = 0; j < num_particles; ++j)
			{
				if (j == i) continue;
				const vec3d & s = p_tmp[j] - p_i;
				const double r2 = s.length2(), r = std::sqrt(r2);
				//const double f = mass[j] / (r2 * r);
				const double f = 1 / (r2 * r); // HACK mass = 1

				const vec3d a = s * f;
				F_i += a;
			}
			//a_out[i] = F_i * (1 / mass[i]); // TODO precompute mass inverses?
			a_out[i] = F_i; // HACK mass = 1
		}
	}

	void step()
	{
		double t_sum = 0;
		int sub_steps = 0;

		do
		{
			if ((sub_steps % 1024) == 0)
			{
				if (reset) { init(clock_gettime(CLOCK_MONOTONIC,&ts) - start_time); reset = false; }
				else if (quit) { return; }
			}

			// Perform a first pass at the current time,
			// keeping track of the maximum force for time step control
			double max_force = 0;
			for (int i = 0; i < num_particles; ++i)
			{
				const vec3d p_i = pos[i];
				vec3d F_i;
				for (int j = 0; j < num_particles; ++j)
				{
					if (j == i) continue;
					const vec3d s = pos[j] - p_i;
					const double r2 = s.length2(), r = std::sqrt(r2);
					//const double f = mass[j] / (r2 * r);
					const double f = 1.0 / (r2 * r); // HACK mass = 1
					const vec3d a = s * f;

					F_i += a;
					//max_force = std::max(max_force, f / mass[i]);
					max_force = std::max(max_force, f);
				}
				//k1_acc[i] = F_i * (1 / mass[i]); // TODO precompute mass inverses?
				k1_acc[i] = F_i; // HACK mass = 1
			}

			// Compute the time step depending on the maximum force
			const double min_step = 1e-11;
			const double step_scale = 32.0 * 64 * 4;
			double step_scaled = std::min(t_step, std::max(min_step, t_step / max_force * step_scale));
			if (t_sum + step_scaled > t_step)
				step_scaled = t_step - t_sum;
			const double step_size = step_scaled;

			// Do 3 prediction steps for Runge Kutta 4th order integration
			// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
			predictionStep(   vel, k1_acc, step_size * 0.5, k2_vel, k2_acc);
			predictionStep(k2_vel, k2_acc, step_size * 0.5, k3_vel, k3_acc);
			predictionStep(k3_vel, k3_acc, step_size * 1.0, k4_vel, k4_acc);

			// Do the final time step using the weighted sum of the predicted velocity and acceleration vectors
			// Separate loops to avoid write data load/store hazards in the velocity vector (might not help perf? might hurt?)
			const double k = step_size / 6.0;
			for (int z = 0; z < num_particles; ++z) { pos[z] += (   vel[z] + (k2_vel[z] + k3_vel[z]) * 2 + k4_vel[z]) * k; }
			for (int z = 0; z < num_particles; ++z) { vel[z] += (k1_acc[z] + (k2_acc[z] + k3_acc[z]) * 2 + k4_acc[z]) * k; }

			t_sum += step_size;
			++sub_steps;
		}
		while (t_sum < t_step);

		//if (num_steps > 5120)
		{
			std::lock_guard<std::mutex> lock(mutex);

			for (int z = 0; z < num_particles; ++z)
				particle_paths[z].push_back(vec3f(
					(float)pos[z].x,
					(float)pos[z].y,
					(float)pos[z].z));
		}

		if (std::log10(sub_steps) > 5.0)
		//if (num_steps % 128 == 0)
		{
			std::cout << "step " << num_steps << " took 10 ^ " << std::log10(sub_steps) << " substeps!" << std::endl;

			//if (false)
			{
				std::pair<double, int> force_particle_pairs[num_particles];
				for (int i = 0; i < num_particles; ++i)
				{
					const vec3d p_i = pos[i];
					vec3d F_i;
					for (int j = 0; j < num_particles; ++j)
					{
						if (j == i) continue;
						const vec3d s = pos[j] - p_i;
						const double r2 = s.length2(), r = std::sqrt(r2);
						const double f = 1 / (r2 * r); // HACK mass = 1
						F_i += s * f;
					}
					force_particle_pairs[i] = std::make_pair(F_i.length(), i);
				}
				std::sort(force_particle_pairs, force_particle_pairs + num_particles, std::greater<std::pair<double, int> >());
				const double log10_max_force = std::log10(force_particle_pairs[0].first);
				for (int i = 0; i < num_particles; ++i)
					if (log10_max_force - std::log10(force_particle_pairs[i].first) < 2)
						std::cout << "particle " << (force_particle_pairs[i].second + 1) << " experiences force 10 ^ " << std::log10(force_particle_pairs[i].first) << std::endl;
			}
		}

		++num_steps;
	}

	const static int num_particles = 9;

	std::vector<vec3f> particle_paths[num_particles];

	double t_current;
	double t_step;

	size_t num_steps;

	vec3d pos[num_particles];
	vec3d vel[num_particles];

	double mass[num_particles];

	vec3d p_tmp[num_particles];

	vec3d                        k1_acc[num_particles];
	vec3d k2_vel[num_particles], k2_acc[num_particles];
	vec3d k3_vel[num_particles], k3_acc[num_particles];
	vec3d k4_vel[num_particles], k4_acc[num_particles];

	std::mutex mutex;
};

nbody grav;

std::thread compute_thread;
void ComputeThread() { while (!quit) { grav.step(); } }

std::vector<vec3f> path_point_pos; // Buffer for passing the moving path points to OpenGL in one call


inline static uint32_t wang_hash(uint32_t x) { x = (x ^ 12345391) * 2654435769; x ^= (x << 6) ^ (x >> 26); x *= 2654435769; x += (x << 5) ^ (x >> 12); return x; }

vec3f cam_lookat;
int lookat_path_idx;

std::vector<vec3f> obj_verts;
std::vector<uint32_t> obj_tris;

void renderScene()
{
	usleep(16); // Max 62.50fps to avoid spinning on the lock too hard

	const double t0 = clock_gettime(CLOCK_MONOTONIC,&ts) - start_time;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	{
		std::lock_guard<std::mutex> lock(grav.mutex);

		const double display_interval = 128 * 64;
		const double t_wrap = std::fmod(t0 * 0.25, display_interval);
		const int num_path_points = grav.particle_paths[0].size();
		const int num_display_pts = (int)(num_path_points / display_interval);
		path_point_pos.resize(num_display_pts);

		const vec3f new_lookat = (grav.particle_paths[lookat_path_idx].size() > 0) ? grav.particle_paths[lookat_path_idx][grav.particle_paths[lookat_path_idx].size() - 1] : vec3f(0);

		const float blend_speed = 0.02f;
		cam_lookat = cam_lookat * (1 - blend_speed) + new_lookat * blend_speed;

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		const double a = t0 * 0.00005;
		vec3f campos = cam_lookat + vec3f(0, 1, 0) * 0.1f + vec3f((float)std::cos(a), 0, (float)std::sin(a)) * 0.1f;
		vec3f upvec = vec3f(0, 1, 0);
		gluLookAt(campos.x, campos.y, campos.z, cam_lookat.x, cam_lookat.y, cam_lookat.z, upvec.x, upvec.y, upvec.z);

		glPushMatrix();

		for (int i = 0; i < nbody::num_particles; ++i)
		{
			const float path_r = wang_hash(i * 3 + 0) / 4294967296.0f;
			const float path_g = wang_hash(i * 3 + 1) / 4294967296.0f;
			const float path_b = wang_hash(i * 3 + 2) / 4294967296.0f;
			glColor3f(path_r, path_g, path_b);

			glLineWidth(2.0f);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, &grav.particle_paths[i][0]);
			glDrawArrays(GL_LINE_STRIP, 0, num_path_points);
			glDisableClientState(GL_VERTEX_ARRAY);

			for (int z = 0; z < num_display_pts; ++z)
			{
				const double i_max = num_path_points - 1;
				const double i_wrap = t_wrap + z * display_interval;
				const int i0 = (int)i_wrap;
				if (i0 >= num_path_points - 1) break;
				const vec3f p0 = grav.particle_paths[i][i0];

				path_point_pos[z] = p0;
			}

			glPointSize(6.0f);
			const float k = 0.35f;
			glColor3f(1 - (1 - path_r) * k, 1 - (1 - path_g) * k, 1 - (1 - path_b) * k);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, &path_point_pos[0]);
			glDrawArrays(GL_POINTS, 0, num_display_pts);
			glDisableClientState(GL_VERTEX_ARRAY);

			glPointSize(8.0f);
			glColor4f(1, 1, 1, 1);
			glBegin(GL_POINTS);
			{
				const vec3f p0 = grav.particle_paths[i][grav.particle_paths[i].size() - 1];
				glVertex3f(p0.x, p0.y, p0.z);
			}
			glEnd();
		}
	}
	glPopMatrix();

	glutSwapBuffers();
}

void changeSize(int w, int h)
{
	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (double)w / h, 0.01, 100.0);
}

void keyboardFunc(unsigned char key, int x, int y)
{
	if (key == 27)
	{
		quit = true;
		compute_thread.join();

		exit(0);
	}
	else if (key >= 49 && key <= 58)
	{
		lookat_path_idx = key - 49;
		return;
	}
	else if (key == 32)
	{
		reset = true;
	}
	//std::cout << "key = " << key * 1.0 << std::endl;
}

int main(int argc, char ** argv)
{
	// init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	//glutInitWindowPosition((3840 - xres) / 2, (2160 - yres) / 2);
	glutInitWindowSize(xres, yres);
	glutCreateWindow("press 1-9 to track particle, space to randomise");

	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);
	glutKeyboardFunc(keyboardFunc);

	glEnable(GL_DEPTH_TEST);	// Enable z-buffering
	glDisable(GL_CULL_FACE);	// Disable backface culling
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH); glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POINT_SMOOTH); glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	cam_lookat = 0;
	lookat_path_idx = 0;


	start_time = clock_gettime(CLOCK_MONOTONIC,&ts);

	compute_thread = std::thread(ComputeThread);

	// Enter GLUT event processing cycle. Anything after this doesn't get called / gets leaked...
	glutMainLoop();

	return 0;
}
