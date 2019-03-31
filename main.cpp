#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <algorithm>
#include <iostream>
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"



#define M_PI 3.141592653589793238

std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0, 1);

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :x(x), y(y), z(z) {};
	double norm2() { return x * x + y * y + z * z; }
	void normalize() {
		double n = sqrt(norm2());
		x /= n;
		y /= n;
		z /= n;
	}
	Vector normalized() {
		Vector result(*this);
		result.normalize();
		return result;
	}

	double operator[](int i) {
		return *(&x + i);
	}

	double x, y, z;
};


// Les différents opérateurs nécessaires pour les vector:

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}


Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector operator-(const Vector& a) {
	return Vector(-a.x, -a.y, -a.z);
}

double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector operator*(const Vector& a, double b) {
	return Vector(a.x * b, a.y * b, a.z * b);
}

Vector operator*(double b, const Vector& a) {
	return Vector(a.x * b, a.y * b, a.z * b);
}

Vector operator/(const Vector& a, double b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}

Vector operator/(const Vector& a, const Vector& b) {
	return Vector(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}




class Ray {
public:
	Ray(const Vector& C, const Vector& u) {
		this->C = C;
		this->u = u;
	};
	Vector C, u;
};

//Les sphères et les maillages héritent de cette classe, cela permet de pouuvoir mettre les deux types d'objets dans la scène 
//puis d'appeler intersect de la même manières pour les deux objets
class Object {
public:
	Object() {}
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double& t, Vector &couleur_texture) = 0;

	bool miroir;
	Vector color;
	bool transparent;
	bool emissif;
};

class Sphere : public Object {
public:
	Sphere(const Vector& O, double R, const Vector& color, bool miroir, bool transparent = false, bool emissif = false) {
		this->O = O;
		this->R = R;
		this->color = color;
		this->miroir = miroir;
		this->transparent = transparent;
		this->emissif = emissif;
	};
	// L'intersection entre un rayon et une sphère
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double& t, Vector &couleur_texture) {
		couleur_texture = color;
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0) {
			double t1 = (-b - sqrt(delta)) / (2.*a);
			double t2 = (-b + sqrt(delta)) / (2.*a);
			if (t1 > 0) {
				P = r.C + r.u*t1;
				t = t1;
			}
			else {
				if (t2 > 0) {
					P = r.C + r.u*t2;
					t = t2;
				}
				else {
					return false;
				}
			}
		}
		N = P - O;
		N.normalize();
		return (delta >= 0);
	}

	Vector O;
	double R;

};

class Scene {
public:
	Scene() {}
	void addSphere(Object* s) {
		spheres.push_back(s);
	}
	// Cette fonction appelle tour à tour les méthodes intersect de chaque objet de la scène puis nous donne l'intersection la plus proche
	bool intersect(const Ray& r, Vector &P, Vector &N, int& object, double& t, Vector &couleur_texture) {
		double smallestT = 1E15;
		bool has_inter = false;
		for (unsigned int i = 0; i < spheres.size(); i++) {
			Vector Plocal, Nlocal, couleur_local;
			double tlocal;
			bool inter = spheres[i]->intersect(r, Plocal, Nlocal, tlocal, couleur_local);
			if (inter && tlocal < smallestT) {
				smallestT = tlocal;
				P = Plocal;
				N = Nlocal;
				couleur_texture = couleur_local;
				has_inter = true;
				object = i;
			}
		}
		t = smallestT;
		return has_inter;
	}

	// Cette fonction permet de trouver une direction aléatoire (mais qui a une plus grande probabilité de se trouver autour de N)
	Vector random_cos(const Vector& N) {
		double r1 = uniform(engine);
		double r2 = uniform(engine);
		double s = sqrt(1 - r2);
		Vector v(cos(2 * M_PI*r1)*s, sin(2 * M_PI*r1)*s, sqrt(r2));


		Vector T1;

		// On cherche T1 et T2 ortogonaux à N
		// On pose l'une des valeurs à 0 mais il faut que ce soit la valeur la plus petite (en valeur absolue)
		if (abs(N.x) <= abs(N.y) && abs(N.x) <= abs(N.z)) {
			T1 = Vector(0, -N.z, N.y);

		}
		else {
			if (abs(N.y) <= abs(N.x) && abs(N.y) <= abs(N.z)) {
				T1 = Vector(-N.z, 0, N.x);
			}
			else {
				T1 = Vector(-N.y, N.x, 0);
			}
		}
		T1.normalize();
		Vector T2 = cross(N, T1);

		return v.x*T1 + v.y*T2 + v.z*N;


	}

	Vector getColor(const Ray& ray, int nbrebonds, bool direct) { //direct a été défini pour éviter le pb de la lumière noire : si direct = true alors soit la lumière est directement en face de la caméra soit elle est vue à travers un mirroir si reflection diffuse on passe a false.
		
		Vector rayColor(0, 0, 0);
		if (nbrebonds <= 0) {
			return rayColor;
		}
		int object;
		Vector P, N, couleur_texture;
		double t;
		if (intersect(ray, P, N, object, t, couleur_texture)) { //On regarde s'il y a une intersection avec un objet de la scene : test en fait les intersections avec chaque objet de la scene
			//std::cout << object << std::endl;
			if (spheres[object]->emissif) { //Si c'est une lumière
				if (direct == true) { // et que l'on regarde à travers un miroir
					// On prend en compte la lumière (éviter point noir)
					rayColor = spheres[object]->color*I / (4.*M_PI*M_PI*dynamic_cast<Sphere*>(spheres[object])->R*dynamic_cast<Sphere*>(spheres[object])->R);
				}
				else { // Sinon la lumière est déjà prise en compte dans la contribution directe
					rayColor = Vector(0, 0, 0);
				}
			}
			else {
				if (spheres[object]->miroir) {// Pour obtenir un rayon réfléchi : r = i - 2<i,N>N avec i rayon incident et N la normale à la surface
					Vector R = ray.u - 2 * dot(ray.u, N)*N;
					Ray rayR(P + 0.001*N, R);
					return getColor(rayR, nbrebonds - 1, direct);
				}

				else {
					if (spheres[object]->transparent) {
						double n1 = 1;
						double n2 = 1.5;
						if (dot(N, ray.u) > 0) {
							N = -N;
							std::swap(n1, n2);

						}
						Vector Tt = (n1 / n2)*(ray.u - dot(ray.u, N)*N);
						double D = 1 - (n1 / n2)*(n1 / n2)*(1 - dot(ray.u, N)*dot(ray.u, N));

						Vector T;
						if (D < 0) {
							T = ray.u - 2.*dot(ray.u, N)*N;
						}
						else {
							Vector Tn = -sqrt(D)*N;
							T = Tn + Tt;
						}

						Ray rayR(P - 0.001*N, T);
						return getColor(rayR, nbrebonds - 1, direct);
					}
					else {

						//Eclairage direct 

						Vector dirLP = (P - L->O).normalized();
						Vector random_dir = random_cos(dirLP);
						Vector random_x = random_dir * L->R + L->O;
						Vector wi = (random_x - P).normalized();
						Vector Np = random_dir;
						double distance_lum2 = (random_x - P).norm2();


						Ray rayon_lumiere(P + 0.01*N, wi);
						Vector P_light, N_light;
						int sphere_id_light;
						double t_light;
						Vector color;
						bool ombre = intersect(rayon_lumiere, P_light, N_light, sphere_id_light, t_light, color);
						if (ombre && t_light*t_light < distance_lum2*0.99) { // Si un autre objet de la scene cache la lumière
							// C'est dans l'ombre
							rayColor = Vector(0, 0, 0);
						}
						else {
							// Sinon on envoie un rayon vers un point aléatoire de la lumière
							double a = std::max(0., dot(N, wi));
							double b = dot(Np, -wi);
							double c = dot(dirLP, random_dir);
							rayColor = (I / (4. *M_PI*distance_lum2) * std::max(0., dot(N, wi)) * dot(Np, -wi) / dot(dirLP, random_dir)) * couleur_texture;

						}

						// contribution indirecte
						// On envoie un rayon dans une direction aléatoire
						Vector reflechi = random_cos(N); //N est la normale autour de laquelle il y a plus de probabilité
						Ray ray_reflechi(P + 0.001*N, reflechi);
						rayColor = rayColor + spheres[object]->color*getColor(ray_reflechi, nbrebonds - 1, false);

					}

				}

			}
		}
		else {
			rayColor = Vector(0, 0, 0);
		}

		return rayColor;

	}

	std::vector<Object*> spheres;
	double I;
	Sphere *L;
};

class Triangle {
public:
	Triangle(const Vector &A, const Vector &B, const Vector &C) {
		this->A = A;
		this->B = B;
		this->C = C;
	}

	// Intersection entre un rayon et un triangle
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double& t, double &alpha, double &beta, double &gamma) {

		N = cross(B - A, C - A);
		if (dot(N, r.u) > 0) {
			N = -N;
		}
		t = dot(A - r.C, N) / dot(r.u, N);
		if (t > 0) {
			P = r.C + r.u*t;

			beta = (dot(P - A, B - A)*(C - A).norm2() - dot(C - A, B - A)*dot(P - A, C - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)*dot(B - A, C - A));
			gamma = ((B - A).norm2()*dot(P - A, C - A) - dot(P - A, B - A)*dot(B - A, C - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)*dot(B - A, C - A));
			alpha = (1 - gamma - beta);
			return(beta < 1 && beta>0 && gamma < 1 && gamma>0 && alpha < 1 && alpha>0);
		}
		else {
			return false;
		}
	}

	Vector A;
	Vector B;
	Vector C;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

// On ajoute cette classe pour diminuer le temps de calcul d'intersection avec un maillage
//Une boite englobante englobe un maillage
class BoiteEnglobante {
public:
	BoiteEnglobante(const Vector& vmin = Vector(-1E9, -1E9, -1E9), const Vector& vmax = Vector(1E9, 1E9, 1E9)) {
		this->vmax = vmax;
		this->vmin = vmin;
	}

	//Intersection entre un rayon et une boite englobante
	virtual bool intersect(const Ray& r) {
		double inter1_x = (vmin.x - r.C.x) / r.u.x;
		double inter2_x = (vmax.x - r.C.x) / r.u.x;
		double inter_min_x = std::min(inter1_x, inter2_x);
		double inter_max_x = std::max(inter1_x, inter2_x);

		double inter1_y = (vmin.y - r.C.y) / r.u.y;
		double inter2_y = (vmax.y - r.C.y) / r.u.y;
		double inter_min_y = std::min(inter1_y, inter2_y);
		double inter_max_y = std::max(inter1_y, inter2_y);

		double inter1_z = (vmin.z - r.C.z) / r.u.z;
		double inter2_z = (vmax.z - r.C.z) / r.u.z;
		double inter_min_z = std::min(inter1_z, inter2_z);
		double inter_max_z = std::max(inter1_z, inter2_z);

		if (std::min(std::min(inter_max_x, inter_max_y), inter_max_z) - std::max(std::max(inter_min_x, inter_min_y), inter_min_z) > 0) {
			return true;
		}
		else {
			return false;
		}
	}

	Vector vmin, vmax;

};

// Création d'un arbre de boites englobantes
// Un noeud est défini par une boite englobante, les indices des triangles qu'elle contient et par ses 2 noeuds fils
class BVHNode {
public:
	BVHNode(BoiteEnglobante b = BoiteEnglobante(), int debut = 0, int fin = 0) {
		this->boite = b;
		this->debut = debut;
		this->fin = fin;
		this->fd = NULL;
		this->fg = NULL;
	}
	BVHNode *fg, *fd;
	BoiteEnglobante boite;
	int debut, fin;
};

// Classe pour les maillages
class Geometry : public Object {
public:
	~Geometry() {


	}
	Geometry(const Vector &translate, double scale, const Vector& color, bool miroir = false, bool transparent = false, bool emissif = false) {
		this->color = color;
		this->miroir = miroir;
		this->transparent = transparent;
		this->emissif = false;
		readOBJ("pig.obj");
		for (int i = 0; i < vertices.size(); i++) {
			//vertices[i] = Vector(vertices[i].z, vertices[i].y, vertices[i].x)*scale + translate; //cochon de profil
			vertices[i] = Vector(vertices[i].x, vertices[i].y, vertices[i].z)*scale + translate; //cochon de face
		}
		//Pour cochon de profil :
		//for (int i = 0; i < indices.size(); i++) {
			//std::swap(indices[i].ni, indices[i].nk);
		//}

		this->boite = createBoite(0, indices.size());
		this->racine = BVHNode(boite, 0, indices.size());
		tri(indices, racine, 0, indices.size());
	};

	BoiteEnglobante createBoite(int debut, int fin) {
		//Boite englobante : 
		Vector vmin(1E9, 1E9, 1E9);
		Vector vmax(-1E9, -1E9, -1E9);
		for (unsigned int i = debut; i < fin; i++) {
			double vertx_min[] = { vmin.x, vertices[indices[i].vtxi].x, vertices[indices[i].vtxj].x, vertices[indices[i].vtxj].x };
			double verty_min[] = { vmin.y, vertices[indices[i].vtxi].y, vertices[indices[i].vtxj].y, vertices[indices[i].vtxj].y };
			double vertz_min[] = { vmin.z, vertices[indices[i].vtxi].z, vertices[indices[i].vtxj].z, vertices[indices[i].vtxj].z };
			vmin.x = *std::min_element(vertx_min, vertx_min + 4);
			vmin.y = *std::min_element(verty_min, verty_min + 4);
			vmin.z = *std::min_element(vertz_min, vertz_min + 4);

			double vertx_max[] = { vmax.x, vertices[indices[i].vtxi].x, vertices[indices[i].vtxj].x, vertices[indices[i].vtxj].x };
			double verty_max[] = { vmax.y, vertices[indices[i].vtxi].y, vertices[indices[i].vtxj].y, vertices[indices[i].vtxj].y };
			double vertz_max[] = { vmax.z, vertices[indices[i].vtxi].z, vertices[indices[i].vtxj].z, vertices[indices[i].vtxj].z };
			vmax.x = *std::max_element(vertx_max, vertx_max + 4);
			vmax.y = *std::max_element(verty_max, verty_max + 4);
			vmax.z = *std::max_element(vertz_max, vertz_max + 4);


		}
		return BoiteEnglobante(vmin, vmax);
	}

	// fonction donnée de lecture d'un fichier .obj
	void readOBJ(const char* obj, bool load_textures = false) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	// intersection entre un rayon et un maillage
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double& t, Vector &couleur_texture) {
		double smallestT = 1E15;
		bool has_inter = false;
		BVHNode bvh = racine;
		// On vérifie déjà si le rayon intersecte la boîte englobante du maillage total
		if (bvh.boite.intersect(r)) {
			std::vector<TriangleIndices> triangles;
			// On trouve la listes des triangles de la plus petite boîte englobante intersectée.
			intersectBoite(r, bvh, triangles);
			//On vérifie l'intersection avec chacun de ces triangles
			for (unsigned int i = 0; i < triangles.size(); i++) {
				Triangle triangle(vertices[triangles[i].vtxi], vertices[triangles[i].vtxj], vertices[triangles[i].vtxk]);
				//std::cout << triangles.size()<<std::endl;
				Vector Plocal, Nlocal;
				double tlocal, alpha, beta, gamma;
				bool inter = triangle.intersect(r, Plocal, Nlocal, tlocal, alpha, beta, gamma);
				if (inter && tlocal < smallestT) {
					smallestT = tlocal;
					P = Plocal;
					//N = Nlocal; // si on ne veut pas prendre en compte le lissage de phong
					//Lissage de phong :
					N = alpha * normals[triangles[i].ni] + beta * normals[triangles[i].nj] + gamma * normals[triangles[i].nk];
					N.normalize();
					has_inter = true;

					if (hasTexture==true) {// Si on veut afficher les textures
						int id_tex = 0; //Une seule texture ici donc on a pas vraiment besoin d'une liste, on prend directement le premier et seul élément
						int width = textures_width[id_tex];
						int height = textures_height[id_tex];
						
						// obtention des coordonnées uv pour la texture (on a un peu trafiquer pour avoir de bonnes couleurs)
						int x = (int) ((uvs[triangles[i].uvi] * alpha + uvs[triangles[i].uvj] * beta + uvs[triangles[i].uvk] * gamma).z*width)+60;
						int y = (int) ((uvs[triangles[i].uvi] * alpha + uvs[triangles[i].uvj] * beta + uvs[triangles[i].uvk] * gamma).x*height/2)+45;
						
						double rouge = (textures[id_tex][3 * (x * width + y)])/255.;
						double vert = (textures[id_tex][3 * (x * width + y) + 1])/255.;
						double bleu = (textures[id_tex][3 * (x * width + y) + 2])/255.;

						couleur_texture = Vector(rouge, vert, bleu);
					}
				}
			}
		}


		t = smallestT;
		return has_inter;
	}

	// intersection avec une boite englobante : récursive sur tout l'arbre bvh pour obtenir la plus petite liste de triangles à tester possible
	void intersectBoite(const Ray & r, const BVHNode & node, std::vector<TriangleIndices> & triangles) {
		if (node.fg) {
			//std::cout << "isParent" << std::endl;
			if (node.fd->boite.intersect(r)) {
				//std::cout << " ------- " << node.fd->debut << " " << node.fd->fin << std::endl << " ------ ";
				intersectBoite(r, *node.fd, triangles);
			}
			if (node.fg->boite.intersect(r)) {
				intersectBoite(r, *node.fg, triangles);
			}
		}
		else {
			//std::cout << " ------- " << node.debut << " " << node.fin << std::endl << " ------ ";
			for (int i = node.debut; i < node.fin; i++) {
				triangles.push_back(indices[i]);
			}
		}
	}

	// Création de l'arbre bvh
	void tri(std::vector<TriangleIndices> & tableau, BVHNode & node, int debut, int fin) {
		BoiteEnglobante boiteEng = node.boite;
		if (fin - debut >= 30) {
			Vector diag = boiteEng.vmax - boiteEng.vmin;
			double diagCoord[] = { diag.x, diag.y, diag.z };
			int dimCoupure = std::distance(diagCoord, std::max_element(diagCoord, diagCoord + 3));
			double valCoupure = boiteEng.vmin[dimCoupure] + diag[dimCoupure] / 2;
			int p = debut;
			for (int i = debut; i < fin; i++) {
				double center = (vertices[tableau[i].vtxi][dimCoupure] + vertices[tableau[i].vtxj][dimCoupure] + vertices[tableau[i].vtxk][dimCoupure]) / 3;
				if (center < valCoupure) {
					std::swap(tableau[i], tableau[p]);
					p++;
				}
			}


			BoiteEnglobante boite1 = createBoite(debut, p);
			//std::cout << debut << " " << p << " " << fin << std::endl; // Pour afficher les limites de chaque noeud
			node.fg = new BVHNode(boite1, debut, p);

			BoiteEnglobante boite2 = createBoite(p, fin);
			node.fd = new BVHNode(boite2, p, fin);
			if (fin != p) {
				tri(tableau, *node.fg, debut, p);
			}
			if (debut != p) {
				tri(tableau, *node.fd, p, fin);
			}

			

		}

	}

	//ajout d'une image de texture
	void addTexture(const char* filename) {
		textures.resize(textures.size() + 1);
		textures_width.resize(textures_width.size() + 1);
		textures_height.resize(textures_height.size() + 1);

		FILE* f;
		fopen_s(&f, filename, "rb");
		unsigned char info[54];
		fread(info, sizeof(unsigned char), 54, f); //Lecture 54-byte header

		textures_width[textures_width.size() - 1] = *(int*)&info[18]; // extraction largeur image
		textures_height[textures_height.size() - 1] = *(int*)&info[22]; //Hauteur

		int size = 3 * textures_width[textures_width.size() - 1] * textures_height[textures_height.size() - 1];
		textures[textures.size() - 1].resize(size); // 3 bits par pixel, 3 couleurs
		fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
		fclose(f);

		for (int i = 0; i < size; i += 3) {
			std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
		}
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	BoiteEnglobante boite;
	BVHNode racine;
	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> textures_width;
	std::vector<int> textures_height;
	bool hasTexture;


};


int main() {
	int W = 512; //Largeur de l'image finale
	int H = 512; //hauteur de l'image finale
	int Nrays = 128; // Le nombre de rayon envoyés par pixel (plus prend plus de temps mais réduit le bruit)

	double fov = 60 * M_PI / 180.; // angle de field of view de la caméra


	Vector C(0, 0, 55); // position de la caméra
	Vector L(-10, 20, 40); // position de la lumière (son centre)

	// L'intensité de la lumière :
	double I = 1000000000;

	//Rayon de la lumière
	double RLum = 5;

	// Distance focus
	double dFocus = 30; // Les objet à une distance dFocus de la caméra seront nets
	double depthOfField_intensity = 1.5; //intensité de de l'effet dof, correspond à l'ouverture du diaphragme

	Vector pixelColor(0, 0, 0);

	// On crée tous nos objets
	Sphere lumiere(Vector(-10, 20, 40), RLum, Vector(1., 1., 1.), false, false, true);
	Sphere s(Vector(-15, 0, 0), 10, Vector(1., 0., 0.), true, false);
	Sphere ssol(Vector(0, -1000, 0), 1000 - 10, Vector(0.3, 0.015, 0.2), false, false);
	Sphere splafond(Vector(0, 1000, 0), 1000 - 60, Vector(1., 1., 1.), false, false);
	Sphere smur1(Vector(1000, 0, 0), 1000 - 60, Vector(1., 1., 1.), false, false);
	Sphere smur2(Vector(-1000, 0, 0), 1000 - 60, Vector(1., 1., 1.), false, false);
	Sphere smur3(Vector(0, 0, 1000), 1000 - 60, Vector(1., 1., 1.), false, false);
	Sphere smur4(Vector(0, 0, -1000), 1000 - 60, Vector(1., 1., 1.), false, false);
	Sphere s2(Vector(8, -5, 30), 5, Vector(1., 1., 1.), false, false);
	Sphere s3(Vector(3, -8, 34), 2, Vector(1., 1., 1.), false, true);
	//Triangle t0(Vector(-10, 0, 20), Vector(0, 0, 20), Vector(0, 10, 20), Vector(1., 1., 1.), false, false, false);
	//Geometry geom(Vector(0,-10,18), 0.6, Vector(1., 1., 1.));
	Geometry geom(Vector(0, -10, 10), 10, Vector(1., 1., 1.));
	geom.hasTexture = true;
	geom.addTexture("TexturePig.bmp");
	Scene scene;
	scene.L = &lumiere;
	scene.I = I;
	// On les ajoute à notre scene
	scene.addSphere(&lumiere);
	scene.addSphere(&s);
	scene.addSphere(&s2);
	scene.addSphere(&s3);
	scene.addSphere(&ssol);
	scene.addSphere(&splafond);
	scene.addSphere(&smur1);
	scene.addSphere(&smur2);
	scene.addSphere(&smur3);
	scene.addSphere(&smur4);
	//scene.addSphere(&t0);
	scene.addSphere(&geom);

	std::vector<unsigned char> image(W*H * 3, 0);
// on fait les calculs en parallèle
#pragma omp parallel for schedule(dynamic,1)
	// On parcours chaque pixel
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector pixelColor;
			for (int rayid = 0; rayid < Nrays; rayid++) { // pour chaque pixels on envoie Nrays rayons
				// pas toujours au centre du pixel (gaussienne autour du centre) (Anti-aliasing)
				double r1 = uniform(engine);
				double r2 = uniform(engine);
				double offset_x = cos(2 * M_PI*r1)*sqrt(-2 * log(r2));
				double offset_y = sin(2 * M_PI*r1)*sqrt(-2 * log(r2));
				Vector u(j - W / 2 + offset_x + 0.5, H / 2 - i + offset_y + 0.5, -W / (2 * tan(fov / 2)));
				u.normalize();
				// Depth of field : envoie du rayon aléatoirement dans l'ouverture du diaphragme
				double Cx_offset = (uniform(engine) - 0.5)*depthOfField_intensity;
				double Cy_offset = (uniform(engine) - 0.5)*depthOfField_intensity;

				Vector focal = C + dFocus * u;
				Vector origine = C + Vector(Cx_offset, Cy_offset, 0);
				Ray ray(origine, (focal - origine).normalized()); // Notre rayon part de C (le point de vue de l'image) dans la direction u qui définit un pixel avec une légère randomisation gaussienne pour éviter les éscaliers (plus flouté)
				pixelColor = pixelColor + scene.getColor(ray, 5, true); // Le 5 veut dire que l'on autorise 5 rebonds
			}
			pixelColor = pixelColor * 1. / Nrays; // On fait la moyenne sur les Nrays rayons envoyés sur le pixel

			image[(i*W + j) * 3 + 0] = std::min(255., std::pow(pixelColor.x, 0.45)); //La puissance 0.45 est la correction gamma, pour mieux visualiser les couleurs sur un écran
			image[(i*W + j) * 3 + 1] = std::min(255., std::pow(pixelColor.y, 0.45));
			image[(i*W + j) * 3 + 2] = std::min(255., std::pow(pixelColor.z, 0.45));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);



	return 0;
}
