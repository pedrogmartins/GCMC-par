#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

struct particle_t {
    double x;           // Position X
    double y;           // Position Y
    double z;           // Position Z
    int index;          // Index of Atom
    std::string atom;   // Atomic identity
    std::vector<int> neighbors;  // List of indices of nearest neighbors

};

class PairPANNAMP {
public:
  PairPANNAMP();  //default constructor
  virtual ~PairPANNAMP(); //original implementation 
  
  int get_parameters(const std::string& directory, const std::string& filename);
  int get_input_line(std::ifstream*, std::string*, std::string*);
  void coeff(const std::string& directory, const std::string& filename);
  void compute_gvect(std::vector<particle_t>& particles, int index, int numneigh, double *G, double* dGdx, int n_atoms);
  float compute(std::vector<particle_t>& particles, int n_atoms); // had virtual on it initially
  double Gradial_d(double, int, double*);
  double Gangular_d(double, double, double, int, int, double*);
  double compute_network(double*, double*, int);

  //virtual 
  void settings(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);

  

  struct parameters{
    int Nspecies;
    // Gvector parameters
    float *eta_rad;
    float Rc_rad;
    float Rs0_rad;
    float Rsst_rad;
    int RsN_rad;
    float *eta_ang;
    float Rc_ang;
    float Rs0_ang;
    float Rsst_ang;
    int RsN_ang;
    int *zeta;
    int ThetasN;
    std::string* species;
    int gsize;
    float * Rs_rad;
    float * Rs_ang;
    float * Thetas;

    // Network parameters
    int *Nlayers;
    int **layers_size;
    int **layers_activation;

    // Useful precalculated quantities
    float cutmax;
    float *twoeta_rad;
    float *zeta_half;
    float iRc_rad;
    float iRc_rad_half;
    float iRc_ang;
    float iRc_ang_half;
    float *Rsi_rad;
    float *Rsi_ang;
    float *Thi_cos;
    float *Thi_sin;
    int **typsh;
  };


protected:
  // Gvect and NN parameters
  struct parameters par;
  // The network [species, layers, array]
  double ***network;


};


void read(std::string filename, std::vector<particle_t>& particles, double& a, double& b, double& c, int& n_atoms);
void fill_ghost_particles(std::vector<particle_t>& particles, double a, double b, double c, int n_atoms);

#endif // COMMON_H
