#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include<array> 


//TO FIX - LATTICE VECTORS PROCESSING NOT SUITABLE FOR NON
//   RECTANGULAR UNIT CELLS 
int n_atoms; //Need to compute the original number of particles

void read(std::string filename, std::vector<particle_t>& particles, double& a, double& b, double& c, int& n_atoms) {
    std::ifstream input_file(filename);
    std::string line;

    // Skip the first two lines of the VASP file
    std::getline(input_file, line);
    std::getline(input_file, line);

    // Read the lattice vectors and calculate their lengths
    double ax, ay, az;
    double bx, by, bz;
    double cx, cy, cz;
    input_file >> ax >> ay >> az;
    input_file >> bx >> by >> bz;
    input_file >> cx >> cy >> cz;
    a = sqrt(ax*ax + ay*ay + az*az);
    b = sqrt(bx*bx + by*by + bz*bz);
    c = sqrt(cx*cx + cy*cy + cz*cz);

    // Skip the next line of the VASP file
    std::getline(input_file, line);
    std::getline(input_file, line);
    std::getline(input_file, line);

    // Read the number of atoms
    std::istringstream iss(line);
    int num_mg, num_h, num_c, num_o;
    iss >> num_mg >> num_h >> num_c >> num_o;

    // Calculate the total number of atoms
    int num_atoms = num_mg + num_h + num_c + num_o;

    // Skip the next two lines of the VASP file
    std::getline(input_file, line);

    // Read the atomic positions and identities
    for (int i = 0; i < num_atoms; i++) {
        double x_frac, y_frac, z_frac;
        std::string atom_line;
        int index = i; 
        input_file >> x_frac >> y_frac >> z_frac;
        std::getline(input_file, atom_line);
        std::istringstream iss(atom_line);
        std::string atom;
        iss >> atom;
        double x_cart = ax * x_frac + bx * y_frac + cx * z_frac;
        double y_cart = ay * x_frac + by * y_frac + cy * z_frac;
        double z_cart = az * x_frac + bz * y_frac + cz * z_frac;
        particle_t particle = {x_cart, y_cart, z_cart, index, atom};
        //std::cout << "Atom: " << atom << std::endl;
        particles.push_back(particle);
    }

    n_atoms = num_atoms;

    fill_ghost_particles(particles, a, b, c, num_atoms);

    int debug = 0;


}


double distance(const particle_t& p1, const particle_t& p2, double a, double b, double c) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;

    //if (dx > 0.5*a) dx -= a;
    //if (dx < -0.5*a) dx += a;
    //if (dy > 0.5*b) dy -= b;
    //if (dy < -0.5*b) dy += b;
    //if (dz > 0.5*c) dz -= c;
    //if (dz < -0.5*c) dz += c;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double distance_ghost(const particle_t& p1, const particle_t& p2, double a, double b, double c) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;

    //if (dx > 0.5*a) dx -= a;
    //if (dx < -0.5*a) dx += a;
    //if (dy > 0.5*b) dy -= b;
    //if (dy < -0.5*b) dy += b;
    //if (dz > 0.5*c) dz -= c;
    //if (dz < -0.5*c) dz += c;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void fill_ghost_particles(std::vector<particle_t>& particles, double a, double b, double c, int n_atoms) {

    std::vector <particle_t> part_18_ghosts;

    //Loop over particles
    for (int i = 0; i < n_atoms; i++) {
        particle_t& pi = particles[i];

        //Populate the total number of particles vector with
        // ghost particles (make sure to update the for loop later
        // as now we don't have to loop over all particles when 
        // computing energies and forces, ghost particles can be ignored) 

        //Now we need to add the neighbors from periodic images
        for (int ix = -1; ix <= 1; ix++) {
            for (int iy = -1; iy <= 1; iy++) {
                for (int iz = -1; iz <= 1; iz++) {
                    if (ix == 0 && iy == 0 && iz == 0) continue;
                    particle_t ghost = pi; // Make a copy of the particle
                    ghost.x = ghost.x + ix*a - iy*-10.9346914291381836;
                    ghost.y += iy*18.9394416809082031;
                    ghost.z += iz*c;
                    //Assign unique index to each ghost particle, we first start 
                    //   at the original index, and add 84 spacings according to 
                    //   to how we traverse the ghost boxes uniquely
                    //int ghost_index = (n_atoms - 1) + ((iz + 1) + 3*(iy + 1) + 9*(iz+1))*84 + (ghost.index + 1);
                    
                    // Was getting exc bad access errors, so let`s create
                    // index based on the numebr of particles already in the 
                    // array
                    ghost.index = particles.size() + 1;
                    particles.push_back(ghost);
                    //if (i == 17 && ix == -1 && iy == 1 && iz == 0) { //}  && ghost.y == -1.59304 && ghost.z == 6.17916) {
                        
                    //    int debug_integer = 0;
                    //}
                    //if (i == 18) {
                    //    part_18_ghosts.push_back(ghost);
                    //}
                    
                    //for (int j = 0; j < particles.size(); j++) {
                    //    particle_t& pj = particles[j];
                    //    double dist = distance(ghost, pj, a, b, c);
                    //    if (dist < rc) {
                    //        pi.neighbors.push_back(j);
                    //    }
                    //}

                    int debug = 0;
                }
            }
        }
    }
    int debug_master = 0;
}


void find_neighbors(std::vector<particle_t>& particles, double rc, double a, double b, double c, int n_atoms) {
    
    //Now we have all the ghost particles, let's loop over each particle
    // and find all of its neighbours

    float dist = 0;

    //First loop over particles in the original cell only
    for (int i = 0; i < n_atoms; i++) {
        particle_t& pi = particles[i];
        //Now loop over all possible neighbors
        for (int j = 0; j < particles.size(); j++) {
            if (i == j) continue;
            particle_t& pj = particles[j];
            if (j < n_atoms) {
                dist = distance(pi, pj, a, b, c);
                if (i == 0 && j == 6) {
                    int asdas = dist;

                }
            }
            else {
                dist = distance_ghost(pi, pj, a, b, c);
            }
            //Now, we need two function to consider the distance,
            //  one if the particle is within the initial unit cell, 
            // another if the particle is a ghost one
            if (dist < rc) {
                pi.neighbors.push_back(j);
            }
        }

    }
}

//Constructor implementation 

PairPANNAMP::PairPANNAMP() {
    // constructor implementation
}

int PairPANNAMP::get_input_line(std::ifstream* file, std::string* key, std::string* value){
  std::string line;
  int parsed = 0; int vc = 1;
  while(!parsed){
    std::getline(*file,line);
    // Exit on EOF
    if(file->eof()) return 0;
    // Exit on bad read
    if(file->bad()) return -1;
    // Remove spaces
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    // Skip empty line
    if(line.length()==0) continue;
    // Skip comments
    if(line.at(0)=='#') continue;
    // Parse headers
    if(line.at(0)=='['){
      *value = line.substr(1,line.length()-2);
      return 1;
    }
    // Check if we have version information:
    if(line.at(0)=='!') { vc=0 ;}
    // Look for equal sign
    std::string eq = "=";
    size_t eqpos = line.find(eq);
    // Parse key-value pair
    if(eqpos != std::string::npos){
      *key = line.substr(0,eqpos);
      *value = line.substr(eqpos+1,line.length()-1);
      if (vc == 0) { vc = 1 ; return 3; }
      return 2;
    }
    std::cout << line << std::endl;
    parsed = 1;
  }
  return -1;
}

int PairPANNAMP::get_parameters(const std::string& directory, const std::string& filename)
{
    //const double panna_pi = 3.14159265358979323846;
    // Parsing the potential parameters
    std::ifstream params_file;
    std::ifstream weights_file;
    std::string key, value;
    std::string dir_string(directory);
    std::string param_string(filename);
    std::string file_string(dir_string+"/"+param_string);
    std::string wfile_string;

    // Initializing some parameters before reading:
    par.Nspecies = -1;
    // Flags to keep track of set parameters
    int Npars = 17;
    int parset[Npars];
    for(int i=0;i<Npars;i++) parset[i]=0;
    int *spset;
    std::string version = "v0"; 
    int gversion=0;
    double tmp_eta_rad; 
    double tmp_eta_ang ; 
    int tmp_zeta ;

    params_file.open(file_string.c_str());
    int section = -1;
    int parseint = get_input_line(&params_file,&key,&value);
    
    while(parseint>0){
        // Parse line
        if(parseint==1){
        // Gvect param section
        if(value=="GVECT_PARAMETERS"){ section = 0; }
        // For now other sections are just species networks
        else {
            // First time after params are read: do checks
            if(section==0){
            if (gversion==0){
                if(parset[5]==0){
                // Set steps if they were omitted
                par.Rsst_rad = (par.Rc_rad - par.Rs0_rad) / par.RsN_rad; parset[5]=1;}
                if(parset[10]==0){
                par.Rsst_ang = (par.Rc_ang - par.Rs0_ang) / par.RsN_ang; parset[10]=1;}}
            else if (gversion==1){parset[5]=1 ; parset[10]=1;}
            // Check that all parameters have been set
            for(int p=0;p<Npars;p++){
                if(parset[p]==0){
                std::cout << "Parameter " << p << " not set!" << std::endl;  return -1; } }
            // Calculate Gsize
            par.gsize = par.Nspecies * par.RsN_rad + (par.Nspecies*(par.Nspecies+1))/2 * par.RsN_ang * par.ThetasN;
            } //section 0 ended
            int match = 0;
            for(int s=0;s<par.Nspecies;s++){
            // If species matches the list, change section
            if(value==par.species[s]){
                section = s+1;
                match = 1;
            }
            }
            if(match==0){
            std::cout << "Species " << value << " not found in species list." << std::endl;
            return -2;
            }
        }
        }// A header is parsed
        else if(parseint==2){
        // Parse param section
        if(section==0){
        std::string comma = ",";
            if(key=="Nspecies"){
            par.Nspecies = std::atoi(value.c_str());
            // Small check
            if(par.Nspecies<1){
                std::cout << "Nspecies needs to be >0." << std::endl;
                return -2; }
            parset[0] = 1;
            // Allocate species list
            par.species = new std::string[par.Nspecies];
            // Allocate network quantities
            par.Nlayers = new int[par.Nspecies];
            par.layers_size = new int*[par.Nspecies];
            par.layers_activation = new int*[par.Nspecies];
            network = new double**[par.Nspecies];
            // Keep track of set species
            spset = new int[par.Nspecies];
            for(int s=0;s<par.Nspecies;s++) {
                par.Nlayers[s] = -1;
                spset[s]=0; } }
            else if(key=="species"){
            //std::string comma = ",";
            size_t pos = 0;
            int s = 0;
            // Parse species list
            while ((pos = value.find(comma)) != std::string::npos) {
                if(s>par.Nspecies-2){
                std::cout << "Species list longer than Nspecies." << std::endl;
                return -2; }
                par.species[s] = value.substr(0, pos);
                value.erase(0, pos+1);  s++; }
            if(value.length()>0){
                par.species[s] = value; s++; };
            if(s<par.Nspecies){
                std::cout << "Species list shorter than Nspecies." << std::endl;
                return -2; }
            parset[1] = 1; }
            // Common features are read.
            // From here on what will be read depends on the gversion
            if(gversion == 0){ // Potentials compatible with OPENKIM
            std::cout << "G Version is " << gversion << std::endl;
            if(key=="eta_rad"){
                tmp_eta_rad = std::atof(value.c_str()); parset[2] = 0; }
            else if(key=="Rc_rad"){
                par.Rc_rad = std::atof(value.c_str()); parset[3] = 1; }
            else if(key=="Rs0_rad"){
                par.Rs0_rad = std::atof(value.c_str());  parset[4] = 1; }
            else if(key=="Rsst_rad"){
                par.Rsst_rad = std::atof(value.c_str()); parset[5] = 1; }
            else if(key=="RsN_rad"){
                par.RsN_rad = std::atoi(value.c_str()); parset[6] = 1; 
                par.eta_rad = new float[par.RsN_rad]; 
                par.twoeta_rad = new float[par.RsN_rad];
                par.Rs_rad  = new float[par.RsN_rad];
                for(int i=0;i<par.RsN_rad;i++) par.eta_rad[i]=tmp_eta_rad;
                for(int i=0;i<par.RsN_rad;i++) par.Rs_rad[i]= par.Rs0_rad + i *(par.Rc_rad - par.Rs0_rad) / par.RsN_rad ; 
                parset[14]=1; parset[2]=1;}
            else if(key=="eta_ang"){
                tmp_eta_ang = std::atof(value.c_str()); parset[7] = 0; }
            else if(key=="Rc_ang"){
                par.Rc_ang = std::atof(value.c_str()); parset[8] = 1; }
            else if(key=="Rs0_ang"){
                par.Rs0_ang = std::atof(value.c_str()); parset[9] = 1; }
            else if(key=="Rsst_ang"){
                par.Rsst_ang = std::atof(value.c_str()); parset[10] = 1; }
            else if(key=="RsN_ang"){
                par.RsN_ang = std::atoi(value.c_str()); parset[11] = 1; 
                par.eta_ang = new float[par.RsN_ang];
                par.Rs_ang  = new float[par.RsN_ang];
                for(int i=0;i<par.RsN_ang;i++) par.eta_ang[i]=tmp_eta_ang;
                for(int i=0;i<par.RsN_ang;i++) par.Rs_ang[i]= par.Rs0_ang + i *(par.Rc_ang - par.Rs0_ang) / par.RsN_ang ; 
                parset[15]=1; parset[7]=1;}
            else if(key=="zeta"){
                tmp_zeta = std::atof(value.c_str()); parset[12] = 0; }
            else if(key=="ThetasN"){
                par.ThetasN = std::atoi(value.c_str()); parset[13] = 1; 
                par.zeta = new int[par.ThetasN];
                par.zeta_half = new float[par.ThetasN];
                par.Thetas = new float[par.ThetasN];
                for(int i=0;i<par.ThetasN;i++) par.zeta[i]=tmp_zeta; parset[12]=1;
                for(int i=0;i<par.ThetasN;i++) par.Thetas[i]= (0.5f+ i)*(M_PI/par.ThetasN); parset[16]=1;}  
            }//gversion = 0
            else if(gversion ==1 ){
            //First read allocation sizes
            if(key=="RsN_rad"){
                par.RsN_rad = std::atoi(value.c_str());
                par.eta_rad = new float[par.RsN_rad];
                par.twoeta_rad = new float[par.RsN_rad];
                par.Rs_rad = new float[par.RsN_rad]; parset[6]=1;}
            // Then param arrays
            else if(key=="eta_rad"){
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.eta_rad[s] = std::atof(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.eta_rad[s] = std::atof(value.c_str()); s++; }; parset[2] = 1; }
            // Then cutoff
            else if(key=="Rc_rad"){
                par.Rc_rad = std::atof(value.c_str()); parset[3] = 1; }
            // Then the bin center arrays
            else if(key=="Rs_rad"){
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.Rs_rad[s] = std::atof(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.Rs_rad[s] = std::atof(value.c_str()); s++; }; parset[14] = 1;
                par.Rs0_rad=par.Rs_rad[0];                           parset[4]=1;}

            else if(key=="RsN_ang"){
                par.RsN_ang = std::atoi(value.c_str());
                par.eta_ang = new float[par.RsN_ang];
                par.Rs_ang = new float[par.RsN_ang]; parset[11]=1;}
            else if(key=="eta_ang") {
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.eta_ang[s] = std::atof(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.eta_ang[s] = std::atof(value.c_str()); s++; }; parset[7] = 1; }

            // Then cutoffs
            else if(key=="Rc_ang"){
                par.Rc_ang = std::atof(value.c_str()); parset[8] = 1; }
            // Then the bin center arrays
            else if(key=="Rs_ang") {
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.Rs_ang[s] = std::atof(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.Rs_ang[s] = std::atof(value.c_str()); s++; }; parset[15] = 1;
                par.Rs0_ang=par.Rs_ang[0];                           parset[9]=1;}

            else if(key=="ThetasN"){
                par.ThetasN = std::atoi(value.c_str());
                par.zeta = new int[par.ThetasN];
                par.zeta_half = new float[par.ThetasN];
                par.Thetas = new float[par.ThetasN]; parset[13]=1;}
            else if(key=="zeta") {
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.zeta[s] = std::atoi(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.zeta[s] = std::atoi(value.c_str()); s++; }; parset[12] = 1; }

            else if(key=="Thetas") {
                size_t pos = 0; int s = 0;
                value=value.substr(1, value.size() - 2); //get rid of [ ]
                while ((pos = value.find(comma)) != std::string::npos) {
                par.Thetas[s] = std::atof(value.substr(0, pos).c_str());
                value.erase(0, pos+1);  s++; }
                if(value.length()>0){ par.Thetas[s] = std::atof(value.c_str()); s++; };  parset[16] = 1; }
            } //gversion = 1
        } //Section 0 (Parameter parsing) is finished.
        // Parse species network
        else if(section<par.Nspecies+1){
            int s=section-1;
            // Read species network
            if(key=="Nlayers"){
            par.Nlayers[s] = std::atoi(value.c_str());
            // This has the extra gvect size
            par.layers_size[s] = new int[par.Nlayers[s]+1];
            par.layers_size[s][0] = par.gsize;
            par.layers_size[s][1] = 0;
            par.layers_activation[s] = new int[par.Nlayers[s]];
            for(int i=0;i<par.Nlayers[s]-1;i++) par.layers_activation[s][i]=1;
            par.layers_activation[s][par.Nlayers[s]-1]=0;
            network[s] = new double*[2*par.Nlayers[s]];
            }
            else if(key=="sizes"){
            if(par.Nlayers[s]==-1){
                std::cout << "Sizes cannot be set before Nlayers." << std::endl;
                return -3;
            }
            std::string comma = ",";
            size_t pos = 0;
            int l = 0;
            // Parse layers list
            while ((pos = value.find(comma)) != std::string::npos) {
                if(l>par.Nlayers[s]-2){
                std::cout << "Layers list longer than Nlayers." << std::endl;
                return -3;
                }
                std::string lsize = value.substr(0, pos);
                par.layers_size[s][l+1] = std::atoi(lsize.c_str());
                value.erase(0, pos+1);
                l++;
            }
            if(value.length()>0){
                par.layers_size[s][l+1] = std::atoi(value.c_str());
                l++;
            };
            if(l<par.Nlayers[s]){
                std::cout << "Layers list shorter than Nlayers." << std::endl;
                return -3;
            }
            }
            else if(key=="activations"){
            if(par.Nlayers[s]==-1){
                std::cout << "Activations cannot be set before Nlayers." << std::endl;
                return -3;
            }
            std::string comma = ",";
            size_t pos = 0;
            int l = 0;
            // Parse layers list
            while ((pos = value.find(comma)) != std::string::npos) {
                if(l>par.Nlayers[s]-2){
                std::cout << "Activations list longer than Nlayers." << std::endl;
                return -3;
                }
                std::string lact = value.substr(0, pos);
                int actnum = std::atoi(lact.c_str());
                if (actnum!=0 && actnum!=1 && actnum!=3 && actnum!=4 ){
                std::cout << "Activations unsupported: " << actnum << std::endl;
                return -3;
                }
                par.layers_activation[s][l] = actnum;
                value.erase(0, pos+1);
                l++;
            }
            if(value.length()>0){
                int actnum = std::atoi(value.c_str());
                if (actnum!=0 && actnum!=1 && actnum!=3 && actnum!=4){
                std::cout << "Activations unsupported: " << actnum << std::endl;
                return -3;
                }
                par.layers_activation[s][l] = actnum;
                l++;
            };
            if(l<par.Nlayers[s]){
                std::cout << "Activations list shorter than Nlayers." << std::endl;
                return -3;
            }
            }
            else if(key=="file"){
            if(par.layers_size[s][1]==0){
                std::cout << "Layers sizes unset before filename for species " << par.species[s] << std::endl;
                return -3;
            }
            // Read filename and load weights
            wfile_string = dir_string+"/"+value;
            weights_file.open(wfile_string.c_str(), std::ios::binary);
            if(!weights_file.is_open()){
                std::cout << "Error reading weights file for " << par.species[s] << std::endl;
                return -3;
            }
            for(int l=0; l<par.Nlayers[s]; l++){
                // Allocate and read the right amount of data
                // Weights
                network[s][2*l] = new double[par.layers_size[s][l]*par.layers_size[s][l+1]];
                for(int i=0; i<par.layers_size[s][l]; i++) {
                for(int j=0; j<par.layers_size[s][l+1]; j++) {
                    float num;
                    weights_file.read(reinterpret_cast<char*>(&num), sizeof(float));
                    if(weights_file.eof()){
                    std::cout << "Weights file " << wfile_string << " is too small." << std::endl;
                    return -3;
                    }
                    network[s][2*l][j*par.layers_size[s][l]+i] = (double)num;
                }
                }
                // Biases
                network[s][2*l+1] = new double[par.layers_size[s][l+1]];
                for(int d=0; d<par.layers_size[s][l+1]; d++) {
                float num;
                weights_file.read(reinterpret_cast<char*>(&num), sizeof(float));
                if(weights_file.eof()){
                    std::cout << "Weights file " << wfile_string << " is too small." << std::endl;
                    return -3;
                }
                network[s][2*l+1][d] = (double)num;
                }
            }
            // Check if we're not at the end
            std::ifstream::pos_type fpos = weights_file.tellg();
            weights_file.seekg(0, std::ios::end);
            std::ifstream::pos_type epos = weights_file.tellg();
            if(fpos!=epos){
                std::cout << "Weights file " << wfile_string << " is too big." << std::endl;
                return -3;
            }
            weights_file.close();
            spset[section-1] = 1;
            }
        }
        else{
            return -3;
        }
        }
        else if(parseint==3){
        // Version information is read:
        if(key == "!version") {
            version = value ;
            std::cout << "Network version " << value << std::endl; }
        else if(key == "!gversion") {
            gversion = std::atoi(value.c_str()) ;
            std::cout << "Gvector version " << value << std::endl;}
        }
        // Get new line
        parseint = get_input_line(&params_file,&key,&value);
    }

    // Derived params - for both gvect types done here
    par.cutmax = par.Rc_rad>par.Rc_ang ? par.Rc_rad : par.Rc_ang;
    for(int i=0; i<par.RsN_rad; i++) {
        par.twoeta_rad[i] = 2.0*par.eta_rad[i];}
    for(int i=0; i<par.ThetasN; i++) {
        par.zeta_half[i] = 0.5f*par.zeta[i];}
    par.iRc_rad = M_PI/par.Rc_rad;
    par.iRc_rad_half = 0.5*par.iRc_rad;
    par.iRc_ang = M_PI/par.Rc_ang;
    par.iRc_ang_half = 0.5*par.iRc_ang;
    //par.Rsi_rad = new float[par.RsN_rad];
    //for(int indr=0; indr<par.RsN_rad; indr++) par.Rsi_rad[indr] = par.Rs0_rad + indr * par.Rsst_rad;
    par.Rsi_rad = par.Rs_rad;
    //par.Rsi_ang = new float[par.RsN_ang];
    //for(int indr=0; indr<par.RsN_ang; indr++) par.Rsi_ang[indr] = par.Rs0_ang + indr * par.Rsst_ang;
    par.Rsi_ang = par.Rs_ang;
    par.Thi_cos = new float[par.ThetasN];
    par.Thi_sin = new float[par.ThetasN];
    for(int indr=0; indr<par.ThetasN; indr++)  {
        //float ti = (indr + 0.5f) * M_PI / par.ThetasN;
        float ti = par.Thetas[indr];
        par.Thi_cos[indr] = cos(ti);
        par.Thi_sin[indr] = sin(ti);
    }
    for(int s=0;s<par.Nspecies;s++){
        if(spset[s]!=1){
        std::cout << "Species network undefined for " << par.species[s] << std::endl;
        return -4;
        }
    }

    // Precalculate gvect shifts for any species pair
    par.typsh = new int*[par.Nspecies];
    for(int s=0; s<par.Nspecies; s++){
        par.typsh[s] = new int[par.Nspecies];
        for(int ss=0; ss<par.Nspecies; ss++){
        if(s<ss) par.typsh[s][ss] = par.Nspecies*par.RsN_rad +
                    (s*par.Nspecies - (s*(s+1))/2 + ss) *
                    par.RsN_ang * par.ThetasN;
        else par.typsh[s][ss] = par.Nspecies*par.RsN_rad +
                    (ss*par.Nspecies - (ss*(ss+1))/2 + s) *
                    par.RsN_ang * par.ThetasN;
        }
    }
    params_file.close();
    delete[] spset;
    return(0);

}

//Function that loads all the Gvector and NN parameters and prints out 
//  success message 
void PairPANNAMP::coeff(const std::string& directory, const std::string& filename)
{

  //if (!allocated) {
  //  allocate();
  //}

  // We now expect a directory and the parameters file name (inside the directory) with all params
  //if (narg != 4) {
  //  error->all(FLERR,"Format of pair_coeff command is\npair_coeff * *  network_directory parameter_file\n");
  //}

  std::cout << "Loading PANNA pair parameters from " << directory << "/" << filename << std::endl;
  int gpout = get_parameters(directory, filename);
  if(gpout==0){
    std::cout << "Network loaded!" << std::endl;
  }
  else{
    std::cout << "Error " << gpout << " while loading network!" << std::endl;
    exit(1);
  }

  //for (int i=1; i<=atom->ntypes; i++) {
  //  for (int j=1; j<=atom->ntypes; j++) {
  //    cutsq[i][j] = par.cutmax * par.cutmax;
  //  }
}

// Radial gvect contribution (and derivative part)
double PairPANNAMP::Gradial_d(double rdiff, int indr, double *dtmp){
  double cent = rdiff - par.Rsi_rad[indr];
  double gauss = exp( - par.eta_rad[indr] * cent * cent);
  double fc = 0.5 * ( 1.0 + cos(rdiff * par.iRc_rad) );
  *dtmp = ( par.iRc_rad_half * sin(rdiff * par.iRc_rad) +
         par.twoeta_rad[indr] * fc * cent ) * gauss / rdiff;
  return gauss * fc;
}

double PairPANNAMP::Gangular_d(double rdiff1, double rdiff2, double cosijk, int Rsi, int Thi, double* dtmp){
  if(cosijk> 0.999999999) cosijk =  0.999999999;
  if(cosijk<-0.999999999) cosijk = -0.999999999;
  double epscorr = 0.001;
  double sinijk = sqrt(1.0 - cosijk*cosijk + epscorr * pow(par.Thi_sin[Thi], 2) );
  double iRij = 1.0/rdiff1;
  double iRik = 1.0/rdiff2;
  double Rcent = 0.5 * (rdiff1 + rdiff2) - par.Rsi_ang[Rsi];
  double fcrad = 0.5 * ( 1.0 + par.Thi_cos[Thi] * cosijk + par.Thi_sin[Thi] * sinijk );
  double fcij = 0.5 * ( 1.0 + cos(rdiff1 * par.iRc_ang) );
  double fcik = 0.5 * ( 1.0 + cos(rdiff2 * par.iRc_ang) );
  double mod_norm = pow( 0.5 * (1.0 + sqrt(1.0 + epscorr * pow(par.Thi_sin[Thi], 2) ) ), par.zeta[Thi]);
  double fact0 = 2.0 * exp( - par.eta_ang[Rsi] * Rcent * Rcent) * pow(fcrad, par.zeta[Thi]-1) / mod_norm;
  double fact1 = fact0 * fcij * fcik;
  double fact2 = par.zeta_half[Thi] * fact1 * ( par.Thi_cos[Thi] - par.Thi_sin[Thi] * cosijk / sinijk );
  double fact3 = par.iRc_ang_half * fact0 * fcrad;
  double G = fact1 * fcrad;
  dtmp[0] = -iRij * ( par.eta_ang[Rsi] * Rcent * G
            + fact2 * cosijk * iRij
            + fact3 * fcik * sin(rdiff1 * par.iRc_ang) );
  dtmp[1] = fact2 * iRij * iRik;
  dtmp[2] = -iRik * ( par.eta_ang[Rsi] * Rcent * G
            + fact2 * cosijk * iRik
            + fact3 * fcij * sin(rdiff2 * par.iRc_ang) );
  return G;
}


void PairPANNAMP::compute_gvect(std::vector<particle_t>& particles, int index, int numneigh, double *G, double* dGdx, int num_atoms)
{
    // int *mask = atom->mask; We don't have a mask object here yet,
    //    let's see if we need one 
    float posx = particles[index].x;
    float posy = particles[index].y;
    float posz = particles[index].z;
    std::string id = particles[index].atom;

    // Elements to store neigh list for angular part
    // We allocate max possible size, so we don't need to reallocate
    int nan = 0;
    int ang_neigh[numneigh];
    int ang_type[numneigh];
    double dists[numneigh];
    double diffx[numneigh];
    double diffy[numneigh];
    double diffz[numneigh];

    int counter_below_rc = 0;
    int num_ele_at_96 = 0;
    int angular_counter = 0;

    //Loop over all the neighbors
    for(int n=0; n < numneigh; n++){ //4; n++) { 

        int nind = particles[index].neighbors[n];
        std::string neig_chem_id = particles[nind].atom;
        double dx = particles[nind].x-posx;
        double dy = particles[nind].y-posy;
        double dz = particles[nind].z-posz;
        double Rij = sqrt(dx*dx+dy*dy+dz*dz);
        int id_index; //shall be used in three if statmnts

        //exclude  atoms that are not in all group 
        //   did not implement groups here, let's see later if 
        //   might affect anything
        //int igroup = group->find("all");
        //int groupbit = group->bitmask[igroup];

        //If neighbor is within cutoff radius
        if (Rij < par.Rc_rad){ //and mask[nind] & groupbit){
            if (index == 9) {
            std::cout << "\nAtom 9 at coordinates: " << posx << " "  << posy << " " << posz << "\n";
            std::cout << "Neightbor to Atom 9: " << nind << "\n";
            std::cout << "With Coordinate x: " << particles[nind].x << "\n";
            std::cout << "With Coordinate y: " << particles[nind].y << "\n";
            std::cout << "With Coordinate z: " << particles[nind].z << "\n";
            std::cout << "At a distance: " << Rij << "\n";
            }

            //counter_below_rc +=1;


            // Add all radial parts
            //Here, we need to match the type with a list of types
            // to find the start index at which we will be adding to
            // the G vector.

            // Get the size of the array with all chemical species list
            int length = par.Nspecies;

            // Loop over each index and get the matching chemical identity
            for (int i = 0; i < length; i++) {
                if (par.species[i] == neig_chem_id) {
                    id_index = i; 
                    //std::cout << id << " ";
                }
            }

            //Now, we know the index of the atom in the 
            //  chemical id's list, so we can appropriately
            //  add the contribution to the gvectors at the 
            //  correct place

            //Gives the index at which we will start adding 
            //  elements to the gvecto r according to the 
            //  neighbor chemical identity 
            int indsh = id_index*par.RsN_rad;

            //Loop over all the radial gaussians
            for(int indr=0; indr<par.RsN_rad; indr++){
                double dtmp; 
                // Getting the simple G and derivative part
                G[indsh+indr] += Gradial_d(Rij, indr, &dtmp);
                
                if (index == 0 && indsh+indr == 96) {
                    num_ele_at_96 +=1;
                }

                //std::cout << "The first element of G in index: " << indsh+indr << "\n";
                //std::cout << "The first element value is:  " << Gradial_d(Rij, indr, &dtmp) << "\n";
                //std::cout << "The first neighbor is atom of index: " << nind << "\n";

                // Filling all derivatives
                int indsh2 = (indsh+indr)*(numneigh+1)*3;
                double derx = dtmp*dx;
                double dery = dtmp*dy;
                double derz = dtmp*dz;
                dGdx[indsh2 + numneigh*3     ] += derx;
                dGdx[indsh2 + numneigh*3 + 1 ] += dery;
                dGdx[indsh2 + numneigh*3 + 2 ] += derz;
                dGdx[indsh2 + n*3     ] -= derx;
                dGdx[indsh2 + n*3 + 1 ] -= dery;
                dGdx[indsh2 + n*3 + 2 ] -= derz;
            }

        }
        // If within radial cutoff, store quantities 
        // so they don't have to be recomputed
        if (Rij < par.Rc_ang){
            ang_neigh[nan] = n;
            ang_type[nan] = id_index + 1;
            dists[nan] = Rij;
            diffx[nan] = dx;
            diffy[nan] = dy;
            diffz[nan] = dz;
            nan++; //Updated with # on angular neighbors
            if (index == 0 && n == 0) {
                angular_counter += 1;
            }
        }
    }

    if (index == 0) {
        std::cout << "Additions to posiiton 96: " << num_ele_at_96 << "\n";
        std::cout << "Angular atoms around first Mg " <<  angular_counter << "\n";
    }

    //if (index == 0) {
    //    std::cout << "Atom 0 neigbors below r_c: " << counter_below_rc << "\n";
    //}

    //Double for loop over all angular neighbors
    for(int n=0; n < nan; n++) {
        for(int m=0; m < nan; m++) { 
            if (ang_neigh[n] == ang_neigh[m] && ang_type[n] == ang_type[m] && dists[n] == dists[m] && diffx[n] == diffx[m]) continue;
            // Compute cosine
            double cos_ijk = (diffx[n]*diffx[m] + diffy[n]*diffy[m] + diffz[n]*diffz[m]) /
                            (dists[n]*dists[m]);
            // Gvect shift due to species
            int indsh = par.typsh[ang_type[n]-1][ang_type[m]-1];
            // Loop over all bins
            for(int Rsi=0; Rsi<par.RsN_ang; Rsi++){
                for(int Thi=0; Thi<par.ThetasN; Thi++){
                double dtmp[3];
                int indsh2 = Rsi * par.ThetasN + Thi;
                // Adding the G part and computing derivative
                G[indsh+indsh2] += Gangular_d(dists[n], dists[m], cos_ijk, Rsi, Thi, dtmp)/2;  //to avoid double
                if (indsh+indsh2 == 137 && index == 0) { 
                    std::cout << "We are going through the atom:" << index << "\n";
                    std::cout << "We are going through its neighbor 1:" << ang_neigh[n] << "\n";
                    std::cout << "Type of neighboor 1:" << ang_type[n] << "\n";
                    std::cout << "We are going through its neighbor 2:" << ang_neigh[m] << "\n";
                    std::cout << "Type of neighboor 2:" << ang_type[m] << "\n";   
                    std::cout << "This was added to the G vector at position 120:" << Gangular_d(dists[n], dists[m], cos_ijk, Rsi, Thi, dtmp) << "\n";                   
                }
                // Computing the derivative contributions
                double dgdxj = dtmp[0]*diffx[n] + dtmp[1]*diffx[m];
                double dgdyj = dtmp[0]*diffy[n] + dtmp[1]*diffy[m];
                double dgdzj = dtmp[0]*diffz[n] + dtmp[1]*diffz[m];
                double dgdxk = dtmp[1]*diffx[n] + dtmp[2]*diffx[m];
                double dgdyk = dtmp[1]*diffy[n] + dtmp[2]*diffy[m];
                double dgdzk = dtmp[1]*diffz[n] + dtmp[2]*diffz[m];
                // Filling all the interested terms
                int indsh3 = (indsh+indsh2)*(numneigh+1)*3;
                dGdx[indsh3 + ang_neigh[n]*3     ] += dgdxj;
                dGdx[indsh3 + ang_neigh[n]*3 + 1 ] += dgdyj;
                dGdx[indsh3 + ang_neigh[n]*3 + 2 ] += dgdzj;
                dGdx[indsh3 + ang_neigh[m]*3     ] += dgdxk;
                dGdx[indsh3 + ang_neigh[m]*3 + 1 ] += dgdyk;
                dGdx[indsh3 + ang_neigh[m]*3 + 2 ] += dgdzk;
                dGdx[indsh3 + numneigh*3     ] -= dgdxj + dgdxk;
                dGdx[indsh3 + numneigh*3 + 1 ] -= dgdyj + dgdyk;
                dGdx[indsh3 + numneigh*3 + 2 ] -= dgdzj + dgdzk;
                }
            }
        }
    }

    if (index == 9) {
        for (int i = 0; i < 360; i++) {
        std::cout << "Element " << i << " of G vector is: " << G[i] << " \n";
        }
    int debug = 0;
    }
    
}

double PairPANNAMP::compute_network(double *G, double *dEdG, int type){
// double PairPANNAMP::compute_network(double *G, double *dEdG, int type, int inum){
  // *1 layer input
  // *2 layer output
  double *lay1, *lay2, *dlay1, *dlay2;
  dlay1 = new double[par.layers_size[type][0]*par.gsize];

  // std::ofstream myfile;
  // std::cout << "computation: " << inum << std::endl;
  // myfile.open("atom_" + std::to_string(inum) +".dat");
  lay1 = G;
  // myfile << "=G=,";
  // for(int i = 0; i <par.gsize; i++){
  //   myfile << lay1[i]<< ',';
  // }
  // myfile << std::endl;

  for(int i=0; i<par.layers_size[type][0]*par.gsize; i++) dlay1[i] = 0.0;
  // dG_i/dG_i = 1
  for(int i=0; i<par.gsize; i++) dlay1[i*par.gsize+i] = 1.0;
  // Loop over layers
  for(int l=0; l<par.Nlayers[type]; l++){
    int size1 = par.layers_size[type][l];
    int size2 = par.layers_size[type][l+1];
    lay2 = new double[size2];
    dlay2 = new double[size2*par.gsize];
    for(int i=0; i<size2*par.gsize; i++) dlay2[i]=0.0;
    // Matrix vector multiplication done by hand for now...
    // We compute W.x+b and W.(dx/dg)
    for(int i=0; i<size2; i++){
      // a_i = b_i
      lay2[i] = network[type][2*l+1][i];
      for(int j=0;j<size1; j++){
        // a_i += w_ij * x_j
        lay2[i] += network[type][2*l][i*size1+j]*lay1[j];
        // lay2[i] += network[type][2*l][j*size2+i]*lay1[j];
        for(int k=0; k<par.gsize; k++)
          // da_i/dg_k += w_ij * dx_j/dg_k
          dlay2[i*par.gsize+k] += network[type][2*l][i*size1+j]*dlay1[j*par.gsize+k];
          // dlay2[i*par.gsize+k] += network[type][2*l][j*size2+i]*dlay1[j*par.gsize+k];
      }
    }

    // here dump of W.x+b
    // myfile << "=L" << std::to_string(l) << "=,";
    // for(int i = 0; i <size2; i++){
    //   myfile << lay2[i]<< ',';
    // }
    // myfile << std::endl;

    // Apply appropriate activation
    // Gaussian
    if(par.layers_activation[type][l]==1){
      for(int i=0; i<size2; i++){
        double tmp = exp(-lay2[i]*lay2[i]);
        for(int k=0; k<par.gsize; k++)
          dlay2[i*par.gsize+k] *= -2.0*lay2[i]*tmp;
        lay2[i] = tmp;
      }
    }
    // ReLU
    else if(par.layers_activation[type][l]==3){
      for(int i=0; i<size2; i++){
        if(lay2[i]<0){
          lay2[i] = 0.0;
          for(int k=0; k<par.gsize; k++) dlay2[i*par.gsize+k] = 0.0;
        }
      }
    }
    // Tanh
    else if(par.layers_activation[type][l]==4){
      for(int i=0; i<size2; i++){
        double tmp = tanh(lay2[i]);
        for(int k=0; k<par.gsize; k++)
          dlay2[i*par.gsize+k] *= (1 - tmp * tmp);
        lay2[i] = tmp;
      }
    }
    // Otherwise it's linear and nothing needs to be done

    if(l!=0) delete[] lay1;
    delete[] dlay1;
    lay1 = lay2;
    dlay1 = dlay2;
  }
  // myfile.close();
  for(int i=0;i<par.gsize;i++) dEdG[i]=dlay1[i];
  double E = lay1[0];
  delete[] lay1;
  delete[] dlay1;
  return E;
}

float PairPANNAMP::compute(std::vector<particle_t>& particles, int n_atoms) {

    //DID NOT INCLUDE FORCE CALCULATIONS

    double total=0.0;

    //The original code loops over the particles here 
    //  and starts the computation

    for (int i = 0; i < n_atoms; i++) {
        
        particle_t& pi = particles[i];

        //initiate the inputs for the g_vector calculation
        double G[par.gsize];
        double dEdG[par.gsize];
        //Now, for dGdX we need to know the number of neighbors 
        //double dGdx[par.gsize*(numneigh[myind]+1)*3];
        int numneigh = std::end(pi.neighbors) - std::begin(pi.neighbors); // / sizeof(pi.neighbors[0]);
        double dGdx[par.gsize*(numneigh+1)*3];

        //Populate arrays with zeros
        for(int j=0; j < par.gsize; j++){
            G[j] = 0.0;
            for(int k=0; k < (numneigh+1)*3; k++)
                dGdx[j*(numneigh+1)*3+k] = 0.0;
        }

        //Now we can call here compute g_vector
        compute_gvect(particles, pi.index, numneigh, G, dGdx, n_atoms);
        
        //Debug printing
        //if (i == 0) {
        //    for (int j = 0; j < 360; j++) {
        //    std::cout << "Element " << j << " of G vector is: " << G[j] << " \n";
        //    }
        //}


        //Now, we just need to get the particle type index 
        // to calculate the energy
        int id_index;
        int length = par.Nspecies;

        for (int n = 0; n < length; n++) {
            if (par.species[n] == pi.atom) {
                id_index = n; 
                //std::cout << id << " ";
            }
        }

        double E = compute_network(G,dEdG,id_index);
        std::cout << "The energy of atom " << i << " (" << pi.atom << ") is: " << E << " eV. \n";
        total += E;

    }
    
    return total;
}


//void PairPANNAMP::allocate()
//{
// This allocates the nytpes of atoms as well as the cutsq
//    array, not sure what this last ones does
//}


//Deconstructor implementation 

PairPANNAMP::~PairPANNAMP() {
    // destructor implementation
}


int main() {
    
    //Task 1. Load atoms from the poscar file
    std::vector<particle_t> particles;
    //Initiate integer for initial number of particles and reserve memory
    //   for all ghost particles
    int n_atoms;
    particles.reserve(99999999999);

    double a, b, c;

    read("CO2-dobpdc_bare.vasp", particles, a, b, c, n_atoms);

    // Debug Task 1. Print out the atoms and their positions
    //for (int i = 0; i < particles.size(); i++) {
    //    std::cout << particles[i].atom << " ";
    //    std::cout << particles[i].x << " ";
    //    std::cout << particles[i].y << " ";
    //    std::cout << particles[i].z << std::endl;
    //}
    //return 0;

    //Task 2. Find nearest neighbors list

    //Save the original number of particles, so this is the for loop
    // max over energy computation, not taking the ghost particles
    // as well. 
    double rc = 6;  // Cutoff radius in angstroms

    PairPANNAMP pairPannampObj; // create an object of the PairPANNAMP class

    //pairPannampObj.fill_ghost_particles(particles, rc, a, b, c, n_atoms);
    find_neighbors(particles, rc, a, b, c, n_atoms);

    //Task 3. Load all the panna.in inputs
    std::string filename = "panna.in";
    std::string directory = "/Users/pedrogm/vs-docs/gcmc_serial/parameters_194k_steps_panna";
    //This executes the function coeff and populates the object
    pairPannampObj.coeff(directory, filename);
    
    //Task4. Do a single point calculation with PANNA from scratch
    float energy = pairPannampObj.compute(particles, n_atoms);

    int debug = 0; 

}