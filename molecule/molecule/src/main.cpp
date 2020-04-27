/**
c++ molecule_meanfield_hubbard.cpp -o main && ./main clars_goblet.xyz -t 1.8 -b 1.9
**/
#include <algorithm> //std::sort
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip> //std::setw()
#include <ios>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>

#include "TBTK/Array.h"
#include "TBTK/AbstractHoppingAmplitudeFilter.h"
#include "TBTK/FileParser.h"
#include "TBTK/FileWriter.h"
#include "TBTK/Functions.h"
#include "TBTK/Model.h"
#include "TBTK/Index.h"
#include "TBTK/IndexTree.h"
#include "TBTK/ParameterSet.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Property/EigenValues.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Vector2d.h"

//#include "molecule_meanfield_hubbard.hpp"

using comperator_t = std::function<bool(double, double)>;
using bonds_t = std::vector<std::pair<int, int>>;
using ss_t = std::string::size_type;
using std::string;
using std::size_t;
using std::vector;
using std::complex;
using std::ofstream;
using namespace TBTK;

void print_help(bool full);
complex<double> H_U(const Index &toIndex, const Index &fromIndex);

std::complex<double> i(0, 1);
double U{ 0.0 };
double spin_and_site_resolved_density_tol{ 1e-5 };
double k_density_tolerance{ 1e-7 };
double k_target_density_per_site{ 1.0 };
double k_mixing_parameter{ 0.5 };
double k_temperature{ 0.001 };
int k_multiplicity{ 0 };
Model model;
Array<double> spin_and_site_resolved_density;
size_t k_num_atoms{ 0 };



ofstream scf_convergence;
int iteration_counter = 0;
int run = 9;


//USAGE:     DEBUG(debug_var++);
unsigned int debug_var = 0;
#define DEBUG(x) std::cout << "---------- Reached position: " << (x) << " ----------"<< std::endl;

#define PRINTVAR(x) std::cout << #x << " = " << (x) << std::endl;
void print_help(bool full);

struct Eigenstate{
    double eigenvalue;
    bool spin_channel;
    std::vector<complex<double>> eigenvector;

    friend std::ostream& operator<<(std::ostream& out, Eigenstate es){
        out << es.spin_channel << " " << es.eigenvalue << " ";
        for(auto val : es.eigenvector){
            out << "(" << val.real() << ", " << val.imag() << ") ";
        }
        return out;
    }
};

struct XYZCoordinate{
    float x;
    float y;
    float z;

    friend std::ostream& operator<<(std::ostream& out, XYZCoordinate &xyz) {
        out << std::setw(7) << xyz.x << "\ty: " << std::setw(7)
            << xyz.y << "\tz: " << std::setw(7) << xyz.z << "\n";
        return out;
    }
};


struct XYZFileLine{
    XYZCoordinate xyz;
    std::string name;

    friend std::istream& operator>>(std::istream& in, XYZFileLine &xyz_file) {
        in >> xyz_file.name >> xyz_file.xyz.x
           >> xyz_file.xyz.y >> xyz_file.xyz.z;
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, XYZFileLine &xyz_file) {
        out << "Atom: " << std::setw(2) << xyz_file.name << "\tx: "
            << std::setw(7) << xyz_file.xyz.x << "\ty: " << std::setw(7)
            << xyz_file.xyz.y << "\tz: " << std::setw(7) << xyz_file.xyz.z;
        return out;
    }
};


/** One atom consists of XYZFileLine (atom name and coordiantes)
 * and an unique id to distingush and access the atoms
**/
class Atom{
public:
    //cartesian grid coordinates are initialized to -1 to denote that
    //the coordiantes are yet to be set.
    //This can only be done once all the atoms of the molecule have been read
    //(done in the Molecule class)
    Atom() : id_(counter_++) {};

    friend std::ostream& operator<<(std::ostream& out, Atom &atom) {
        out << atom.name_coords_ << "\n";
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Atom &atom) {
        in >> atom.name_coords_.name >> atom.name_coords_.xyz.x
           >> atom.name_coords_.xyz.y >> atom.name_coords_.xyz.z;
        return in;
    }

    std::ostream& print_atom(std::ostream& out){
        out << "Id: " << std::setw(3) << id_ << "\tAtom: " << std::setw(7)
            << name_coords_ << std::setw(7);
        return out;
    }

    XYZFileLine name_coords_;

private:
    friend class Molecule;
    const int id_;
    static int counter_;
};
int Atom::counter_ = 0;

class Molecule{
public:
    void construct_molecule(std::istream& in){
        std::string n_atom;
        getline(in, n_atom);
        std::cout << "Number of atoms in file: " << n_atom << '\n';

        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        Atom atom;
        while (in >> atom) {
            if(atom.name_coords_.name != "C") continue;
            atoms_.push_back(atom);
            val_x.push_back(atom.name_coords_.xyz.x);
            val_y.push_back(atom.name_coords_.xyz.y);
        }
    }

    // double find_unit_cell_size(){
    //     double max_val=0, min_val=0;
    //
    //     comperator_t compFunctor = [](double elem1, double elem2){
	// 			return elem1 < elem2;
	// 	};
    //     std::sort(val_y.begin(), val_y.end(), compFunctor);
    //
    //     min_val = val_y[0];
    //     max_val = val_y[val_y.size()-1];
    //
    //     return max_val - min_val;
    // }

    size_t size(){
        return atoms_.size();
    }

    friend std::ostream& operator<<(std::ostream& out, const Molecule& atoms){
        for(auto i : atoms.atoms_){
            out << i;
        }
        return out;
    }

private:
    std::vector<Atom> atoms_;
    std::vector<double> val_x;
    std::vector<double> val_y;
    friend class Bonds;
    friend class Atom;
};

class Bonds{
public:

    void print_bonds(){
        int bonds_counter = 0;
        std::cout << "-------------------------------------------------------\n"
                  << "Bonds added: " << '\n';
        for(auto i : bonds_){
            std::cout << "{" << i.first << "," << i.second << "}  ";
            ++bonds_counter;
        }
        std::cout << "\nTotal bonds created: " << bonds_counter
                  << "\n-------------------------------------------------------\n";
    }

    void add_bonds(Molecule mol, double threshold){
        bond_threshold_ = threshold;
        std::cout << "Size of the considered molecule (no H): " << mol.size() << "\n";
        for(size_t i = 0; i < mol.size(); ++i){
            double x_1 = mol.atoms_[i].name_coords_.xyz.x;
            double y_1 = mol.atoms_[i].name_coords_.xyz.y;
            double z_1 = mol.atoms_[i].name_coords_.xyz.z;

            for (size_t j = i; j < mol.size(); j++){
                if(i == j) continue;
                double x_2 = mol.atoms_[j].name_coords_.xyz.x;
                double y_2 = mol.atoms_[j].name_coords_.xyz.y;
                double z_2 = mol.atoms_[j].name_coords_.xyz.z;
                double x_diff = x_1 > x_2 ? x_1 - x_2 : x_2 - x_1;
                double y_diff = y_1 > y_2 ? y_1 - y_2 : y_2 - y_1;
                double z_diff = z_1 > z_2 ? z_1 - z_2 : z_2 - z_1;
                double result = std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
                //static int counter_same_dist = 0;

                //if two atoms are inside the sum of their atomic radii
                //(plus a threshold that is provided by the user, default 1.3)
                //there is a bond between them
                if(result < bond_threshold_){
                    bonds_.push_back(std::make_pair(i, j));
                }
            }
        }
    }

    void print_hoppings(){
        int hoppings_counter = 0;
        for(auto i : bonds_){
            std::cout << "model << HoppingAmplitude(-t, {" << std::setw(3)
                      << i.first << ",s}, {" << std::setw(3)
                      << i.second << ",s})+HC;" << std::endl;
            ++hoppings_counter;
        }
        std::cout << "Amount of hoppings created: " << hoppings_counter << '\n';
    }

    bonds_t get_bonds(){
        return bonds_;
    }

    std::pair<int, int> get_bond_by_index(int i){
        return bonds_[i];
    }

    size_t size(){
        return bonds_.size();
    }

    bonds_t bonds_;
private:
    //adjacency list for all bonds
    double bond_threshold_;
    friend class Molecule;

};

Model create_hamiltonian(Molecule mol, Bonds bonds, complex<double> t, bool hubbard){
    Model model;
    for(int s = 0; s < 2; ++s){
        for(unsigned int i = 0; i < bonds.size(); ++i){
            model << HoppingAmplitude(
                -t,
                {bonds.get_bond_by_index(i).first, s},
                {bonds.get_bond_by_index(i).second, s}
            ) + HC;
        }
    }
    if(hubbard){
        for(unsigned int s = 0; s < 2; ++s){
            for(unsigned int i = 0; i < mol.size(); ++i){
                model << HoppingAmplitude(
                    H_U,
                    {i, s},
                    {i, s}
                );
            }
        }
    }
    model.construct();
	model.setTemperature(k_temperature);

    std::cout << "Model size in create_hamiltonian: " << model.getBasisSize() << std::endl;
    return model;
}

class Postprocessor {
public:
    Postprocessor(string eigenvector) : eigenvector_(eigenvector){};

    bool is_close(complex<double> num, complex<double> lim){
        return (num.real() > lim.real() - threshold && num.real() < lim.real() + threshold && num.imag() > lim.imag() - threshold && num.imag() < lim.imag() + threshold);
    }

    void add_entries(ss_t start_of_num, ss_t end_of_num, ss_t end_of_cplx, bool spin_is_up){
        while(start_of_num != string::npos){
            string number_1 { "" };
            string number_2 { "" };
            start_of_num = eigenvector_.find('(', end_of_cplx+1);
            end_of_num = eigenvector_.find(',', end_of_cplx+1);
            end_of_cplx = eigenvector_.find(')', end_of_cplx+1);
            if(start_of_num == string::npos){
                break;
            }
            if(spin_is_up){
                number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
                number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
            }
            start_of_num = eigenvector_.find('(', end_of_cplx+1);
            end_of_num = eigenvector_.find(',', end_of_cplx+1);
            end_of_cplx = eigenvector_.find(')', end_of_cplx+1);
            if(!spin_is_up){
                number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
                number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
            }
            complex<double> num {std::stod(number_1), std::stod(number_2)};
            if(spin_is_up){
                spin_up_.push_back(num);
            }
            else{
                spin_down_.push_back(num);
            }
        }
    }

    void process(){
        ss_t start_of_num = eigenvector_.find('(');
        ss_t eigenvalue_end = start_of_num-1;
        eigenvalue_ = eigenvector_.substr(0, eigenvalue_end);
        eigenvector_ = eigenvector_.erase(0, eigenvalue_end);
        start_of_num = eigenvector_.find('(');
        ss_t end_of_num = eigenvector_.find(',');
        ss_t end_of_cplx = eigenvector_.find(')');

        if(start_of_num != string::npos){
            auto number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
            //std::cout << "number_1: " << number_1 << '\n';
            auto number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
            //std::cout << "number_2: " << number_2 << '\n';

            complex<double> first_num {std::stod(number_1), std::stod(number_2)};

            start_of_num = eigenvector_.find('(', end_of_cplx+1);
            end_of_num = eigenvector_.find(',', end_of_cplx+1);
            end_of_cplx = eigenvector_.find(')', end_of_cplx+1);

            number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
            //std::cout << "number_1: " << number_1 << '\n';
            number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
            //std::cout << "number_2: " << number_2 << '\n';

            complex<double> second_num {std::stod(number_1), std::stod(number_2)};

            //std::cout << "The first two numbers: " << first_num << " " << second_num << '\n';

            bool spin_is_up = true;
            complex<double> zero {0.0,0.0};

            complex<double> temp_1 { std::stod(number_1) };
            complex<double> temp_2 { std::stod(number_2) };

            complex<double> sum_even { std::abs(temp_1) };
            complex<double> sum_odd { std::abs(temp_2) };
            ss_t temp_end_of_cplx { end_of_cplx };

            //std::cout << "eigenvector_: " << eigenvector_ << '\n';
            while (start_of_num != string::npos) {
                start_of_num = eigenvector_.find('(', end_of_cplx+1);
                end_of_num = eigenvector_.find(',', end_of_cplx+1);
                end_of_cplx = eigenvector_.find(')', end_of_cplx+1);
                if(start_of_num == string::npos){
                    break;
                }
                number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
                number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
                complex<double> odd_number {std::stod(number_1), std::stod(number_2)};

                start_of_num = eigenvector_.find('(', end_of_cplx+1);
                end_of_num = eigenvector_.find(',', end_of_cplx+1);
                end_of_cplx = eigenvector_.find(')', end_of_cplx+1);

                number_1 = eigenvector_.substr(start_of_num+1, end_of_num-start_of_num-1);
                number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);
                complex<double> even_number {std::stod(number_1), std::stod(number_2)};

                sum_even += std::abs(even_number);
                sum_odd += std::abs(odd_number);

            }

            //std::cout << "Sum even: " << sum_even << '\n';
            //std::cout << "Sum odd: " << sum_odd << '\n';
            start_of_num = eigenvector_.find('(', 4);
            end_of_num = eigenvector_.find(',', start_of_num);
            end_of_cplx = eigenvector_.find(')', end_of_num);
            bool odd_greater_than_even { std::abs(sum_odd.real()) > std::abs(sum_even.real()) };

            if(odd_greater_than_even){
                spin_up_.push_back(first_num);
            }
            else{
                spin_down_.push_back(second_num);
            }
            add_entries(start_of_num, end_of_num, end_of_cplx, odd_greater_than_even);
        }
    }

    void print_results() {
        std::cout << "Spin up: " << '\n';
        for(auto i : spin_up_){
            std::cout << i << '\n';
        }
        std::cout << "Spin down: " << '\n';
        for(auto i : spin_down_){
            std::cout << i << '\n';
        }
    }

    string stringyfy(){
        string eigenvector {""};
        if(spin_up_.size()){
            eigenvector += "su " + eigenvalue_;
            for(auto i : spin_up_){
                eigenvector += " (";
                eigenvector += std::to_string(i.real());
                eigenvector += ",";
                eigenvector += std::to_string(i.imag());
                eigenvector += ")";
            }
        }
        else{
            eigenvector += "sd " + eigenvalue_;
            for(auto i : spin_down_){
                eigenvector += " (" + std::to_string(i.real());
                eigenvector += ",";
                eigenvector += std::to_string(i.imag());
                eigenvector += ")";
            }
        }
        return eigenvector;
    }

private:
    string eigenvalue_;
    string eigenvector_;
    vector<complex<double>> spin_up_;
    vector<complex<double>> spin_down_;

    static double threshold;
};
double Postprocessor::threshold = 1e-9;


/** Start self consistet mean-field Hubbard **/

//Initialize the spin and site resolved density with random numbers between 0
//and 1.
void init_spin_and_site_resolved_density(Array<double>& spin_and_site_resolved_density){
	srand(time(nullptr));
	for(unsigned int spin = 0; spin < 2; spin++){
        if(spin){
            for(unsigned int site = 0; site < k_num_atoms; site++){
                spin_and_site_resolved_density[
                    {spin, site}
                ] = (rand()%100)/100.0; //0.0; //
            }
        }
        else{
            for(unsigned int site = 0; site < k_num_atoms; site++){
                spin_and_site_resolved_density[
                    {spin, site}
                ] = (rand()%100)/100.0; //1.0; //
            }
        }
    }
}


//Callback function responsible for returning the current value of H_U for the
//given indices. Since H_U is supposed to be diagonal, toIndex and fromIndex is
//assumed to be equal and only the fromIndex is used to determine the spin and
//site. Compare to the model specification where H_C only is passed to the
//model with both the 'to' and 'from' indices equal.
complex<double> H_U(const Index &toIndex, const Index &fromIndex){
	unsigned int spin = fromIndex[1];
	unsigned int site = fromIndex[0];
	return U*spin_and_site_resolved_density[{(spin + 1)%2, site}];
}

//Adjusts the models chemical potential to fix the density. The algorithm is
//implemented as a binary search. The density is calculated, and depending on
//whether the density is too small or too large, the chemical potential is
//increased or decreased, respectively. This is iterated with a doubled step
//length until the calculation overshoots, at which point the procedure
//continues as before, but now cutting the step in two each iteration. The
//procedure stops once a density is found that differes from the
//k_target_density_per_site by at most k_density_tolerance.
void fixDensity(PropertyExtractor::Diagonalizer &propertyExtractor){
    double step_length = 1;

	//Get the eigenvalues.
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

	//Perform binary search. The flags step_direction and has_overshot
	//indicates in which direction the previous step was taken and whether
	//the calculation has yet overshot. The initial step can be taken in
	//either direction, therefore step_direction is initialized to zero.
	//Once the step direction has been set in the first iteration of the
	//loop, the first change in step direction will indicate that the step
	//has overshot the target value. At this point it is time to start
	//halving the step size since we are in the vicinity of the target
	//density. Before the overshot occurs, the step size is instead doubled
	//each step to rapidly cause the overshot independently of the initial
	//step_length.
	int step_direction = 0;
	bool has_overshot = false;
	while(true){
		//Calculate the density per unit cell.
		double density_per_unit_cell = 0;
		for(int n = 0; n < model.getBasisSize(); n++){
			density_per_unit_cell
				+= Functions::fermiDiracDistribution(
					eigenValues(n),
					model.getChemicalPotential(),
					model.getTemperature()
				)/(model.getBasisSize()/4);
		}

		//Exit the loop if the target density is met within the given
		//tolerance.
		if(
			abs(density_per_unit_cell - 2*k_target_density_per_site)
			< k_density_tolerance
		){
            //std::cout << "Density per site: " << density_per_unit_cell << '\n';
            //std::cout << "Final chemical potetntial: " << model.getChemicalPotential() << '\n';
			break;
		}

		//Determine whether an overshot has occured and step the chemical
		//potential.
		if(density_per_unit_cell < 2*k_target_density_per_site){
			if(step_direction == -1){
				has_overshot = true;
            }
			step_direction = 1;
		}
		else{
			if(step_direction == 1){
				has_overshot = true;
            }
			step_direction = -1;
		}
		model.setChemicalPotential(
			model.getChemicalPotential() + step_direction*step_length
		);

		if(has_overshot){
			step_length /= 2.0;
        }
		else{
			step_length *= 2.0;
	    }
    }
}

//Callback function that is to be called each time the model Hamiltonian has
//been diagonalized. The function first fixes the density to the target
//density, then calculates the spin and site resolved density. The new spin and
//site resolved density is mixed with the previous values to stabilize the
//self-consistent calculation.
bool self_consistency_callback(Solver::Diagonalizer &solver){
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

    auto eigen_values = solver.getEigenValues();
    auto eigen_vectors = solver.getEigenVectors();

    std::vector<Eigenstate> eigenstates;
    const int basis_size = model.getBasisSize();

    for(int i = 0; i < basis_size; ++i){

        bool state {false};
        double eval {eigen_values[i]};
        std::vector<complex<double>> evec;

        complex<double> spin_up {0.0, 0.0};
        complex<double> spin_down {0.0, 0.0};

        for (size_t j = 0; j < basis_size; j+=2) {
            spin_down += std::abs(eigen_vectors[i*basis_size + j]);
        }
        for (size_t j = 1; j < basis_size; j+=2) {
            spin_up += std::abs(eigen_vectors[i*basis_size + j]);
        }
        bool up_greater_than_down { std::abs(spin_up.real()) > std::abs(spin_down.real()) };
        //std::cout << "domination spin:" << odd_greater_than_even << '\n';
        state = up_greater_than_down;

        if(up_greater_than_down){
            for (size_t j = 1; j < basis_size; j+=2) {
                evec.push_back(eigen_vectors[i*basis_size + j]);
            }
        }
        else{
            for (size_t j = 0; j < basis_size; j+=2) {
                evec.push_back(eigen_vectors[i*basis_size + j]);
            }
        }
        eigenstates.push_back({eval, state, evec});
    }

    // std::cout << "eigenstates: " << '\n';
    // std::cout << eigenstates.size() << '\n';
    // std::cout << eigenstates[0] << '\n';
    // std::cout << eigenstates[1] << '\n';
    // std::cout << eigenstates[2] << '\n';


    // Multiplicity m = 2*S + 1
    // m = 1 --> S = 0
    // m = 3 --> S = 1

    int atoms_per_spin_channel {basis_size / 4};

    //TODO: fix for even multiplicity (doublet, quartet, ...)

    int electrons_spin_up {atoms_per_spin_channel};
    int electrons_spin_down {atoms_per_spin_channel};

    int spin_change {(k_multiplicity-1)/2};
    electrons_spin_up += spin_change;
    electrons_spin_down -= spin_change;

    static int print_cap{ 0 };
    if(!print_cap++){
        std::cout << "\nelectrons spin up: " << electrons_spin_up << '\n';
        std::cout << "electrons spin down: " << electrons_spin_down << '\n';
    }

    int counter_up {0};
    int counter_down {0};
    std::vector<double> density_up(basis_size/2, 0.0);
    std::vector<double> density_down(basis_size/2, 0.0);

    for (size_t j = 0; j < basis_size; j++){
        if(!eigenstates[j].spin_channel){
            counter_down += 1;
            if(counter_down > electrons_spin_down){
                continue;
            }
            for (size_t k = 0; k < basis_size/2; k++){
                density_down[k] += std::norm(eigenstates[j].eigenvector[k]);
            }
        }
        else{
            counter_up += 1;
            if(counter_up > electrons_spin_up){
                continue;
            }
            for (size_t k = 0; k < basis_size/2; k++){
                density_up[k] += std::norm(eigenstates[j].eigenvector[k]);
            }
        }
    }

    // std::cout << "Current density:" << '\n';
    // for (size_t i = 0; i < basis_size/2; i++) {
    //     std::cout << curr_density[i] << '\n';
    // }


	if(!k_multiplicity){
        fixDensity(propertyExtractor);
    }

    Array<double> old_spin_and_site_resolved_density
		= spin_and_site_resolved_density;

	//Calculate the spin and site resolved density. Note that the k-indices
	//are summed over using the IDX_SUM_ALL specifier, while the density is
	//stored separately for each spin and site because of the IDX_ALL
	//specifier.
	// Property::Density density = propertyExtractor.calculateDensity({
	// 	{IDX_ALL, IDX_ALL}
	// });

    // for(unsigned int spin = 0; spin < 2; spin++){
	// 	for(unsigned int site = 0; site < k_num_atoms; site++){
    //         std::cout << "Density of spin: " << spin << " and site: " << site << " is: " << density({site, spin}) << '\n';
    //     }
    // }

    // static int print_count = 0;
    // if(print_count++ < 3){
    //     for (size_t spin_ = 0; spin_ < 2; spin_++) {
    //         double total_spin_density = 0;
    //         for (size_t site = 0; site < k_num_atoms; site++) {
    //             total_spin_density += spin_and_site_resolved_density[{spin_, site}];
    //
    //         }
    //         // std::cout << std::boolalpha;
    //         // std::cout << "Total spin for spin " << spin_ << " is " << total_spin_density << '\n';
    //     }
    // }


    //Update the spin and site resolved density. Mix with the previous
    //value to stabilize the self-consistent calculation.
    for(unsigned int site = 0; site < k_num_atoms; site++){
        spin_and_site_resolved_density[{0, site}]
            = k_mixing_parameter*spin_and_site_resolved_density[
                {0, site}
            ] + (1 - k_mixing_parameter)*density_down[site]; ///(model.getBasisSize()/k_num_atoms);

    }
    for(unsigned int site = 0; site < k_num_atoms; site++){
        spin_and_site_resolved_density[{1, site}]
    		= k_mixing_parameter*spin_and_site_resolved_density[
    			{1, site}
    		] + (1 - k_mixing_parameter)*density_up[site];///(model.getBasisSize()/k_num_atoms);

    }


    //
	// //Update the spin and site resolved density. Mix with the previous
	// //value to stabilize the self-consistent calculation.
	// for(unsigned int spin = 0; spin < 2; spin++){
	// 	for(unsigned int site = 0; site < k_num_atoms; site++){
    //         spin_and_site_resolved_density[{spin, site}]
	// 			= k_mixing_parameter*spin_and_site_resolved_density[
	// 				{spin, site}
	// 			] + (1 - k_mixing_parameter)*density({
	// 				static_cast<int>(site),
    //                 static_cast<int>(spin)
	// 			});///(model.getBasisSize()/k_num_atoms);
    //
	// 	}
	// }

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
    unsigned int max_spin {0}, max_site{0};
	double max_difference = 0;
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms; site++){
			double difference = abs(
				spin_and_site_resolved_density[{spin, site}]
				- old_spin_and_site_resolved_density[{spin, site}]
			);
			if(difference > max_difference){
				max_difference = difference;
                max_spin = spin; max_site = site;
            }
		}
	}

    scf_convergence << run << "," << iteration_counter++ << ",";
    scf_convergence << std::abs(max_difference - spin_and_site_resolved_density_tol) << std::endl;
    // std::cout << "Conv. difference:"
    //             << std::abs(max_difference - spin_and_site_resolved_density_tol)
    //             << "\tabs. value: "
    //             << spin_and_site_resolved_density[{max_spin, max_site}]
    //             << "\tlocation: " << max_spin << " " << max_site << '\n';

	//Return whether the spin and site resolved density has converged. The
	//self-consistent loop will stop once true is returned.
	if(max_difference > spin_and_site_resolved_density_tol)
		return false;
	else
		return true;
}

int main(int argc, char *argv[]){
    string periodicity_direction { "" };
    double periodicity_distance { 0.0 };
    complex<double> t { 1.0 };
    double threshold { 1.7 };
    bool hubbard { false };
    int multiplicity { 0 };

    scf_convergence.open("scf_convergence.txt");
    //scf_convergence << "run,iteration,error" << std::endl;

    if(argc == 1){
        std::cout << "You need to provide at least a file." << '\n';
        print_help(false); exit(0);
    }
    string file = argv[1];
    std::ifstream in(file);

    for(int i = 0; i < argc; ++i){
        try{
            if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--periodic")){
                periodicity_direction = argv[++i];
                periodicity_distance = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bond_threshold")){
                threshold = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-t") || !strcmp(argv[i], "--hopping_amplitude")){
                t = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-H") || !strcmp(argv[i], "--Hubbard")){
                U = std::stod(argv[++i]);
                hubbard = U;
            }
            if(!strcmp(argv[i], "-T") || !strcmp(argv[i], "--temperature")){
                k_temperature = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-M") || !strcmp(argv[i], "--multiplicity")){
                multiplicity = std::stod(argv[++i]);
                k_multiplicity = multiplicity;
                if(k_temperature != 0.001){
                    std::cout << "Setting temperature to 0." << '\n';
                }
                k_temperature = 0;
            }
            if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
                print_help(true); exit(0);
            }
        }
        catch(...){
            std::cout << "Invalid amount of parameters. Don't forget to provide a value for the flags you're setting!" << '\n';
            print_help(false); exit(0);
        }
    }
    if(multiplicity && k_temperature){
        std::cout << "Invalid parameter combination. Cannot use both temperature and multiplicity as inputs at the same time." << '\n';
        print_help(false); exit(0);
    }
    bool periodic = false;
    if(periodicity_direction != ""){
        periodic = true;
    }

    std::cout << "Variable values: " << '\n';
    PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
    PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(hubbard); //PRINTVAR(periodic);
    PRINTVAR(k_temperature); PRINTVAR(multiplicity);

    //Usage example
    Molecule example_mol;

    //Add molecules from xyz file
    example_mol.construct_molecule(in);
    //std::cout << example_mol;

    k_num_atoms = example_mol.size();

    Bonds bonds;

    //Create bonds in molecule
    bonds.add_bonds(example_mol, threshold);
    //bonds.print_bonds();
    //bonds.print_hoppings();

    Array<double> dummy_array({2, example_mol.size()});
    spin_and_site_resolved_density = dummy_array;

    init_spin_and_site_resolved_density(spin_and_site_resolved_density);

    model = create_hamiltonian(example_mol, bonds, t, hubbard);

    //Setup the solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
    if(hubbard){
        solver.setSelfConsistencyCallback(self_consistency_callback);
	    solver.setMaxIterations(1000);
    }
	//Run the solver. This will run a self-consistent loop where the
	//Hamiltonian first is diagonalized, and then the
	//self_consistency_callback is called. The procedure is repeated until
	//either the self-consistency callback returns true, or the maximum
	//number of iterations is reached.
	solver.run();

    if(hubbard) std::cout << "Final chemical potential: " << model.getChemicalPotential() << '\n';
	//Create PropertyExtractor
	PropertyExtractor::Diagonalizer pe(solver);

    auto eigen_vectors = solver.getEigenVectors();
	auto eigen_values = solver.getEigenValues();

	int basisSize = model.getBasisSize();
	ofstream myfile;
	myfile.open("eigenval_eigenvec.txt");
	std::cout << "Writing Eigenvalues and Eigenvectors to file" << '\n';
	for(int i = 0; i < basisSize; ++i){
		myfile << eigen_values[i] << " ";
		for(int j = 0; j < basisSize; ++j)
			myfile << eigen_vectors[i*basisSize + j] << " ";
		myfile << '\n';
	}
	std::cout << "Writing done." << '\n';
	myfile.close();

    std::ifstream read_back;
    read_back.open("eigenval_eigenvec.txt");


    ofstream processed_data;
    processed_data.open("processed_ev.txt");
    std::cout << "Writing processed data to file." << '\n';
    string file_line { "" };
    while(std::getline(read_back, file_line)){
        Postprocessor processor_down(file_line);
        processor_down.process();
        processed_data << processor_down.stringyfy() << "\n";
    }

    std::cout << "Writing done." << '\n';
    return 0;
}

void print_help(bool full){
    if(full){
        printf("/**********************************************************************/\n");
        printf("// Tight-Binding and Mean-Field Hubbard approximation                 //\n");
        printf("// Bachelor's thesis, Spring Semester 2019, ETH Zurich                //\n");
        printf("// Author: Robin Worreby                                              //\n");
        printf("// License: Use if you like, but give me credit.                      //\n");
        printf("/**********************************************************************/\n");
        printf("\n");
    }
    printf("Usage: \n./Application molecule.xyz [parameters] \nPossible parameters:\n");
    printf("  -p or --periodic \t\t- sets the periodicity direction and distance, \n\
            \t\t\t  two parameters needed [X, Y, Z] and [distance]. \n\
            \t\t\t  Currently only X-direction is supported.\n");
    printf("  -t or --hopping_amplitude \t- sets the Hamiltonian parameter \n\
            \t\t\t  t value (Hopping amplitude), default is 1.0\n");
    printf("  -b or --bond_threshold \t- sets the bond threshold in Ångström\n");
    printf("  -H or --Hubbard \t\t- sets the Hamiltonian parameter U, \n\
            \t\t\t  default is 0.0 \n");
    printf("  -T or --temperature \t\t- sets the temperature for the system in K, \n\
            \t\t\t  default is 0.01 K. \n\
            \t\t\t  This is mutually exclusive with the -M (multiplicity) parameter. \n");
    printf("  -M or --multiplicity \t\t- sets the multiplicity, \n\
            \t\t\t  i.e. the solution that is to be found. \n\
            \t\t\t  Default is no specific solution desired. \n\
            \t\t\t  This is mutually exclusive with the -T (temperature) parameter. \n");
    printf("  -h or --help \t\t\t- prints this help info.\n");
}

/**
Triplet
su -8.53539 (0.281844,0.000000) (0.091730,0.000000) (0.029534,0.000000) (0.023610,0.000000) (0.063844,0.000000) (0.201811,0.000000) (0.019072,0.000000) (0.011999,0.000000) (0.022535,0.000000) (0.366645,0.000000) (0.470127,0.000000) (0.257849,0.000000) (0.152495,0.000000) (0.171976,0.000000) (0.257327,0.000000) (0.427018,0.000000) (0.366900,0.000000) (0.055401,0.000000) (0.032083,0.000000) (0.049855,0.000000) (0.014944,0.000000) (0.012289,0.000000)
su -7.41155 (0.439671,0.000000) (0.169286,0.000000) (0.064506,0.000000) (0.055785,0.000000) (0.132906,0.000000) (0.356348,0.000000) (0.041514,0.000000) (0.023992,0.000000) (0.037839,0.000000) (0.321237,0.000000) (0.085459,0.000000) (0.049859,0.000000) (0.158763,0.000000) (-0.111327,0.000000) (-0.343315,0.000000) (-0.544927,0.000000) (-0.226454,0.000000) (0.063943,0.000000) (0.015828,0.000000) (-0.025974,0.000000) (0.014741,0.000000) (0.021806,0.000000)
sd -7.41006 (-0.190082,0.000000) (-0.081559,0.000000) (-0.059502,0.000000) (-0.097105,0.000000) (-0.210851,0.000000) (-0.365069,0.000000) (-0.104316,0.000000) (-0.150050,0.000000) (-0.124683,0.000000) (-0.136424,0.000000) (-0.242222,0.000000) (-0.440045,0.000000) (-0.405930,0.000000) (-0.258558,0.000000) (-0.106772,0.000000) (-0.066009,0.000000) (-0.100608,0.000000) (-0.193323,0.000000) (-0.165613,0.000000) (-0.143324,0.000000) (-0.194632,0.000000) (-0.256733,0.000000)
sd -6.84265 (-0.112053,-0.000000) (-0.052032,0.000000) (-0.036239,-0.000000) (-0.049165,0.000000) (-0.090310,-0.000000) (-0.182420,0.000000) (0.007929,0.000000) (0.343703,0.000000) (0.114147,-0.000000) (-0.085709,0.000000) (-0.141531,-0.000000) (-0.217312,0.000000) (-0.162479,-0.000000) (-0.098092,-0.000000) (-0.048033,0.000000) (-0.037826,-0.000000) (-0.062854,0.000000) (0.034239,0.000000) (0.164444,-0.000000) (0.024133,0.000000) (0.434598,0.000000) (0.687288,-0.000000)
su -6.25197 (-0.047263,0.000000) (-0.017782,0.000000) (0.001941,0.000000) (0.023720,0.000000) (0.071815,0.000000) (0.106079,0.000000) (0.084482,0.000000) (0.410365,0.000000) (0.239370,0.000000) (-0.159035,0.000000) (-0.144412,0.000000) (0.046379,0.000000) (0.147599,0.000000) (0.100514,0.000000) (0.089988,0.000000) (0.032852,0.000000) (-0.069747,0.000000) (0.221241,0.000000) (0.253020,0.000000) (0.108870,0.000000) (0.403843,0.000000) (0.604568,0.000000)
sd -5.86842 (0.235849,0.000000) (0.166337,0.000000) (0.178196,0.000000) (0.266977,0.000000) (0.412663,0.000000) (0.418487,0.000000) (0.187612,0.000000) (0.050489,0.000000) (0.083671,0.000000) (0.003949,0.000000) (-0.225590,0.000000) (-0.362283,0.000000) (0.037466,0.000000) (-0.351501,0.000000) (-0.195834,0.000000) (-0.131978,0.000000) (-0.143446,0.000000) (0.014555,0.000000) (-0.077028,0.000000) (-0.179410,0.000000) (-0.038267,0.000000) (0.015840,0.000000)
su -5.72157 (-0.272866,0.000000) (-0.154657,0.000000) (-0.090942,0.000000) (-0.105729,0.000000) (-0.217038,0.000000) (-0.439442,0.000000) (-0.059794,0.000000) (0.120136,0.000000) (0.008535,0.000000) (0.239312,0.000000) (0.514282,0.000000) (0.156165,0.000000) (-0.136919,0.000000) (-0.058868,0.000000) (-0.278769,0.000000) (-0.299518,0.000000) (0.153070,0.000000) (-0.036485,0.000000) (0.030011,0.000000) (-0.009459,0.000000) (0.127102,0.000000) (0.214931,0.000000)
su -5.07462 (0.297908,0.000000) (0.174186,0.000000) (0.070101,0.000000) (0.009731,0.000000) (-0.044088,0.000000) (-0.099448,0.000000) (-0.022609,0.000000) (0.156691,0.000000) (-0.029332,0.000000) (0.241228,0.000000) (-0.112360,0.000000) (-0.470823,0.000000) (-0.371844,0.000000) (-0.363218,0.000000) (-0.087810,0.000000) (0.271370,0.000000) (0.136673,0.000000) (-0.209039,0.000000) (-0.112363,0.000000) (-0.169170,0.000000) (0.101265,0.000000) (0.283249,0.000000)
sd -4.19318 (-0.276409,0.000000) (-0.225564,0.000000) (-0.145110,0.000000) (-0.037262,0.000000) (0.085764,0.000000) (0.063139,0.000000) (0.102361,0.000000) (-0.085824,0.000000) (0.121525,0.000000) (-0.356101,0.000000) (-0.427738,0.000000) (-0.014818,0.000000) (0.254966,0.000000) (0.159895,0.000000) (-0.001254,0.000000) (-0.162213,0.000000) (-0.315103,0.000000) (0.275024,0.000000) (0.286296,0.000000) (0.252357,0.000000) (0.049940,0.000000) (-0.237428,0.000000)
su -3.99763 (-0.233765,0.000000) (-0.153219,0.000000) (-0.028830,0.000000) (0.089080,0.000000) (0.231435,0.000000) (0.266318,0.000000) (0.141927,0.000000) (0.088181,0.000000) (0.172846,0.000000) (-0.267789,0.000000) (0.134603,0.000000) (-0.048947,0.000000) (0.212164,0.000000) (-0.415341,0.000000) (-0.407541,0.000000) (0.151619,0.000000) (0.374363,0.000000) (0.142653,0.000000) (-0.091461,0.000000) (-0.210084,0.000000) (-0.121512,0.000000) (-0.065123,0.000000)
sd -3.88522 (-0.028556,0.000000) (0.176800,0.000000) (0.338783,0.000000) (0.398167,0.000000) (0.249956,0.000000) (-0.161645,0.000000) (0.108716,0.000000) (-0.086565,0.000000) (-0.042199,0.000000) (-0.065468,0.000000) (-0.093431,0.000000) (-0.144071,0.000000) (-0.367635,0.000000) (0.352293,0.000000) (0.328320,0.000000) (0.217115,0.000000) (0.070348,0.000000) (-0.118582,0.000000) (0.137585,0.000000) (0.296172,0.000000) (0.084161,0.000000) (-0.064831,0.000000)
sd -3.48836 (0.348338,0.000000) (0.272090,0.000000) (0.089099,0.000000) (-0.133972,0.000000) (-0.267501,0.000000) (0.129161,0.000000) (-0.325340,0.000000) (-0.346643,0.000000) (-0.306408,0.000000) (0.161278,0.000000) (-0.071532,0.000000) (-0.128852,0.000000) (0.017025,0.000000) (-0.023838,0.000000) (-0.100363,0.000000) (-0.135471,0.000000) (-0.128477,0.000000) (0.016839,0.000000) (0.325567,0.000000) (0.200212,0.000000) (0.354473,0.000000) (-0.071238,0.000000)
su -2.84235 (-0.066700,0.000000) (-0.007921,0.000000) (0.056513,0.000000) (0.109467,0.000000) (0.145620,0.000000) (0.010444,0.000000) (0.130702,0.000000) (0.429959,0.000000) (0.170761,0.000000) (-0.018121,0.000000) (0.067743,0.000000) (0.159637,0.000000) (-0.075161,0.000000) (0.162761,0.000000) (0.216118,0.000000) (-0.115383,0.000000) (-0.141499,0.000000) (-0.265462,0.000000) (-0.528276,0.000000) (-0.184185,0.000000) (-0.415647,0.000000) (0.170512,0.000000)
sd -2.74334 (0.092597,0.000000) (-0.206941,0.000000) (-0.368193,0.000000) (-0.262224,0.000000) (0.091362,0.000000) (0.212160,0.000000) (0.137613,0.000000) (-0.055426,0.000000) (0.113420,0.000000) (0.118765,0.000000) (0.078473,0.000000) (-0.400962,0.000000) (-0.081751,0.000000) (-0.129872,0.000000) (0.213705,0.000000) (0.410123,0.000000) (0.365916,0.000000) (0.129024,0.000000) (0.209987,0.000000) (0.065075,0.000000) (0.116585,0.000000) (-0.158509,0.000000)
sd -2.62541 (-0.149542,0.000000) (0.030499,0.000000) (0.188828,0.000000) (0.201864,0.000000) (0.015464,0.000000) (-0.361586,0.000000) (0.173865,0.000000) (-0.040306,0.000000) (0.235670,0.000000) (0.137388,0.000000) (0.341435,0.000000) (0.049850,0.000000) (-0.024322,0.000000) (-0.302737,0.000000) (-0.342974,0.000000) (-0.132054,0.000000) (0.162112,0.000000) (0.295014,0.000000) (0.328312,0.000000) (0.021538,0.000000) (0.154858,0.000000) (-0.266698,0.000000)
su -2.61087 (0.054704,0.000000) (-0.029396,0.000000) (-0.089986,0.000000) (-0.124580,0.000000) (-0.129373,0.000000) (-0.253158,0.000000) (0.166178,0.000000) (0.254360,0.000000) (0.517382,0.000000) (0.290656,0.000000) (-0.096360,0.000000) (-0.170172,0.000000) (0.005230,0.000000) (-0.059876,0.000000) (0.062299,0.000000) (0.068193,0.000000) (-0.112246,0.000000) (0.429519,0.000000) (0.140635,0.000000) (0.042532,0.000000) (-0.253754,0.000000) (-0.337296,0.000000)
sd -2.39273 (-0.301747,0.000000) (-0.341886,0.000000) (-0.109167,0.000000) (0.216958,0.000000) (0.310037,0.000000) (0.107896,0.000000) (-0.068014,0.000000) (-0.402539,0.000000) (-0.402417,0.000000) (-0.130852,0.000000) (0.130259,0.000000) (0.136526,0.000000) (0.029677,0.000000) (-0.132328,0.000000) (-0.103512,0.000000) (0.010025,0.000000) (0.116378,0.000000) (-0.226576,0.000000) (-0.022203,0.000000) (-0.140322,0.000000) (0.336931,0.000000) (0.127222,0.000000)
su -2.359 (0.071810,0.000000) (-0.422201,0.000000) (-0.539166,0.000000) (-0.450075,0.000000) (-0.211339,0.000000) (0.270740,0.000000) (-0.146544,0.000000) (0.022425,0.000000) (-0.117157,0.000000) (0.155399,0.000000) (-0.108577,0.000000) (-0.007793,0.000000) (0.188535,0.000000) (-0.086146,0.000000) (0.055228,0.000000) (0.088367,0.000000) (-0.128195,0.000000) (-0.057441,0.000000) (-0.154717,0.000000) (-0.133407,0.000000) (-0.034875,0.000000) (0.130942,0.000000)
sd -1.78239 (-0.154724,0.000000) (0.177540,0.000000) (0.327977,0.000000) (0.123651,0.000000) (-0.241446,0.000000) (-0.026877,0.000000) (-0.242213,0.000000) (-0.034515,0.000000) (-0.032784,0.000000) (-0.302765,0.000000) (-0.173623,0.000000) (-0.053053,0.000000) (0.392788,0.000000) (-0.217901,0.000000) (0.191848,0.000000) (0.401208,0.000000) (0.232381,0.000000) (0.227345,0.000000) (-0.015113,0.000000) (-0.266239,0.000000) (0.021912,0.000000) (0.016990,0.000000)
su -1.71256 (0.235251,0.000000) (0.428955,0.000000) (0.136880,0.000000) (-0.240279,0.000000) (-0.480020,0.000000) (-0.141945,0.000000) (-0.243028,0.000000) (0.133119,0.000000) (-0.006569,0.000000) (-0.330435,0.000000) (-0.077956,0.000000) (0.255896,0.000000) (0.253061,0.000000) (-0.033168,0.000000) (-0.180844,0.000000) (0.069196,0.000000) (0.107140,0.000000) (0.101301,0.000000) (-0.123770,0.000000) (-0.100212,0.000000) (-0.152028,0.000000) (0.056528,0.000000)
su -0.605557 (-0.243459,0.000000) (0.011804,0.000000) (0.248860,0.000000) (0.229192,0.000000) (-0.015522,0.000000) (-0.168551,0.000000) (-0.074496,0.000000) (-0.178890,0.000000) (-0.103090,0.000000) (0.301507,0.000000) (-0.023683,0.000000) (0.156865,0.000000) (0.337934,0.000000) (-0.291558,0.000000) (0.085533,0.000000) (0.239450,0.000000) (-0.438757,0.000000) (0.160576,0.000000) (-0.106149,0.000000) (-0.344022,0.000000) (0.097517,0.000000) (0.109299,0.000000)
sd -0.0113227 (0.116228,0.000000) (0.223113,0.000000) (-0.044855,0.000000) (-0.234882,0.000000) (0.034556,0.000000) (-0.182419,0.000000) (0.415449,0.000000) (-0.258563,0.000000) (0.163296,0.000000) (-0.002676,0.000000) (-0.117375,0.000000) (0.215623,0.000000) (-0.054082,0.000000) (0.024881,0.000000) (0.238036,0.000000) (0.046418,0.000000) (-0.219403,0.000000) (-0.018026,0.000000) (-0.124736,0.000000) (-0.455428,0.000000) (0.415117,0.000000) (-0.112009,0.000000)

su 0.485843 (0.002217,0.000000) (-0.129361,0.000000) (-0.009112,0.000000) (0.124220,0.000000) (0.085367,0.000000) (-0.084415,0.000000) (0.001882,0.000000) (0.124967,0.000000) (-0.083131,0.000000) (0.211562,0.000000) (-0.275184,0.000000) (0.251996,0.000000) (-0.013920,0.000000) (0.223697,0.000000) (-0.602261,0.000000) (0.386660,0.000000) (-0.124416,0.000000) (-0.168087,0.000000) (0.030281,0.000000) (0.337812,0.000000) (-0.157452,0.000000) (0.028279,0.000000)
sd 0.815587 (-0.047818,0.000000) (0.368999,0.000000) (0.052850,0.000000) (-0.371319,0.000000) (0.044590,0.000000) (0.135618,0.000000) (0.219655,0.000000) (-0.008617,0.000000) (-0.007255,0.000000) (-0.505614,0.000000) (-0.014010,0.000000) (0.144043,0.000000) (-0.110199,0.000000) (-0.017824,0.000000) (-0.363162,0.000000) (0.020270,0.000000) (0.365091,0.000000) (-0.214985,0.000000) (-0.001826,0.000000) (0.225844,0.000000) (-0.011154,0.000000) (0.011603,0.000000)
su 1.00723 (-0.218859,0.000000) (-0.113906,0.000000) (0.234783,0.000000) (0.201027,0.000000) (-0.150200,0.000000) (0.149822,0.000000) (-0.395192,0.000000) (0.171337,0.000000) (-0.242962,0.000000) (0.224946,0.000000) (-0.114816,0.000000) (-0.332890,0.000000) (0.209387,0.000000) (0.056118,0.000000) (0.103505,0.000000) (-0.181002,0.000000) (0.271618,0.000000) (0.150248,0.000000) (0.064245,0.000000) (0.215424,0.000000) (-0.352041,0.000000) (0.134670,0.000000)
sd 1.44852 (0.005466,0.000000) (-0.094412,0.000000) (0.015379,0.000000) (0.090133,0.000000) (-0.060160,0.000000) (0.236945,0.000000) (-0.291327,0.000000) (0.539502,0.000000) (0.078934,0.000000) (-0.143700,0.000000) (0.010648,0.000000) (0.150251,0.000000) (-0.199025,0.000000) (0.005002,0.000000) (0.012850,0.000000) (-0.008101,0.000000) (-0.011722,0.000000) (-0.223731,0.000000) (0.048406,0.000000) (-0.166161,0.000000) (0.386358,0.000000) (-0.477646,0.000000)
su 1.65253 (0.108527,0.000000) (-0.047924,0.000000) (-0.090373,0.000000) (0.035988,0.000000) (0.096914,0.000000) (-0.277012,0.000000) (0.246473,0.000000) (-0.399107,0.000000) (0.089386,0.000000) (0.169642,0.000000) (-0.400712,0.000000) (0.056930,0.000000) (0.155990,0.000000) (0.205345,0.000000) (-0.014619,0.000000) (-0.184212,0.000000) (0.440422,0.000000) (0.158350,0.000000) (-0.250897,0.000000) (-0.142474,0.000000) (-0.009148,0.000000) (0.258253,0.000000)
su 2.37243 (0.006084,0.000000) (-0.165248,0.000000) (0.100572,0.000000) (0.151715,0.000000) (-0.113448,0.000000) (0.014587,0.000000) (-0.142433,0.000000) (0.365886,0.000000) (0.043765,0.000000) (0.140333,0.000000) (-0.285204,0.000000) (0.161804,0.000000) (0.084443,0.000000) (0.045706,0.000000) (0.038311,0.000000) (-0.111301,0.000000) (0.248636,0.000000) (-0.232323,0.000000) (-0.058162,0.000000) (-0.234596,0.000000) (0.483986,0.000000) (-0.460022,0.000000)
sd 2.6771 (-0.307828,0.000000) (-0.107516,0.000000) (0.380489,0.000000) (-0.171510,0.000000) (-0.217234,0.000000) (0.290459,0.000000) (0.108994,0.000000) (-0.117141,0.000000) (0.160614,0.000000) (0.022875,0.000000) (0.294854,0.000000) (-0.379758,0.000000) (0.081873,0.000000) (0.259556,0.000000) (0.127485,0.000000) (-0.348309,0.000000) (0.079510,0.000000) (-0.015200,0.000000) (-0.240440,0.000000) (-0.024619,0.000000) (0.166779,0.000000) (-0.020739,0.000000)
su 2.81215 (0.494828,0.000000) (-0.402188,0.000000) (-0.169744,0.000000) (0.452673,0.000000) (0.057603,0.000000) (-0.258809,0.000000) (-0.215364,0.000000) (0.068793,0.000000) (-0.127894,0.000000) (-0.259585,0.000000) (0.063758,0.000000) (0.150979,0.000000) (-0.103596,0.000000) (-0.129432,0.000000) (0.079689,0.000000) (-0.019988,0.000000) (-0.024903,0.000000) (0.193320,0.000000) (0.141720,0.000000) (-0.111945,0.000000) (-0.146043,0.000000) (0.038427,0.000000)
sd 3.10825 (-0.167957,0.000000) (0.219993,0.000000) (-0.015846,0.000000) (-0.205842,0.000000) (0.244650,0.000000) (0.141344,0.000000) (-0.231280,0.000000) (-0.008903,0.000000) (-0.087573,0.000000) (-0.222218,0.000000) (0.329477,0.000000) (0.044011,0.000000) (-0.314929,0.000000) (-0.095316,0.000000) (0.317124,0.000000) (-0.176101,0.000000) (-0.184349,0.000000) (0.266897,0.000000) (0.323952,0.000000) (-0.244232,0.000000) (-0.245452,0.000000) (0.099626,0.000000)
sd 3.23521 (0.274637,0.000000) (-0.061547,0.000000) (-0.220320,0.000000) (0.268657,0.000000) (-0.090939,0.000000) (0.038676,0.000000) (-0.193113,0.000000) (-0.225077,0.000000) (0.231175,0.000000) (-0.217524,0.000000) (-0.106301,0.000000) (-0.003666,0.000000) (-0.250705,0.000000) (0.363906,0.000000) (-0.186694,0.000000) (-0.195341,0.000000) (0.343161,0.000000) (0.336802,0.000000) (-0.095430,0.000000) (-0.273073,0.000000) (0.006387,0.000000) (0.084108,0.000000)
sd 3.69433 (-0.080479,0.000000) (0.370807,0.000000) (-0.309820,0.000000) (-0.026880,0.000000) (0.345534,0.000000) (-0.325168,0.000000) (-0.140701,0.000000) (0.224928,0.000000) (-0.219434,0.000000) (0.038490,0.000000) (0.044147,0.000000) (-0.289715,0.000000) (0.353600,0.000000) (0.196817,0.000000) (0.046473,0.000000) (-0.246680,0.000000) (0.193062,0.000000) (0.030342,0.000000) (-0.149682,0.000000) (-0.040873,0.000000) (0.145961,0.000000) (-0.133887,0.000000)
su 3.70628 (0.064735,0.000000) (-0.324235,0.000000) (0.304715,0.000000) (0.132698,0.000000) (-0.381533,0.000000) (0.261429,0.000000) (-0.125366,0.000000) (-0.271436,0.000000) (0.382133,0.000000) (-0.079065,0.000000) (0.131583,0.000000) (-0.133948,0.000000) (-0.223156,0.000000) (0.286106,0.000000) (-0.170492,0.000000) (0.090033,0.000000) (-0.106100,0.000000) (0.130573,0.000000) (-0.262851,0.000000) (-0.052742,0.000000) (0.129155,0.000000) (0.060765,0.000000)
sd 4.64835 (0.168831,0.000000) (0.123866,0.000000) (-0.342976,0.000000) (0.378053,0.000000) (-0.292894,0.000000) (0.006463,0.000000) (0.136657,0.000000) (-0.070163,0.000000) (0.122132,0.000000) (-0.366475,0.000000) (0.306582,0.000000) (-0.052719,0.000000) (0.109482,0.000000) (-0.289244,0.000000) (0.326262,0.000000) (-0.176096,0.000000) (-0.093045,0.000000) (-0.173415,0.000000) (-0.081660,0.000000) (0.246194,0.000000) (0.029959,0.000000) (0.012871,0.000000)
su 4.71157 (0.173706,0.000000) (-0.261711,0.000000) (0.221944,0.000000) (0.039567,0.000000) (-0.259580,0.000000) (-0.110793,0.000000) (0.350730,0.000000) (-0.080814,0.000000) (0.127316,0.000000) (-0.072861,0.000000) (0.034336,0.000000) (-0.027344,0.000000) (0.355956,0.000000) (-0.340400,0.000000) (0.149289,0.000000) (-0.044546,0.000000) (0.004149,0.000000) (-0.406019,0.000000) (-0.009102,0.000000) (0.429758,0.000000) (-0.013182,0.000000) (0.034636,0.000000)
su 4.91668 (0.079163,0.000000) (0.128581,0.000000) (-0.283317,0.000000) (0.176515,0.000000) (0.102003,0.000000) (-0.229125,0.000000) (-0.064972,0.000000) (0.095263,0.000000) (-0.072566,0.000000) (-0.108437,0.000000) (0.238699,0.000000) (-0.457169,0.000000) (0.394781,0.000000) (0.235420,0.000000) (-0.113729,0.000000) (0.066474,0.000000) (-0.120287,0.000000) (0.052796,0.000000) (-0.387883,0.000000) (0.171457,0.000000) (0.255077,0.000000) (-0.125580,0.000000)
sd 5.0965 (-0.221007,0.000000) (0.225001,0.000000) (-0.132669,0.000000) (-0.008829,0.000000) (0.148985,0.000000) (0.071397,0.000000) (-0.352399,0.000000) (-0.167602,0.000000) (0.349852,0.000000) (0.049408,0.000000) (0.148711,0.000000) (-0.007866,0.000000) (-0.100894,0.000000) (-0.027589,0.000000) (-0.199249,0.000000) (0.344845,0.000000) (-0.314692,0.000000) (0.155654,0.000000) (-0.409391,0.000000) (0.261269,0.000000) (0.166091,0.000000) (0.000459,0.000000)
su 6.10097 (0.159565,0.000000) (-0.296664,0.000000) (0.441586,0.000000) (-0.372560,0.000000) (0.104517,0.000000) (-0.119904,0.000000) (0.326141,0.000000) (0.176185,0.000000) (-0.395339,0.000000) (-0.074651,0.000000) (0.092005,0.000000) (-0.187553,0.000000) (0.089912,0.000000) (0.256813,0.000000) (-0.094191,0.000000) (0.034532,0.000000) (-0.042524,0.000000) (0.123738,0.000000) (0.097245,0.000000) (-0.266641,0.000000) (-0.019923,0.000000) (-0.048402,0.000000)
sd 6.12976 (-0.252353,0.000000) (0.255009,0.000000) (-0.246082,0.000000) (0.240131,0.000000) (-0.289560,0.000000) (0.183351,0.000000) (0.250632,0.000000) (0.061967,0.000000) (-0.161136,0.000000) (0.053066,0.000000) (0.154398,0.000000) (-0.170981,0.000000) (0.027690,0.000000) (0.323031,0.000000) (-0.307436,0.000000) (0.284138,0.000000) (-0.224768,0.000000) (-0.083123,0.000000) (0.250930,0.000000) (-0.279270,0.000000) (-0.090985,0.000000) (0.007902,0.000000)
sd 6.52917 (0.415503,0.000000) (-0.286662,0.000000) (0.187206,0.000000) (-0.117708,0.000000) (0.092770,0.000000) (-0.194485,0.000000) (0.082496,0.000000) (0.113794,0.000000) (-0.253321,0.000000) (-0.389459,0.000000) (0.361023,0.000000) (-0.185716,0.000000) (0.065947,0.000000) (0.149157,0.000000) (-0.176567,0.000000) (0.225672,0.000000) (-0.279512,0.000000) (0.201944,0.000000) (-0.127924,0.000000) (-0.009638,0.000000) (0.057750,0.000000) (-0.044906,0.000000)
su 7.07529 (-0.063966,0.000000) (0.132015,0.000000) (-0.251184,0.000000) (0.339296,0.000000) (-0.368599,0.000000) (0.084921,0.000000) (0.295362,0.000000) (-0.017148,0.000000) (-0.001362,0.000000) (0.003065,0.000000) (0.052531,0.000000) (-0.172882,0.000000) (0.151207,0.000000) (0.263058,0.000000) (-0.083293,0.000000) (0.024635,0.000000) (-0.023128,0.000000) (-0.275567,0.000000) (0.413220,0.000000) (-0.400472,0.000000) (-0.164965,0.000000) (0.050738,0.000000)
sd 7.80752 (-0.129751,0.000000) (0.073289,0.000000) (-0.059039,0.000000) (0.082190,0.000000) (-0.175370,0.000000) (0.165619,0.000000) (0.269435,0.000000) (0.176099,0.000000) (-0.476565,0.000000) (0.094392,0.000000) (-0.103144,0.000000) (0.138254,0.000000) (-0.262287,0.000000) (-0.128917,0.000000) (0.069744,0.000000) (-0.052163,0.000000) (0.060373,0.000000) (0.529286,0.000000) (-0.338129,0.000000) (0.174491,0.000000) (0.117263,0.000000) (-0.068326,0.000000)
su 8.21976 (0.072089,0.000000) (-0.069956,0.000000) (0.124564,0.000000) (-0.216571,0.000000) (0.362840,0.000000) (-0.181232,0.000000) (-0.464334,0.000000) (-0.144392,0.000000) (0.415590,0.000000) (-0.027311,0.000000) (0.041376,0.000000) (-0.128588,0.000000) (0.242342,0.000000) (0.117986,0.000000) (-0.032957,0.000000) (0.009817,0.000000) (-0.013614,0.000000) (-0.375536,0.000000) (0.268583,0.000000) (-0.182984,0.000000) (-0.101957,0.000000) (0.061386,0.000000)

Singlet
su -8.8927 (0.357992,0.000000) (0.175986,0.000000) (0.096807,0.000000) (0.101104,0.000000) (0.225349,0.000000) (0.509104,0.000000) (0.118251,0.000000) (0.063782,0.000000) (0.186053,0.000000) (0.131854,0.000000) (0.120171,0.000000) (0.276726,0.000000) (0.508891,0.000000) (0.096115,0.000000) (0.028823,0.000000) (0.019089,0.000000) (0.040063,0.000000) (0.264073,0.000000) (0.081548,0.000000) (0.042281,0.000000) (0.024253,0.000000) (0.021629,0.000000)
sd -8.09034 (0.238405,0.000000) (0.106090,-0.000000) (0.029997,0.000000) (0.016144,-0.000000) (0.043254,0.000000) (0.131618,-0.000000) (0.032747,-0.000000) (0.029252,-0.000000) (0.062103,0.000000) (0.216329,-0.000000) (0.330451,0.000000) (0.474807,-0.000000) (0.200079,0.000000) (0.526122,-0.000000) (0.240160,-0.000000) (0.128439,-0.000000) (0.161803,-0.000000) (0.106931,-0.000000) (0.114542,0.000000) (0.285031,-0.000000) (0.032549,-0.000000) (0.016890,0.000000)
su -6.96783 (-0.471306,0.000000) (-0.380061,0.000000) (-0.239933,0.000000) (-0.135649,0.000000) (-0.095590,0.000000) (-0.206729,0.000000) (0.101522,0.000000) (0.175824,0.000000) (0.376414,0.000000) (-0.152738,0.000000) (0.026295,0.000000) (0.215465,0.000000) (0.270759,0.000000) (0.113976,0.000000) (0.040747,0.000000) (0.019837,0.000000) (0.016696,0.000000) (0.356846,0.000000) (0.147712,0.000000) (0.075002,0.000000) (0.061690,0.000000) (0.070747,0.000000)
sd -6.9319 (0.642294,0.000000) (0.358839,0.000000) (0.111583,0.000000) (0.047973,0.000000) (0.085507,0.000000) (0.237553,0.000000) (0.034623,0.000000) (0.001599,0.000000) (0.011033,0.000000) (0.351267,0.000000) (0.130681,0.000000) (-0.101936,0.000000) (0.039751,0.000000) (-0.353548,0.000000) (-0.177416,0.000000) (-0.053882,0.000000) (0.031906,0.000000) (-0.010944,0.000000) (-0.084646,0.000000) (-0.240945,0.000000) (-0.025505,0.000000) (-0.007401,0.000000)
su -6.21535 (-0.050274,0.000000) (-0.131594,0.000000) (-0.159313,0.000000) (-0.166431,0.000000) (-0.205965,0.000000) (-0.002396,0.000000) (-0.292734,0.000000) (-0.277277,0.000000) (-0.522194,0.000000) (0.069116,0.000000) (0.232385,0.000000) (0.438255,0.000000) (0.253475,0.000000) (0.228032,0.000000) (0.105833,0.000000) (0.090029,0.000000) (0.129777,0.000000) (-0.164258,0.000000) (-0.051392,0.000000) (0.055022,0.000000) (-0.048153,0.000000) (-0.105709,0.000000)
sd -5.75327 (-0.013084,0.000000) (-0.019389,0.000000) (-0.019186,0.000000) (-0.042184,0.000000) (-0.135705,0.000000) (-0.101779,0.000000) (-0.304895,0.000000) (-0.424739,0.000000) (-0.581335,0.000000) (0.107575,0.000000) (0.202847,0.000000) (0.084072,0.000000) (-0.135638,0.000000) (0.047117,0.000000) (0.081915,0.000000) (0.105239,0.000000) (0.156347,0.000000) (-0.348486,0.000000) (-0.209123,0.000000) (-0.117216,0.000000) (-0.128327,0.000000) (-0.197986,0.000000)
sd -5.11309 (-0.298588,0.000000) (-0.283788,0.000000) (-0.106448,0.000000) (-0.032597,0.000000) (-0.005510,0.000000) (-0.079586,0.000000) (0.095263,0.000000) (0.163841,0.000000) (0.206959,0.000000) (0.123968,0.000000) (0.487876,0.000000) (0.251592,0.000000) (0.100561,0.000000) (-0.305961,0.000000) (-0.119867,0.000000) (0.111437,0.000000) (0.345741,0.000000) (0.075643,0.000000) (-0.124435,0.000000) (-0.375886,0.000000) (-0.023718,0.000000) (0.054813,0.000000)
su -5.04168 (0.152416,0.000000) (0.369394,0.000000) (0.275339,0.000000) (0.025992,0.000000) (-0.229592,0.000000) (-0.439531,0.000000) (-0.001168,0.000000) (0.172632,0.000000) (0.227195,0.000000) (0.200560,0.000000) (0.288850,0.000000) (0.238163,0.000000) (-0.238889,0.000000) (0.234465,0.000000) (0.160382,0.000000) (0.177817,0.000000) (0.227678,0.000000) (0.049239,0.000000) (0.086250,0.000000) (0.115545,0.000000) (0.061865,0.000000) (0.088696,0.000000)
su -4.46457 (0.432909,0.000000) (0.008092,0.000000) (-0.425268,0.000000) (-0.527876,0.000000) (-0.391006,0.000000) (0.050990,0.000000) (-0.145805,0.000000) (0.117440,0.000000) (0.122871,0.000000) (0.218832,0.000000) (0.001784,0.000000) (-0.179278,0.000000) (-0.016134,0.000000) (-0.161433,0.000000) (-0.097154,0.000000) (-0.067547,0.000000) (-0.035820,0.000000) (0.121462,0.000000) (0.051213,0.000000) (-0.043022,0.000000) (0.045192,0.000000) (0.066924,0.000000)
su -4.04435 (0.066846,0.000000) (-0.341525,0.000000) (-0.336177,0.000000) (-0.017046,0.000000) (0.312470,0.000000) (0.191212,0.000000) (0.274822,0.000000) (0.118808,0.000000) (0.150157,0.000000) (0.182822,0.000000) (0.267863,0.000000) (0.028456,0.000000) (-0.312446,0.000000) (0.068055,0.000000) (0.144156,0.000000) (0.249268,0.000000) (0.307761,0.000000) (-0.303228,0.000000) (-0.184837,0.000000) (-0.048532,0.000000) (-0.065682,0.000000) (0.023357,0.000000)
sd -4.00743 (0.066126,0.000000) (0.111125,0.000000) (0.046970,0.000000) (0.009246,0.000000) (-0.019000,0.000000) (-0.129668,0.000000) (0.069862,0.000000) (0.209581,0.000000) (0.138126,0.000000) (0.044484,0.000000) (-0.016420,0.000000) (-0.426343,0.000000) (-0.325644,0.000000) (0.037977,0.000000) (0.442344,0.000000) (0.498731,0.000000) (0.364308,0.000000) (-0.112594,0.000000) (-0.038894,0.000000) (-0.001247,0.000000) (0.028508,0.000000) (0.110901,0.000000)
sd -3.37729 (0.014575,0.000000) (-0.090710,0.000000) (-0.085724,0.000000) (-0.108967,0.000000) (-0.218490,0.000000) (-0.207255,0.000000) (-0.214190,0.000000) (0.064415,0.000000) (-0.096747,0.000000) (0.300280,0.000000) (0.250877,0.000000) (-0.199245,0.000000) (-0.192886,0.000000) (-0.153601,0.000000) (-0.291914,0.000000) (-0.132458,0.000000) (0.108589,0.000000) (0.055487,0.000000) (0.388262,0.000000) (0.467332,0.000000) (0.238398,0.000000) (0.158253,0.000000)
su -2.83677 (-0.090223,0.000000) (0.043219,0.000000) (0.104976,0.000000) (0.021799,0.000000) (-0.084410,0.000000) (-0.173503,0.000000) (0.068169,0.000000) (0.044131,0.000000) (0.168674,0.000000) (0.126758,0.000000) (0.265597,0.000000) (0.225248,0.000000) (0.191555,0.000000) (-0.372093,0.000000) (-0.332736,0.000000) (-0.211524,0.000000) (0.043853,0.000000) (-0.086189,0.000000) (-0.420347,0.000000) (-0.404504,0.000000) (-0.270633,0.000000) (-0.123960,0.000000)
sd -2.69611 (0.039353,0.000000) (0.122016,0.000000) (0.025568,0.000000) (-0.068911,0.000000) (-0.200569,0.000000) (-0.244745,0.000000) (-0.122648,0.000000) (0.481068,0.000000) (0.051003,0.000000) (0.119054,0.000000) (0.035856,0.000000) (0.181323,0.000000) (-0.245617,0.000000) (0.251024,0.000000) (-0.026234,0.000000) (-0.270114,0.000000) (-0.279463,0.000000) (-0.321583,0.000000) (-0.295866,0.000000) (-0.179480,0.000000) (-0.004361,0.000000) (0.286968,0.000000)
su -2.51964 (0.092397,0.000000) (0.000744,0.000000) (-0.092230,0.000000) (-0.047036,0.000000) (0.053379,0.000000) (-0.004156,0.000000) (0.097749,0.000000) (0.229674,0.000000) (0.055968,0.000000) (-0.003830,0.000000) (-0.097246,0.000000) (0.372558,0.000000) (-0.144883,0.000000) (0.339055,0.000000) (-0.075771,0.000000) (-0.463057,0.000000) (-0.502246,0.000000) (-0.325332,0.000000) (-0.099796,0.000000) (0.129920,0.000000) (0.026385,0.000000) (0.149763,0.000000)
sd -2.14437 (-0.122539,0.000000) (0.493672,0.000000) (0.284325,0.000000) (0.038772,0.000000) (-0.193786,0.000000) (-0.073858,0.000000) (-0.346863,0.000000) (0.107793,0.000000) (-0.158324,0.000000) (-0.383330,0.000000) (-0.041287,0.000000) (0.201205,0.000000) (0.208646,0.000000) (-0.162686,0.000000) (-0.030324,0.000000) (0.146818,0.000000) (0.166481,0.000000) (0.157069,0.000000) (0.157144,0.000000) (-0.121828,0.000000) (0.201094,0.000000) (0.212026,0.000000)
su -2.00716 (0.087157,0.000000) (0.168961,0.000000) (-0.081396,0.000000) (-0.194365,0.000000) (-0.042257,0.000000) (-0.061494,0.000000) (0.227024,0.000000) (-0.343379,0.000000) (0.253129,0.000000) (-0.130841,0.000000) (-0.227978,0.000000) (-0.006523,0.000000) (-0.020008,0.000000) (0.247527,0.000000) (0.243964,0.000000) (0.105425,0.000000) (-0.132376,0.000000) (0.077762,0.000000) (-0.202774,0.000000) (0.027094,0.000000) (-0.409812,0.000000) (-0.495537,0.000000)
sd -1.96351 (-0.246202,0.000000) (0.399837,0.000000) (0.350453,0.000000) (0.232968,0.000000) (0.177959,0.000000) (-0.237798,0.000000) (0.343664,0.000000) (-0.149095,0.000000) (0.147883,0.000000) (-0.072245,0.000000) (0.220165,0.000000) (0.025412,0.000000) (-0.262518,0.000000) (0.041241,0.000000) (-0.248599,0.000000) (-0.154681,0.000000) (0.115512,0.000000) (-0.127881,0.000000) (-0.045712,0.000000) (0.207991,0.000000) (-0.145795,0.000000) (-0.212174,0.000000)
su -1.9308 (0.096474,0.000000) (-0.172148,0.000000) (-0.097475,0.000000) (0.144482,0.000000) (0.185306,0.000000) (-0.186539,0.000000) (0.163265,0.000000) (-0.421093,0.000000) (-0.038274,0.000000) (0.330085,0.000000) (0.249451,0.000000) (0.025377,0.000000) (-0.200994,0.000000) (-0.047389,0.000000) (-0.251898,0.000000) (-0.309916,0.000000) (-0.067370,0.000000) (0.264746,0.000000) (0.335098,0.000000) (0.177216,0.000000) (0.052523,0.000000) (-0.247086,0.000000)
sd -1.23645 (0.045414,0.000000) (0.173363,0.000000) (-0.046896,0.000000) (-0.245413,0.000000) (-0.443659,0.000000) (-0.344905,0.000000) (-0.134939,0.000000) (-0.177652,0.000000) (0.352054,0.000000) (0.142750,0.000000) (-0.032408,0.000000) (-0.019954,0.000000) (0.011381,0.000000) (0.027273,0.000000) (0.121489,0.000000) (-0.004550,0.000000) (-0.124179,0.000000) (0.376547,0.000000) (0.007308,0.000000) (-0.118929,0.000000) (-0.249083,0.000000) (-0.380821,0.000000)
sd -0.16546 (0.198378,0.000000) (-0.013911,0.000000) (-0.192741,0.000000) (-0.205762,0.000000) (-0.136939,0.000000) (0.163456,0.000000) (-0.127231,0.000000) (0.248994,0.000000) (0.101034,0.000000) (-0.354006,0.000000) (-0.090210,0.000000) (0.015992,0.000000) (0.057065,0.000000) (0.021795,0.000000) (-0.359530,0.000000) (0.053571,0.000000) (0.369948,0.000000) (-0.143485,0.000000) (-0.242458,0.000000) (0.320992,0.000000) (-0.364444,0.000000) (-0.159482,0.000000)
su -0.0851608 (0.176903,0.000000) (-0.470003,0.000000) (0.141644,0.000000) (0.413381,0.000000) (-0.172924,0.000000) (-0.019709,0.000000) (-0.388576,0.000000) (0.171068,0.000000) (0.088602,0.000000) (0.316341,0.000000) (-0.061622,0.000000) (-0.142834,0.000000) (0.018029,0.000000) (0.135221,0.000000) (0.208891,0.000000) (0.018289,0.000000) (-0.202552,0.000000) (0.140928,0.000000) (-0.151953,0.000000) (-0.017802,0.000000) (-0.243483,0.000000) (-0.089612,0.000000)

su 0.697249 (0.207198,0.000000) (0.027591,0.000000) (-0.233893,0.000000) (0.133686,0.000000) (0.185038,0.000000) (-0.468716,0.000000) (0.275957,0.000000) (-0.082237,0.000000) (-0.205122,0.000000) (0.178022,0.000000) (-0.193911,0.000000) (0.028247,0.000000) (0.266982,0.000000) (-0.099377,0.000000) (0.253534,0.000000) (0.212226,0.000000) (-0.241475,0.000000) (0.043010,0.000000) (-0.088156,0.000000) (-0.288448,0.000000) (0.201158,0.000000) (0.229438,0.000000)
sd 0.814497 (-0.255469,0.000000) (-0.082029,0.000000) (0.318479,0.000000) (0.329419,0.000000) (0.089765,0.000000) (0.153848,0.000000) (-0.404712,0.000000) (0.312104,0.000000) (-0.057086,0.000000) (0.284204,0.000000) (0.065478,0.000000) (-0.220273,0.000000) (0.221403,0.000000) (-0.050604,0.000000) (0.113278,0.000000) (-0.014255,0.000000) (-0.110876,0.000000) (0.125601,0.000000) (-0.136059,0.000000) (0.177709,0.000000) (-0.358831,0.000000) (-0.129453,0.000000)
su 1.73311 (0.230311,0.000000) (0.141876,0.000000) (-0.422012,0.000000) (0.311020,0.000000) (0.189027,0.000000) (-0.105927,0.000000) (-0.337961,0.000000) (0.138222,0.000000) (-0.034771,0.000000) (-0.416761,0.000000) (-0.101526,0.000000) (0.282838,0.000000) (-0.229720,0.000000) (-0.040665,0.000000) (-0.164429,0.000000) (0.030561,0.000000) (0.154441,0.000000) (0.253208,0.000000) (0.012540,0.000000) (-0.105536,0.000000) (-0.146184,0.000000) (-0.059128,0.000000)
sd 1.75296 (0.030074,0.000000) (-0.084929,0.000000) (0.064685,0.000000) (0.112692,0.000000) (0.035804,0.000000) (0.097411,0.000000) (-0.191216,0.000000) (-0.107439,0.000000) (0.046099,0.000000) (-0.064846,0.000000) (0.035815,0.000000) (-0.348223,0.000000) (-0.064469,0.000000) (0.523212,0.000000) (-0.199811,0.000000) (-0.339356,0.000000) (0.374943,0.000000) (0.255989,0.000000) (-0.013014,0.000000) (-0.364958,0.000000) (0.108181,0.000000) (0.055460,0.000000)
su 1.92108 (-0.102416,0.000000) (0.225106,0.000000) (-0.217417,0.000000) (0.023359,0.000000) (0.198292,0.000000) (-0.221178,0.000000) (0.044633,0.000000) (0.359671,0.000000) (-0.221772,0.000000) (0.172543,0.000000) (0.037085,0.000000) (-0.343582,0.000000) (0.315448,0.000000) (0.123172,0.000000) (-0.112495,0.000000) (-0.122253,0.000000) (0.160962,0.000000) (-0.047836,0.000000) (-0.042748,0.000000) (0.408508,0.000000) (-0.362768,0.000000) (-0.047607,0.000000)
su 2.2454 (0.370119,0.000000) (-0.294209,0.000000) (0.083236,0.000000) (0.189089,0.000000) (-0.260760,0.000000) (-0.171457,0.000000) (0.215132,0.000000) (-0.067522,0.000000) (0.121748,0.000000) (-0.216539,0.000000) (-0.262120,0.000000) (0.054907,0.000000) (0.230095,0.000000) (-0.050592,0.000000) (-0.311090,0.000000) (0.090501,0.000000) (0.264340,0.000000) (-0.357928,0.000000) (0.072221,0.000000) (0.281798,0.000000) (0.070995,0.000000) (-0.063066,0.000000)
sd 2.43376 (-0.171332,0.000000) (0.291575,0.000000) (-0.227509,0.000000) (-0.331857,0.000000) (0.015267,0.000000) (0.277000,0.000000) (0.059061,0.000000) (0.027957,0.000000) (-0.055456,0.000000) (-0.227052,0.000000) (0.459286,0.000000) (-0.285286,0.000000) (0.090227,0.000000) (-0.072404,0.000000) (0.369867,0.000000) (-0.361191,0.000000) (-0.092393,0.000000) (-0.021709,0.000000) (-0.026635,0.000000) (0.060020,0.000000) (-0.033205,0.000000) (0.021979,0.000000)
sd 2.71407 (0.316625,0.000000) (-0.262534,0.000000) (0.069746,0.000000) (0.267642,0.000000) (0.073641,0.000000) (0.009873,0.000000) (-0.264881,0.000000) (0.098742,0.000000) (0.134103,0.000000) (-0.411351,0.000000) (0.247768,0.000000) (0.167572,0.000000) (-0.393638,0.000000) (-0.151771,0.000000) (0.214290,0.000000) (-0.121688,0.000000) (-0.108174,0.000000) (-0.005713,0.000000) (0.262270,0.000000) (-0.063001,0.000000) (-0.008783,0.000000) (-0.262589,0.000000)
su 3.33312 (0.056101,0.000000) (0.011963,0.000000) (-0.079355,0.000000) (0.120224,0.000000) (-0.081950,0.000000) (-0.014622,0.000000) (0.000563,0.000000) (-0.312244,0.000000) (0.081359,0.000000) (-0.123347,0.000000) (0.055110,0.000000) (-0.273332,0.000000) (0.060689,0.000000) (0.405589,0.000000) (0.039138,0.000000) (-0.426377,0.000000) (0.352884,0.000000) (0.138359,0.000000) (-0.361710,0.000000) (-0.134554,0.000000) (0.167632,0.000000) (0.315795,0.000000)
su 3.42413 (-0.214823,0.000000) (0.297867,0.000000) (-0.374205,0.000000) (0.338091,0.000000) (-0.090807,0.000000) (0.060269,0.000000) (-0.277659,0.000000) (-0.192510,0.000000) (0.391438,0.000000) (0.131613,0.000000) (0.091724,0.000000) (-0.083126,0.000000) (0.159997,0.000000) (-0.090354,0.000000) (0.083079,0.000000) (0.043427,0.000000) (-0.124471,0.000000) (-0.376923,0.000000) (0.059684,0.000000) (0.085240,0.000000) (0.261383,0.000000) (-0.140089,0.000000)
sd 3.82294 (0.082965,0.000000) (-0.189289,0.000000) (0.273352,0.000000) (0.097046,0.000000) (-0.261217,0.000000) (0.084628,0.000000) (-0.119208,0.000000) (-0.429364,0.000000) (0.403669,0.000000) (-0.103402,0.000000) (0.101373,0.000000) (-0.047176,0.000000) (0.114595,0.000000) (-0.112800,0.000000) (0.094212,0.000000) (-0.046117,0.000000) (-0.035056,0.000000) (-0.134510,0.000000) (-0.398650,0.000000) (0.236277,0.000000) (0.179751,0.000000) (0.331369,0.000000)
sd 4.19028 (-0.270253,0.000000) (0.170759,0.000000) (-0.074415,0.000000) (-0.135523,0.000000) (0.075906,0.000000) (0.422454,0.000000) (-0.315409,0.000000) (-0.059810,0.000000) (0.343915,0.000000) (0.121307,0.000000) (0.037491,0.000000) (0.008907,0.000000) (-0.180896,0.000000) (0.122716,0.000000) (-0.238346,0.000000) (0.311756,0.000000) (-0.203966,0.000000) (-0.253538,0.000000) (0.096936,0.000000) (-0.095473,0.000000) (0.267368,0.000000) (-0.233388,0.000000)
sd 4.61948 (-0.068832,0.000000) (-0.073919,0.000000) (0.229785,0.000000) (-0.071412,0.000000) (-0.217647,0.000000) (0.098006,0.000000) (0.089661,0.000000) (0.131384,0.000000) (0.084052,0.000000) (0.168839,0.000000) (-0.281974,0.000000) (0.080120,0.000000) (0.183846,0.000000) (-0.100723,0.000000) (0.239556,0.000000) (-0.374034,0.000000) (0.350575,0.000000) (-0.388074,0.000000) (0.191687,0.000000) (-0.036983,0.000000) (0.233141,0.000000) (-0.347732,0.000000)
su 4.63839 (-0.131102,0.000000) (0.036299,0.000000) (0.042995,0.000000) (-0.128705,0.000000) (0.191908,0.000000) (0.004276,0.000000) (-0.216963,0.000000) (-0.251642,0.000000) (0.140580,0.000000) (0.317266,0.000000) (-0.308327,0.000000) (0.109981,0.000000) (-0.073062,0.000000) (0.118430,0.000000) (-0.470817,0.000000) (0.359253,0.000000) (-0.033165,0.000000) (0.101160,0.000000) (-0.277027,0.000000) (0.195910,0.000000) (-0.031846,0.000000) (0.301145,0.000000)
su 5.08904 (-0.069117,0.000000) (0.093618,0.000000) (-0.173740,0.000000) (0.308784,0.000000) (-0.441369,0.000000) (0.208136,0.000000) (0.341913,0.000000) (-0.112531,0.000000) (-0.139667,0.000000) (-0.101563,0.000000) (0.226737,0.000000) (0.009184,0.000000) (-0.120794,0.000000) (-0.129436,0.000000) (-0.037537,0.000000) (0.173785,0.000000) (-0.235265,0.000000) (0.158988,0.000000) (-0.095350,0.000000) (0.230209,0.000000) (-0.281995,0.000000) (0.355983,0.000000)
sd 5.18351 (-0.100816,0.000000) (0.216522,0.000000) (-0.415872,0.000000) (0.133379,0.000000) (0.365340,0.000000) (-0.080519,0.000000) (-0.324324,0.000000) (-0.151866,0.000000) (0.185654,0.000000) (0.167631,0.000000) (-0.282497,0.000000) (0.232754,0.000000) (-0.163384,0.000000) (-0.180414,0.000000) (0.199589,0.000000) (-0.256830,0.000000) (0.259275,0.000000) (0.068476,0.000000) (-0.117669,0.000000) (0.111703,0.000000) (-0.037789,0.000000) (0.150855,0.000000)
su 5.57262 (0.210226,0.000000) (-0.122360,0.000000) (0.129106,0.000000) (-0.199789,0.000000) (0.304664,0.000000) (-0.191286,0.000000) (-0.256319,0.000000) (-0.201817,0.000000) (0.176823,0.000000) (-0.332907,0.000000) (0.366054,0.000000) (-0.071431,0.000000) (0.099545,0.000000) (-0.270095,0.000000) (0.204549,0.000000) (-0.008214,0.000000) (-0.190185,0.000000) (-0.065225,0.000000) (-0.118715,0.000000) (0.336472,0.000000) (-0.116513,0.000000) (0.247270,0.000000)
sd 6.20819 (0.114180,0.000000) (0.082977,0.000000) (-0.343679,0.000000) (0.336614,0.000000) (0.088399,0.000000) (-0.376230,0.000000) (-0.059617,0.000000) (0.108599,0.000000) (0.035508,0.000000) (-0.093962,0.000000) (0.136339,0.000000) (-0.197620,0.000000) (0.412791,0.000000) (0.057627,0.000000) (-0.044521,0.000000) (0.056804,0.000000) (-0.078524,0.000000) (-0.140438,0.000000) (-0.199347,0.000000) (0.046495,0.000000) (0.410825,0.000000) (-0.317360,0.000000)
su 6.60796 (-0.085541,0.000000) (0.031390,0.000000) (-0.013548,0.000000) (0.007611,0.000000) (-0.005894,0.000000) (0.040426,0.000000) (-0.033252,0.000000) (-0.151843,0.000000) (0.081108,0.000000) (0.224067,0.000000) (-0.388251,0.000000) (0.243960,0.000000) (-0.053922,0.000000) (-0.319087,0.000000) (0.364790,0.000000) (-0.317125,0.000000) (0.311425,0.000000) (-0.086070,0.000000) (0.213857,0.000000) (0.068376,0.000000) (-0.343052,0.000000) (0.296196,0.000000)
sd 7.51161 (-0.004439,0.000000) (0.069007,0.000000) (-0.219734,0.000000) (0.305338,0.000000) (-0.159228,0.000000) (-0.039690,0.000000) (-0.010039,0.000000) (-0.149877,0.000000) (0.184939,0.000000) (-0.012120,0.000000) (0.042605,0.000000) (-0.102909,0.000000) (0.247745,0.000000) (0.075293,0.000000) (-0.031133,0.000000) (0.019754,0.000000) (-0.021193,0.000000) (-0.405694,0.000000) (0.482330,0.000000) (-0.157931,0.000000) (-0.435930,0.000000) (0.276395,0.000000)
su 7.66972 (-0.031238,0.000000) (0.011214,0.000000) (-0.008570,0.000000) (0.016827,0.000000) (-0.041033,0.000000) (0.044627,0.000000) (0.057608,0.000000) (0.130823,0.000000) (-0.111928,0.000000) (0.064495,0.000000) (-0.130499,0.000000) (0.156066,0.000000) (-0.105740,0.000000) (-0.312124,0.000000) (0.200332,0.000000) (-0.116045,0.000000) (0.092747,0.000000) (0.229787,0.000000) (-0.516220,0.000000) (0.428694,0.000000) (0.415313,0.000000) (-0.264590,0.000000)
sd 8.61075 (0.120258,0.000000) (-0.113623,0.000000) (0.295109,0.000000) (-0.509268,0.000000) (0.544273,0.000000) (-0.356899,0.000000) (-0.229125,0.000000) (-0.054687,0.000000) (0.135825,0.000000) (-0.044314,0.000000) (0.037324,0.000000) (-0.076694,0.000000) (0.236803,0.000000) (0.029592,0.000000) (-0.010858,0.000000) (0.007979,0.000000) (-0.013525,0.000000) (-0.186884,0.000000) (0.124954,0.000000) (-0.039246,0.000000) (-0.083686,0.000000) (0.054768,0.000000)

**/
