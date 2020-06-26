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
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Range.h"
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
complex<double> t{ 0.0 };
double spin_and_site_resolved_density_tol{ 1e-5 }; //5
double k_density_tolerance{ 1e-7 }; //7
double k_target_density_per_site{ 1.0 };
double k_mixing_parameter{ 0.5 };
double k_temperature{ 0.001 };
int k_multiplicity{ 0 };
int k_initial_guess{ 0101 };
Model model;
Array<double> spin_and_site_resolved_density;
size_t k_num_atoms{ 0 };
size_t k_num_unit_cells{ 0 };


ofstream scf_convergence;
int iteration_counter{ 0 };
int run{ 9 };

int verbosity{ 1 };
int k_num_iterations{ 1000 };

//USAGE:     DEBUG(debug_var++);
unsigned int debug_var = 0;
#define DEBUG(x) std::cout << "---------- Reached position: " << (x) \
                           << " ----------"<< std::endl;

#define PRINTVAR(x) std::cout << #x << " = " << (x) << std::endl;
void print_help(bool full);

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}

struct Eigenstate{
    double eigenvalue;
    bool spin_channel;
    double occupation;
    std::vector<complex<double>> eigenvector;

    friend std::ostream& operator<<(std::ostream& out, Eigenstate es){
        out << es.spin_channel << " " << es.eigenvalue << " "
            << es.occupation << " ";
        for(auto val : es.eigenvector){
            out << "(" << val.real() << ", " << val.imag() << ") ";
        }
        out << "\n";
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
        std::cout << "----------------------------------------------------\n"
                  << "Bonds added: " << '\n';
        for(auto i : bonds_){
            std::cout << "{" << i.first << "," << i.second << "}  ";
            ++bonds_counter;
        }
        std::cout << "\nTotal bonds created: " << bonds_counter
                  << "\n----------------------------------------------------\n";
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
                double result = std::sqrt(x_diff * x_diff + y_diff * y_diff
                                          + z_diff * z_diff);

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

    std::cout << "Model size in create_hamiltonian: "
              << model.getBasisSize() << std::endl;
    return model;
}

class Postprocessor {
public:
    Postprocessor(string eigenvector) : eigenvector_(eigenvector){};

    bool is_close(complex<double> num, complex<double> lim){
        return (num.real() > lim.real() - threshold && num.real() < lim.real()
                + threshold && num.imag() > lim.imag() - threshold && num.imag()
                < lim.imag() + threshold);
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
            number_2 = eigenvector_.substr(end_of_num+1, end_of_cplx-end_of_num-1);

            complex<double> second_num {std::stod(number_1), std::stod(number_2)};

            bool spin_is_up = true;
            complex<double> zero {0.0,0.0};

            complex<double> temp_1 { std::stod(number_1) };
            complex<double> temp_2 { std::stod(number_2) };

            complex<double> sum_even { std::abs(temp_1) };
            complex<double> sum_odd { std::abs(temp_2) };
            ss_t temp_end_of_cplx { end_of_cplx };

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
        static double total_energy_stringyfy{ 0.0 };
        static int total_energy_counter{ 0 };
        total_energy_counter++;
        total_energy_stringyfy += std::stod(eigenvalue_);
        if(total_energy_counter == 19){
            std::cout << "Total energy: " << total_energy_stringyfy << '\n';
        }

        if(spin_up_.size()){
            eigenvector += "1 " + eigenvalue_;
            for(auto i : spin_up_){
                eigenvector += " (";
                eigenvector += std::to_string(i.real());
                eigenvector += ",";
                eigenvector += std::to_string(i.imag());
                eigenvector += ")";
            }
        }
        else{
            eigenvector += "0 " + eigenvalue_;
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


double get_initial_guess(bool spin){
    switch(k_initial_guess){
        case 0:
            return 0;
            break;
        case 1:
            return 1;
            break;
        case 10:
            return spin;
            break;
        case 1010:
            return (rand()%100)/100.0;
            break;
        default:
            return (rand()%100)/100.0;
    }
}


/** Start self consistet mean-field Hubbard **/


void init_spin_and_site_resolved_density(Array<double>& spin_and_site_resolved_density){
	srand(time(nullptr));
	for(unsigned int spin = 0; spin < 2; spin++){
        if(spin){
            for(unsigned int site = 0; site < k_num_atoms; site++){
                spin_and_site_resolved_density[
                    {spin, site}
                ] = get_initial_guess(spin);
            }

        }
        else{
            for(unsigned int site = 0; site < k_num_atoms; site++){
                spin_and_site_resolved_density[
                    {spin, site}
                ] = get_initial_guess(spin);
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

		if(
			abs(density_per_unit_cell - 2*k_target_density_per_site)
			< k_density_tolerance
		){
			break;
		}

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

    if(k_temperature){
        fixDensity(propertyExtractor);
    }

    Array<double> old_spin_and_site_resolved_density
        = spin_and_site_resolved_density;

    //Calculate the spin and site resolved density. Note that the k-indices
    //are summed over using the IDX_SUM_ALL specifier, while the density is
    //stored separately for each spin and site because of the IDX_ALL
    //specifier.
    Property::Density density_temp = propertyExtractor.calculateDensity({
    	{IDX_ALL, IDX_ALL}
    });

    const int basis_size = model.getBasisSize();

    std::vector<Eigenstate> eigenstates;
    std::vector<double> density_up(basis_size/2, 0.0);
    std::vector<double> density_down(basis_size/2, 0.0);

    if(k_multiplicity){
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
            eigenstates.push_back({eval, state, 0.0, evec});
        }


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

        for (size_t j = 0; j < basis_size; j++){
            if(!eigenstates[j].spin_channel){
                counter_down += 1;
                if(counter_down > electrons_spin_down){
                    continue;
                }
                eigenstates[j].occupation = 1.0;
                for (size_t k = 0; k < basis_size/2; k++){
                    density_down[k] += std::norm(eigenstates[j].eigenvector[k]);
                }
            }
            else{
                counter_up += 1;
                if(counter_up > electrons_spin_up){
                    continue;
                }
                eigenstates[j].occupation = 1.0;
                for (size_t k = 0; k < basis_size/2; k++){
                    density_up[k] += std::norm(eigenstates[j].eigenvector[k]);
                }
            }
        }



        // std::cout << "Current density:" << '\n';
        // for (size_t i = 0; i < basis_size/2; i++) {
        //     std::cout << curr_density[i] << '\n';
        // }


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
    }
    else if(k_temperature){
    	// //Update the spin and site resolved density. Mix with the previous
    	// //value to stabilize the self-consistent calculation.
    	for(unsigned int spin = 0; spin < 2; spin++){
    		for(unsigned int site = 0; site < k_num_atoms; site++){
                spin_and_site_resolved_density[{spin, site}]
    				= k_mixing_parameter*spin_and_site_resolved_density[
    					{spin, site}
    				] + (1 - k_mixing_parameter)*density_temp({
    					static_cast<int>(site),
                        static_cast<int>(spin)
    			});///(model.getBasisSize()/k_num_atoms);

    		}
    	}
    }

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
	double max_difference{ 0 };
    double diff_total{ 0 };
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms; site++){
			double difference = abs(
				spin_and_site_resolved_density[{spin, site}]
				- old_spin_and_site_resolved_density[{spin, site}]
			);
            diff_total += difference;

			if(difference > max_difference){
				max_difference = difference;
            }
		}
	}

    // scf_convergence << run << "," << iteration_counter++ << ",";
    // scf_convergence << diff_total << std::endl;
    // std::cout << "Conv. difference:"
    //             << std::abs(max_difference - spin_and_site_resolved_density_tol)
    //             << "\tabs. value: "
    //             << spin_and_site_resolved_density[{max_spin, max_site}]
    //             << "\tlocation: " << max_spin << " " << max_site << '\n';

	//Return whether the spin and site resolved density has converged. The
	//self-consistent loop will stop once true is returned.
    if(max_difference > spin_and_site_resolved_density_tol)
		return false;
	else{
        if(k_multiplicity){
            double total_energy{ 0.0 };
            double magnetization{ 0.0 };
            std::cout << "Writing results to ev_with_occupations.txt" << '\n';
            std::ofstream afile("ev_with_occupations.txt", std::ios::out);
            if (afile.is_open()) {
                for (size_t i = 0; i < eigenstates.size(); i++) {
                    afile << eigenstates[i];
                    if(eigenstates[i].occupation){
                        total_energy += eigenstates[i].eigenvalue;
                        //spin_density = spin_density + eigenstates[i].eigenvector;
                    }
                }
                std::cout << '\n';
                for (size_t i = 0; i < density_down.size(); i++) {
                    magnetization += std::abs(density_down[i] - density_up[i]);
                }
                afile.close();
                std::cout << "Total Magnetization: " << magnetization << '\n';
                std::cout << "Total Energy: " << total_energy << '\n';
            }

            // std::ofstream scaling_results("7agnr_scaling.txt", std::ios::app);
            // if (scaling_results.is_open()) {
            //     scaling_results << k_num_unit_cells;
            //     scaling_results << ",";
            //     scaling_results << k_multiplicity;
            //     scaling_results << ",";
            //     scaling_results << (U/t.real());
            //     scaling_results << ",";
            //     scaling_results << magnetization;
            //     scaling_results << ",";
            //     scaling_results << total_energy;
            //     scaling_results << "\n";
            //     scaling_results.close();
            // }

        }
        return true;
    }
}


int main(int argc, char *argv[]){
    string periodicity_direction { "" };
    double periodicity_distance { 0.0 };
    double threshold { 1.7 };
    bool hubbard { false };
    int multiplicity { 0 };
    string ig{ "" };

    // scf_convergence.open("scf_convergence.txt");
    // scf_convergence << "run,iteration,error" << std::endl;

    if(argc == 1){
        std::cout << "You need to provide at least a file." << '\n';
        print_help(false); exit(0);
    }
    string file = argv[1];

    // Extracting the number of unit cells from the input file
    // string length{ "" };
    // for(int i = 24; i < file.size()-4; ++i){
    //     length += file[i];
    // }
    // k_num_unit_cells = std::stoi(length);
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
            if(!strcmp(argv[i], "-it") || !strcmp(argv[i], "--iterations")){
                k_num_iterations = std::stod(argv[++i]);
            }

            if(!strcmp(argv[i], "-M") || !strcmp(argv[i], "--multiplicity")){
                multiplicity = std::stod(argv[++i]);
                k_multiplicity = multiplicity;
                if(k_temperature != 0.001){
                    std::cout << "Setting temperature to 0." << '\n';
                }
                k_temperature = 0;
            }
            if(!strcmp(argv[i], "-ig") || !strcmp(argv[i], "--initial_guess")){
                ig = argv[++i];
                if(ig == "zero"){
                    k_initial_guess = 0;
                }
                else if(ig == "one"){
                    k_initial_guess = 1;
                }
                else if(ig == "random"){
                    k_initial_guess = 1010;
                }
                else if(ig == "zeroone"){
                    k_initial_guess = 10;
                }
                else{
                    k_initial_guess = 1010;
                }
            }
            if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")){
                verbosity = 2;
            }
            if(!strcmp(argv[i], "-q") || !strcmp(argv[i], "--quiet")){
                verbosity = 0;
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

    if(verbosity){
        std::cout << "Variable values: " << '\n';
        PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
        PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(U);
        PRINTVAR(k_temperature); PRINTVAR(k_multiplicity);
        std::cout << "initial_guess = " << ig << "\n\n";
    }


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

    if(verbosity == 2){
        bonds.print_hoppings();
    }

    Array<double> dummy_array({2, example_mol.size()});
    spin_and_site_resolved_density = dummy_array;

    init_spin_and_site_resolved_density(spin_and_site_resolved_density);

    model = create_hamiltonian(example_mol, bonds, t, hubbard);

    //Setup the solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
    if(hubbard){
        solver.setSelfConsistencyCallback(self_consistency_callback);
	    solver.setMaxIterations(k_num_iterations);
    }

	//Run the solver. This will run a self-consistent loop where the
	//Hamiltonian first is diagonalized, and then the
	//self_consistency_callback is called. The procedure is repeated until
	//either the self-consistency callback returns true, or the maximum
	//number of iterations is reached.
	solver.run();

    if(verbosity){
        if(hubbard && k_temperature) std::cout << "Final chemical potential: " \
                << model.getChemicalPotential() << '\n';
    }

	PropertyExtractor::Diagonalizer pe(solver);

    auto eigen_vectors = solver.getEigenVectors();
	auto eigen_values = solver.getEigenValues();

    for(int i = 0; i < model.getBasisSize(); ++i){
        std::cout << eigen_values[i] << " ";
    }

	int basisSize = model.getBasisSize();
	ofstream myfile;
	myfile.open("eigenval_eigenvec.txt");
	if(verbosity){
        std::cout << "Writing raw Eigenvalues and Eigenvectors to file eigenval_eigenvec.txt ..." << '\n';
	}
    for(int i = 0; i < basisSize; ++i){
		myfile << eigen_values[i] << " ";
		for(int j = 0; j < basisSize; ++j)
			myfile << eigen_vectors[i*basisSize + j] << " ";
		myfile << '\n';
	}
    if(verbosity > 0){
	   std::cout << "Writing done." << '\n';
    }

    myfile.close();

    std::ifstream read_back;
    read_back.open("eigenval_eigenvec.txt");

    ofstream processed_data;
    processed_data.open("processed_ev.txt");
    std::cout << "Writing processed data (eigenvalue, spin channel, occupation, eigenvector) to file processed_ev.txt ..." << '\n';
    string file_line { "" };
    while(std::getline(read_back, file_line)){
        Postprocessor processor_down(file_line);
        processor_down.process();
        processed_data << processor_down.stringyfy() << "\n";
    }
    if(verbosity > 0){
       std::cout << "Writing done." << '\n';
    }
    return 0;
}

void print_help(bool full){
    if(full){
        printf("/**********************************************************************/\n");
        printf("// Tight-Binding and Mean-Field Hubbard approximation                 //\n");
        printf("// Bachelor's thesis, 2019/2020, ETH Zurich                           //\n");
        printf("// Author: Robin Worreby, BSc RW/CSE                                  //\n");
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
    printf("  -it or --iterations \t\t- sets the max number of iterations for the self-consistent loop. \n");
    printf("  -ig or --initial_guess \t- lets you set the inital spin and site resolved density.\n\
            \t\t\t  Possibilities are: zero, one, zeroone, random. \n\
            \t\t\t  Where 'zeroone' means 1 on one spin channel and 0 on the other.\n\
            \t\t\t  Default is random.\n");
    printf("  -v or --verbose \t\t- prints additional information.\n");
    printf("  -q or --quiet \t\t- prints only the bare minimum.\n");
    printf("  -h or --help \t\t\t- prints this help info.\n");
}
