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
using cartesian_grid_t = std::multiset<int, std::pair<double, double>>;
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
double U { 0.0 };
double spin_and_site_resolved_density_tol { 1e-5 };
double DENSITY_TOLLERANCE { 1e-7 };
double TARGET_DENSITY_PER_SITE { 1.0 };
double MIXING_PARAMETER { 0.5 };
Model model;
Array<double> spinAndSiteResolvedDensity;
size_t k_num_atoms { 0 };

//USAGE:     DEBUG(debug_var++);
unsigned int debug_var = 0;
#define DEBUG(x) std::cout << "---------- Reached position: " << (x) << " ----------"<< std::endl;

#define PRINTVAR(x) std::cout << #x << " = " << (x) << std::endl;
void print_help(bool full);

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

    size_t size(){
        return bonds_.size();
    }

    bonds_t bonds_;
private:
    //adjacency list for all bonds
    double bond_threshold_;
    friend class Molecule;

};

Model create_hamiltonian(Molecule mol, Bonds bonds, complex<double> t_, bool hubbard){
    Model model;
    for(int s = 0; s < 2; ++s){
        for(unsigned int i = 0; i < bonds.bonds_.size(); ++i){
            model << HoppingAmplitude(-t_, {bonds.bonds_[i].first,s}, {bonds.bonds_[i].second,s})+HC;
        }
    }
    if(hubbard){
        for(unsigned int s = 0; s < 2; ++s){
            for(unsigned int i = 0; i < mol.size(); ++i){
                model << HoppingAmplitude(H_U, {i, s}, {i, s});
            }
        }
    }
    model.construct();
	model.setTemperature(0.01);

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
void initSpinAndSiteResolvedDensity(Array<double>& spinAndSiteResolvedDensity){
	srand(time(nullptr));
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms; site++){
			spinAndSiteResolvedDensity[
				{spin, site}
			] = 1.0; //(rand()%100)/100.;
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
	return U*spinAndSiteResolvedDensity[{(spin + 1)%2, site}];
}

//Adjusts the models chemical potential to fix the density. The algorithm is
//implemented as a binary search. The density is calculated, and depending on
//whether the density is too small or too large, the chemical potential is
//increased or decreased, respectively. This is iterated with a doubled step
//length until the calculation overshoots, at which point the procedure
//continues as before, but now cutting the step in two each iteration. The
//procedure stops once a density is found that differes from the
//TARGET_DENSITY_PER_SITE by at most DENSITY_TOLLERANCE.
void fixDensity(PropertyExtractor::Diagonalizer &propertyExtractor){
    double stepLength = 1;

	//Get the eigenvalues.
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

	//Perform binary search. The flags stepDirection and hasOvershot
	//indicates in which direction the previous step was taken and whether
	//the calculation has yet overshot. The initial step can be taken in
	//either direction, therefore stepDirection is initialized to zero.
	//Once the step direction has been set in the first iteration of the
	//loop, the first change in step direction will indicate that the step
	//has overshot the target value. At this point it is time to start
	//halving the step size since we are in the vicinity of the target
	//density. Before the overshot occurs, the step size is instead doubled
	//each step to rapidly cause the overshot independently of the initial
	//stepLength.
	int stepDirection = 0;
	bool hasOvershot = false;
	while(true){
		//Calculate the density per unit cell.
		double densityPerUnitCell = 0;
		for(int n = 0; n < model.getBasisSize(); n++){
			densityPerUnitCell
				+= Functions::fermiDiracDistribution(
					eigenValues(n),
					model.getChemicalPotential(),
					model.getTemperature()
				)/(model.getBasisSize()/4);
		}

		//Exit the loop if the target density is met within the given
		//tollerance.
		if(
			abs(densityPerUnitCell - 2*TARGET_DENSITY_PER_SITE)
			< DENSITY_TOLLERANCE
		){
            //std::cout << "Density per site: " << densityPerUnitCell << '\n';
            //std::cout << "Final chemical potetntial: " << model.getChemicalPotential() << '\n';
			break;
		}

		//Determine whether an overshot has occured and step the chemical
		//potential.
		if(densityPerUnitCell < 2*TARGET_DENSITY_PER_SITE){
			if(stepDirection == -1)
				hasOvershot = true;

			stepDirection = 1;
		}
		else{
			if(stepDirection == 1)
				hasOvershot = true;

			stepDirection = -1;
		}
		model.setChemicalPotential(
			model.getChemicalPotential() + stepDirection*stepLength
		);


		//Scale the stepLength depending on whether the overshot has
		//occurred or not.
		if(hasOvershot)
			stepLength /= 2.0;
		else
			stepLength *= 2.0;
	}
}

//Callback function that is to be called each time the model Hamiltonian has
//been diagonalized. The function first fixes the density to the target
//density, then calculates the spin and site resolved density. The new spin and
//site resolved density is mixed with the previous values to stabilize the
//self-consistent calculation.
bool selfConsistencyCallback(Solver::Diagonalizer &solver){
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

	//Fix the density.
	fixDensity(propertyExtractor);

	//Save the old result for later comparison.
	Array<double> oldSpinAndSiteResolvedDensity
		= spinAndSiteResolvedDensity;

	//Calculate the spin and site resolved density. Note that the k-indices
	//are summed over using the IDX_SUM_ALL specifier, while the density is
	//stored separately for each spin and site because of the IDX_ALL
	//specifier.
	Property::Density density = propertyExtractor.calculateDensity({
		{IDX_ALL, IDX_ALL}
	});

	//Update the spin and site resolved density. Mix with the previous
	//value to stabilize the self-consistent calculation.
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms; site++){
			spinAndSiteResolvedDensity[{spin, site}]
				= MIXING_PARAMETER*spinAndSiteResolvedDensity[
					{spin, site}
				] + (1 - MIXING_PARAMETER)*density({
					(int)site,
                    (int)spin
				});///(model.getBasisSize()/k_num_atoms);
		}
	}

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
    unsigned int max_spin {0}, max_site{0};
	double maxDifference = 0;
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms; site++){
			double difference = abs(
				spinAndSiteResolvedDensity[{spin, site}]
				- oldSpinAndSiteResolvedDensity[{spin, site}]
			);
			if(difference > maxDifference){
				maxDifference = difference;
                max_spin = spin; max_site = site;
            }

		}
	}

    //std::cout << "Conv. difference:"
                //<< std::abs(maxDifference - spin_and_site_resolved_density_tol)
                //<< "\tabs. value: "
                //<< spinAndSiteResolvedDensity[{max_spin, max_site}]
                //<< "\tlocation: " << max_spin << " " << max_site << '\n';

	//Return whether the spin and site resolved density has converged. The
	//self-consistent loop will stop once true is returned.
	if(maxDifference > spin_and_site_resolved_density_tol)
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
    //double h_u { 0.0 };

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
            if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
                print_help(true); exit(0);
            }
        }
        catch(...){
            std::cout << "Invalid amount of parameters. Don't forget to provide a value for the flags you're setting!" << '\n';
            print_help(false); exit(0);
        }
    }
    bool periodic = false;
    if(periodicity_direction != ""){
        periodic = true;
    }

    std::cout << "Variable values: " << '\n';
    PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
    PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(hubbard); PRINTVAR(periodic);

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
    spinAndSiteResolvedDensity = dummy_array;

    initSpinAndSiteResolvedDensity(spinAndSiteResolvedDensity);


    model = create_hamiltonian(example_mol, bonds, t, hubbard);

    //Setup the solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
    if(hubbard){
        solver.setSelfConsistencyCallback(selfConsistencyCallback);
	    solver.setMaxIterations(1000);
    }
	//Run the solver. This will run a self-consistent loop where the
	//Hamiltonian first is diagonalized, and then the
	//selfConsistencyCallback is called. The procedure is repeated until
	//either the self-consistency callback returns true, or the maximum
	//number of iterations is reached.
	solver.run();

    if(hubbard) std::cout << "Final chemical potential: " << model.getChemicalPotential() << '\n';
	//Create PropertyExtractor
	PropertyExtractor::Diagonalizer pe(solver);

	//Setup energy window
	const double k_lower_bound = -5.0;
	const double k_upper_bound = 5.0;
	const int k_resolution = 1000;
	pe.setEnergyWindow(k_lower_bound, k_upper_bound, k_resolution);

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
        printf("// Author: Robin Worreby                                             //\n");
        printf("// License: Use if you like, but give me credit.                      //\n");
        printf("/**********************************************************************/\n");
        printf("\n");
    }
    printf("Usage: \n./Application molecule.xyz [parameters]\nParameters:\n");
    printf("  -p or --periodic \t\t- sets the periodicity direction and distance, two parameters needed [X, Y, Z] and [distance]\n");
    printf("  -t or --hopping_amplitude \t- sets the Hamiltonian parameter t value (Hopping amplitude), default is 1.0\n");
    printf("  -b or --bond_threshold \t- sets the bond threshold in Ångström\n");
    printf("  -H or --Hubbard \t\t- sets the Hamiltonian parameter U, default is 0.0 \n");
    printf("  -h or --help \t\t\t- prints this help info.\n");
}
