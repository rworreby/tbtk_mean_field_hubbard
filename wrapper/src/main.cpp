// c++ xyz_periodic_handler.cpp -o main && ./main gnr_7_periodic -p X 4.26 -t 1.8
// ./build/Application gnr_7_periodic.xyz -p X 4.26 -t 1.8

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
#include "TBTK/BrillouinZone.h"
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
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Vector3d.h"

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
double U { 0.0 };
double spin_and_site_resolved_density_tol { 1e-5 };
double DENSITY_TOLLERANCE { 1e-7 };
double TARGET_DENSITY_PER_SITE { 1.0 };
double MIXING_PARAMETER { 0.5 };
Model model;
Array<double> spinAndSiteResolvedDensity;
size_t k_num_atoms { 0 };
size_t k_num_atoms_unit_cell { 0 };
const double k_eps { 0.0001 };
const double atomic_radii_C { 0.68 };
const double threshold { 1e-4 };
bool periodic_hubbard { false };

//USAGE:     DEBUG(debug_var++);
size_t debug_var = 0;
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
		set_min_max_x();
    }

	double find_unit_cell_size(){
        double max_val=0, min_val=0;
        std::vector<double> temp_y = val_y;
        comperator_t compFunctor = [](double elem1, double elem2){
				return elem1 < elem2;
		};
        std::sort(temp_y.begin(), temp_y.end(), compFunctor);

        min_y = min_val = temp_y[0];
        max_y = max_val = temp_y[temp_y.size()-1];

        return max_val - min_val;
    }

    void set_min_max_x(){
        double max_val=0, min_val=0;
        std::vector<double> temp_x = val_x;
        comperator_t compFunctor = [](double elem1, double elem2){
				return elem1 < elem2;
		};
        std::sort(temp_x.begin(), temp_x.end(), compFunctor);

        min_x = min_val = temp_x[0];
        max_x = max_val = temp_x[temp_x.size()-1];
    }

    size_t size(){
        return atoms_.size();
    }
	double get_x_coords(int index){
        return val_x.at(index);
    }
    double get_y_coords(int index){
        return val_y.at(index);
    }
    double get_x_max(){
        return max_x;
    }
    double get_x_min(){
        return min_x;
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
	double min_x;
	double max_x;
	double min_y;
	double max_y;
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

    void add_bonds(Molecule mol, double threshold, complex<double> t){
        bond_threshold_ = threshold;
		t_ = t;
        std::cout << "Size of the considered molecule (no H): " << mol.size() << "\n";
        for(int i = 0; i < mol.size(); ++i){
            double x_1 = mol.atoms_[i].name_coords_.xyz.x;
            double y_1 = mol.atoms_[i].name_coords_.xyz.y;
            double z_1 = mol.atoms_[i].name_coords_.xyz.z;

            for (int j = i; j < mol.size(); j++){
                if(i == j) continue;
                double x_2 = mol.atoms_[j].name_coords_.xyz.x;
                double y_2 = mol.atoms_[j].name_coords_.xyz.y;
                double z_2 = mol.atoms_[j].name_coords_.xyz.z;
                double x_diff = x_1 > x_2 ? x_1 - x_2 : x_2 - x_1;
                double y_diff = y_1 > y_2 ? y_1 - y_2 : y_2 - y_1;
                double z_diff = z_1 > z_2 ? z_1 - z_2 : z_2 - z_1;
                double result = std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

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
	complex<double> t_;
    friend class Molecule;

};

bool check_straight_bond_orientation(double x_1, double x_2,
                                     double y_1, double y_2){
    bool same_x = (x_1 > x_2 - k_eps) && (x_1 < x_2 + k_eps);
    bool same_y = (y_1 > y_2 - k_eps) && (y_1 < y_2 + k_eps);
    return same_x || same_y;
}

Model create_hamiltonian_molecule(Molecule mol, Bonds bonds, complex<double> t_, bool hubbard){
    std::cout << "Hamiltonian Molecule: Not crashed yet 1" << '\n';

	for(int s = 0; s < 2; ++s){
        for(int i = 0; i < bonds.bonds_.size(); ++i){
            model << HoppingAmplitude(
				-t_,
				{0, s, bonds.bonds_[i].first},
				{0, s, bonds.bonds_[i].second}
			) + HC;
        }
    }
    std::cout << "Hamiltonian Molecule: Not crashed yet 2" << '\n';

    if(hubbard){
        for(int s = 0; s < 2; ++s){
            for(int i = 0; i < mol.size(); ++i){
                model << HoppingAmplitude(
					H_U,
					{0, s, i},
					{0, s, i}
				);
            }
        }
    }
    std::cout << "Hamiltonian Molecule: Not crashed yet 3" << '\n';

    model.construct();
	model.setTemperature(0.01);

    std::cout << "Model size in create_hamiltonian: " << model.getBasisSize() << std::endl;
    return model;
}

Model create_hamiltonian_periodic(Molecule molecule,
		std::vector<int> atoms_in_unit_cell,
        std::vector<int> atoms_on_unit_cell_border, Bonds bonds,
        bonds_t bonds_in_unit_cell, BrillouinZone brillouinZone,
		vector<vector<double>> mesh, vector<unsigned int> numMeshPoints,
		Vector3d r[3], complex<double> t, bool hubbard){
	for(int m = 0; m < mesh.size(); m++){
        Index kIndex = brillouinZone.getMinorCellIndex(
            mesh[m],
            numMeshPoints
        );

        Vector3d k({mesh[m][0], mesh[m][1], mesh[m][2]});

        complex<double> h_within_unit_cell = -t;
        complex<double> h_next_unit_cell = -t * exp(-i*Vector3d::dotProduct(k, r[0]));

        for(auto unit_cell_bond : bonds_in_unit_cell){
			for(int s = 0; s < 2; s++){
	            model << HoppingAmplitude(
	    			h_within_unit_cell,
	    			{kIndex[0], kIndex[1], kIndex[2], unit_cell_bond.first, s},
	    			{kIndex[0], kIndex[1], kIndex[2], unit_cell_bond.second, s}
	    		) + HC;
			}
        }
        // TODO: implement!
        for(auto border_atom : atoms_on_unit_cell_border){

        }
        if(hubbard){
    		for(auto atom : atoms_in_unit_cell){
    			for (int s = 0; s < 2; s++) {
					model << HoppingAmplitude(
						H_U,
						{kIndex[0], kIndex[1], kIndex[2], atom, s},
						{kIndex[0], kIndex[1], kIndex[2], atom, s}
					);
				}
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
	for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms; site++){
			spinAndSiteResolvedDensity[
				{spin, site}
			] = (rand()%100)/100.0; //1.0;
		}
	}
}

//Callback function responsible for returning the current value of H_U for the
//given indices. Since H_U is supposed to be diagonal, toIndex and fromIndex is
//assumed to be equal and only the fromIndex is used to determine the spin and
//site. Compare to the model specification where H_C only is passed to the
//model with both the 'to' and 'from' indices equal.
complex<double> H_U(const Index &toIndex, const Index &fromIndex){
	int spin = fromIndex[1];
	int site = fromIndex[0];
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
void fixDensity(PropertyExtractor::BlockDiagonalizer &propertyExtractor){
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
		for(size_t n = 0; n < model.getBasisSize(); n++){
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
bool selfConsistencyCallback(Solver::BlockDiagonalizer &solver){
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

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
	    {0, IDX_ALL, IDX_ALL}
    });

	//Update the spin and site resolved density. Mix with the previous
	//value to stabilize the self-consistent calculation.
	for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms; site++){
			spinAndSiteResolvedDensity[{spin, site}]
				= MIXING_PARAMETER*spinAndSiteResolvedDensity[
					{spin, site}
				] + (1 - MIXING_PARAMETER)*density({
                    0,
                    (int)spin,
                    (int)site,
				});///(model.getBasisSize()/k_num_atoms);
		}
	}

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
    int max_spin {0}, max_site{0};
	double maxDifference = 0;
	for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms; site++){
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

//Callback function that is to be called each time the model Hamiltonian has
//been diagonalized. The function first fixes the density to the target
//density, then calculates the spin and site resolved density. The new spin and
//site resolved density is mixed with the previous values to stabilize the
//self-consistent calculation.
bool selfConsistencyCallbackPeriodic(Solver::BlockDiagonalizer &solver){
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

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
	    {IDX_SUM_ALL, IDX_SUM_ALL, IDX_SUM_ALL, IDX_ALL, IDX_ALL}
    });

	//Update the spin and site resolved density. Mix with the previous
	//value to stabilize the self-consistent calculation.
    for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms_unit_cell; site++){
			spinAndSiteResolvedDensity[{spin, site}]
				= MIXING_PARAMETER*spinAndSiteResolvedDensity[
					{spin, site}
				] + (1 - MIXING_PARAMETER)*density({
                    IDX_SUM_ALL,
                    IDX_SUM_ALL,
                    IDX_SUM_ALL,
					(int)site,
                    (int)spin
				}); ///(model.getBasisSize()/k_num_atoms);
		}
	}

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
    int max_spin {0}, max_site{0};
	double maxDifference = 0;
	for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms; site++){
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
    double periodicity_distance { 1.0 };
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
                periodic_hubbard = true;
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
	size_t BRILLOUIN_ZONE_RESOLUTION = 1000;
	const size_t K_POINTS_PER_PATH = BRILLOUIN_ZONE_RESOLUTION/10;
	vector<unsigned int> numMeshPoints = {
		BRILLOUIN_ZONE_RESOLUTION,
		1,
		1
	};

	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	std::cout << "Parameter values: " << '\n' << std::boolalpha;
    PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
    PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(hubbard); PRINTVAR(periodic);

    //Usage example
    Molecule molecule;

    //Add molecules from xyz file
    molecule.construct_molecule(in);
    //std::cout << molecule;

    k_num_atoms = molecule.size();

    Bonds bonds;

    //Create bonds in molecule
    bonds.add_bonds(molecule, threshold, t);


	std::vector<int> atoms_in_unit_cell;
    std::vector<int> atoms_on_unit_cell_border;
	double x_min = molecule.get_x_min();
	double x_max = molecule.get_x_min() + periodicity_distance;

	if(periodic){
		for(size_t i = 0; i < molecule.size(); ++i){
			double mol_x_coord = molecule.get_x_coords(i);
            if(mol_x_coord + k_eps < x_max && mol_x_coord >= x_min - k_eps){
                atoms_in_unit_cell.push_back(i);
                if(mol_x_coord > x_max - 0.8){
                    atoms_on_unit_cell_border.push_back(i);
                }
            }
		}

		std::cout << "Selected atom in unit cell: " << '\n';
	    for(size_t j = 0; j < atoms_in_unit_cell.size(); ++j){
	        std::cout << atoms_in_unit_cell[j] << " ";
	    }
	    std::cout << '\n';
        std::cout << "Selected atom at unit cell border: " << '\n';
        for(size_t j = 0; j < atoms_on_unit_cell_border.size(); ++j){
            std::cout << atoms_on_unit_cell_border[j] << " ";
        }
        std::cout << '\n';

	    std::cout << "The box spans in x-dir from " << molecule.get_x_min()
	              << " to " << (molecule.get_x_min() + periodicity_distance) << '\n';
        k_num_atoms_unit_cell = atoms_in_unit_cell.size();
	}

    bonds_t bonds_in_unit_cell;
    bonds_t bonds_crossing_unit_cell_borders;

	if(periodic){
	    std::cout << "First elements of the bonds: " << '\n';
	    for(auto el : bonds.bonds_){
	        bool first_element_in_unit_cell = false;
	        bool second_element_in_unit_cell = false;
	        for(size_t i = 0; i < atoms_in_unit_cell.size(); ++i){
	            if(el.first == atoms_in_unit_cell[i]){
	                first_element_in_unit_cell = true;
	            }
	            if(el.second == atoms_in_unit_cell[i]){
	                second_element_in_unit_cell = true;
	            }
	            if(first_element_in_unit_cell && second_element_in_unit_cell){
	                bonds_in_unit_cell.push_back(el);
	                break;
	            }
	        }
		}

        // TODO: Correct this function so only bonds ranging over the edges of
        //       a unit cell are considered.
	    std::cout << '\n';
        for(auto el : bonds.bonds_){
	        bool one_element_in_unit_cell = false;
	        bool one_element_in_next_unit_cell = false;
	        for(size_t i = 0; i < molecule.size(); ++i){
	            if(el.first == atoms_on_unit_cell_border[i]){
	                one_element_in_unit_cell = true;
	            }
	            if(el.second == atoms_on_unit_cell_border[i]){
	                one_element_in_next_unit_cell = true;
	            }
	            if(one_element_in_unit_cell && one_element_in_next_unit_cell){
	                break;
	            }
                // TODO: add check for right hand side border
                if(one_element_in_unit_cell || one_element_in_next_unit_cell){
                    //for(size_t i = 0; i < atoms_on_unit_cell_border.size(); ++i){

                    bonds_crossing_unit_cell_borders.push_back(el);
                    //}
                }
	        }
            std::cout << "Val: " << one_element_in_unit_cell << " " <<  one_element_in_next_unit_cell << '\n';
		}
        std::cout << "Found bonds going from unit cell to neighboring cell:" << '\n';
        for(auto el : bonds_crossing_unit_cell_borders){
            std::cout << "(" << el.first << "," << el.second << "), ";
        }
        std::cout << '\n';
	    std::cout << "Chosen pairs in unit cell: " << '\n';
	    for(auto el : bonds_in_unit_cell){
	        std::cout << "(" << el.first << "," << el.second << "), ";
		}
	    std::cout << '\n';
	    std::cout << std::boolalpha;
	    std::cout << "Checking if two atoms are straight:" << '\n';

	    std::cout << "Atom 1 and 2: "
	            << check_straight_bond_orientation(molecule.get_x_coords(1),
	            molecule.get_x_coords(2), molecule.get_y_coords(1),
	            molecule.get_y_coords(2)) << '\n';

	    std::cout << "Atom 2 and 3: "
	            << check_straight_bond_orientation(molecule.get_x_coords(2),
	            molecule.get_x_coords(3), molecule.get_y_coords(2),
	            molecule.get_y_coords(3)) << '\n';
	}
	//std::cout << molecule << '\n';
    std::cout << "X coords of atom 2: " << molecule.get_x_coords(2) << '\n';
    std::cout << "Y coords of atom 2: " << molecule.get_y_coords(2) << '\n';

    std::cout << "Not crashed yet -3.5" << '\n';
    double unit_cell_size_y = molecule.find_unit_cell_size();
    std::cout << "Unit cell size x: " << periodicity_distance << '\n';
    std::cout << "Unit cell size y: " << unit_cell_size_y << '\n';

    std::cout << "X coords of atom 2: " << molecule.get_x_coords(2) << '\n';
    std::cout << "Y coords of atom 2: " << molecule.get_y_coords(2) << '\n';

    std::cout << "Not crashed yet -3" << '\n';

	// Setup lattice vector. By using three three-dimensional vectors
	// instead of two two-dimensional vectors, the cross product expression
	// for the reciprocal vectors can be expressed in terms of cross
	// products.
    Vector3d r[3];
    r[0] = Vector3d({periodicity_distance, 0, 0});
    r[1] = Vector3d({-periodicity_distance/2, periodicity_distance*sqrt(3)/2, 0});
    r[2] = Vector3d({0, 0, periodicity_distance});

    std::cout << "Not crashed yet -2" << '\n';

    // Calculate the reciprocal lattice vectors.
    Vector3d r_AB[3];
    r_AB[0] = (r[0] + 2*r[1])/3.;
	r_AB[1] = -r[1] + r_AB[0];
	r_AB[2] = -r[0] - r[1] + r_AB[0];
	r[1] = Vector3d({0,	periodicity_distance,	0});

    Vector3d k[3];
	for(int n = 0; n < 3; n++){
		k[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	}
	std::cout << "What is k[0].x: " << k[0].x << '\n';
    std::cout << "What is k[0].y: " << k[0].y << '\n';
    std::cout << "What is k[0].z: " << k[0].z << '\n';

	//Setup the BrillouinZone.
    BrillouinZone brillouinZone(
		{
			{k[0].x,k[0].y,k[0].z},
			{k[1].x,k[1].y,k[1].z},
			{k[2].x,k[2].y,k[2].z},
		},
		SpacePartition::MeshType::Nodal
	);

    std::cout << "Not crashed yet 0" << '\n';
    //Create mesh.
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(
		numMeshPoints
	);

    std::cout << "Not crashed yet 1" << '\n';


    Array<double> dummy_array({2, molecule.size()});
    spinAndSiteResolvedDensity = dummy_array;

    initSpinAndSiteResolvedDensity(spinAndSiteResolvedDensity);

    std::cout << "Not crashed yet 2" << '\n';
	if(!periodic){
    	model = create_hamiltonian_molecule(molecule, bonds, t, hubbard);
	}
	else{
        model = create_hamiltonian_periodic(molecule, atoms_in_unit_cell,
            atoms_on_unit_cell_border, bonds, bonds_in_unit_cell,
            brillouinZone, mesh, numMeshPoints, r, t, hubbard
            );
	}

    std::cout << "Not crashed yet 3" << '\n';

	Solver::BlockDiagonalizer solver;

	solver.setModel(model);

    std::cout << "Geometry of model: " << solver.getModel().getGeometry() << '\n';
	if(hubbard){
        if(periodic){
            solver.setSelfConsistencyCallback(selfConsistencyCallbackPeriodic);
    	    solver.setMaxIterations(1000);
        }
        else{
            solver.setSelfConsistencyCallback(selfConsistencyCallback);
    	    solver.setMaxIterations(1000);
        }
    }
    std::cout << "Not crashed yet 4" << '\n';

	//Run the solver. This will run a self-consistent loop where the
	//Hamiltonian first is diagonalized, and then the
	//selfConsistencyCallback is called. The procedure is repeated until
	//either the self-consistency callback returns true, or the maximum
	//number of iterations is reached.
    solver.run();

    std::cout << "Not crashed yet 5" << '\n';

    if(hubbard) std::cout << "Final chemical potential: " << model.getChemicalPotential() << '\n';
	//Create PropertyExtractor
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

	//Setup energy window
	const double k_lower_bound = -5.0;
	const double k_upper_bound = 5.0;
	const size_t k_resolution = 1000;
	propertyExtractor.setEnergyWindow(k_lower_bound, k_upper_bound, k_resolution);

    std::cout << "Not crashed yet 6" << '\n';

	Vector3d lowest({-M_PI / periodicity_distance, 0, 0});
    Vector3d highest({M_PI / periodicity_distance, 0, 0});
    vector<vector<Vector3d>> paths = {
        {lowest, highest}
    };
	Range interpolator(0, 1, K_POINTS_PER_PATH);

//	size_t basisSize = model.getBasisSize();
	ofstream myfile;
	myfile.open("eigenval_eigenvec.txt");
	std::cout << "Writing Eigenvalues and Eigenvectors to file" << '\n';

    //Property::EigenValues eigen_values = propertyExtractor.getEigenValues();
	//FileWriter::writeEigenValues(eigen_values);

    if(periodic){
        for (size_t i = 0; i < atoms_in_unit_cell.size(); i++) {
            myfile << "value" << i;
            if(i != atoms_in_unit_cell.size()-1){
                myfile << ", ";
            }
            else{
                myfile << std::endl;
            }
        }
        for(unsigned int p = 0; p < 1; p++){
			Vector3d startPoint = paths[p][0];
			Vector3d endPoint = paths[p][1];

	        std::cout << "startPoint: " << startPoint << '\n';
	        std::cout << "endPoint: " << endPoint << '\n';
			for(unsigned int n = 0; n < K_POINTS_PER_PATH; n++){
				Vector3d k = (
					interpolator[n]*endPoint
					+ (1 - interpolator[n])*startPoint
				);

				Index kIndex = brillouinZone.getMinorCellIndex(
					{k.x, k.y, k.z},
					numMeshPoints
				);

	            for(size_t i = 0; i < atoms_in_unit_cell.size(); ++i){
	                myfile << propertyExtractor.getEigenValue(kIndex, atoms_in_unit_cell[i]);
	                if(i != atoms_in_unit_cell.size()-1)
	                    myfile << ", ";
	                else{
	                    myfile << std::endl;
	                }
	            }
			}
		}
    }
    else{
        /**
        std::cout << "All Eigenvalues: " << '\n';
        for (size_t k = 0; k < k_num_atoms*2; k++) {
            std::cout << propertyExtractor.getEigenValue(k) << " ";
        }
        **/
        std::cout << '\n';
        std::cout << "Eigenvalues polarisation test: " << '\n';
        std::cout << propertyExtractor.getEigenValue(k_num_atoms/2-2) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms/2-1) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms/2) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms/2+1) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms+k_num_atoms/2-2) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms+k_num_atoms/2-1) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms+k_num_atoms/2) << "\n";
        std::cout << propertyExtractor.getEigenValue(k_num_atoms+k_num_atoms/2+1) << "\n";


            /**
            tring spinchannel = "";
            switch(s){
                case 0:
                    spinchannel = "su";
                    break;
                case 1:
                    spinchannel = "sd";
                    break;
                default:
                    break;
            }
            **/
        for(int spin = 0; spin < 2; spin++){
            string spinchannel = "su";
            if(spin == 1)
                spinchannel = "sd";

            for(int i = 0; i < k_num_atoms; ++i){
                myfile << spinchannel << " ";
                myfile << propertyExtractor.getEigenValue(k_num_atoms*spin + i) << " ";
                for(int j = 0; j < k_num_atoms; ++j){
                    myfile << propertyExtractor.getAmplitude(i + k_num_atoms*spin, {0, spin, j}) << " ";
                }
                myfile << std::endl;
            }
        }
    }

        /**
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
        **/


	myfile.close();
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
