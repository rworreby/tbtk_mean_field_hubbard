// c++ xyz_periodic_handler.cpp -o main && ./main gnr_7_periodic -p X 4.26 -t 1.8
// ./build/Application gnr_7_periodic.xyz -p X 4.26 -t 1.8

#include <algorithm> //std::sort()
#include <complex>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip> //std::setw()
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>

#include "TBTK/Array.h"
#include "TBTK/BrillouinZone.h"
#include "TBTK/Model.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Vector3d.h"

using namespace TBTK;
using comperator_t = std::function<bool(double, double)>;
using bonds_t = std::vector<std::pair<int, int>>;
using std::string;
using std::size_t;
using std::vector;
using std::complex;
using std::ofstream;

const double k_eps { 0.0001 };
const double atomic_radii_C { 0.68 };
const double threshold { 1e-4 };
std::complex<double> i(0, 1);
Model model;

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
 *  and an unique id to distingush and access the atoms
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
    void construct_molecule_from_file(std::istream& in){
        std::string n_atom;
        getline(in, n_atom);
        std::cout << "Number of atoms in file: " << n_atom << '\n';

        //Read comment line here and parse:
        std::string comment_line;
        getline(in, comment_line);

        //in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        Atom atom;
        while (in >> atom) {
            if(atom.name_coords_.name != "C") continue;
            atoms_.push_back(atom);
            val_x.push_back(atom.name_coords_.xyz.x);
            val_y.push_back(atom.name_coords_.xyz.y);
        }

        set_min_max_x();
    }

    void construct_molecule_from_indexlist(Molecule other_mol,
            std::vector<int> indices, double periodicity_distance){
        for(auto index : indices){
            atoms_.push_back(other_mol.atoms_[index]);
            val_x.push_back(other_mol.atoms_[index].name_coords_.xyz.x);
            val_y.push_back(other_mol.atoms_[index].name_coords_.xyz.y);

            atoms_.push_back(other_mol.atoms_[index]);
            atoms_[atoms_.size()-1].name_coords_.xyz.x += periodicity_distance;
            val_x.push_back(other_mol.atoms_[index].name_coords_.xyz.x);
            val_y.push_back(other_mol.atoms_[index].name_coords_.xyz.y);
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

    friend std::ostream& operator<<(std::ostream& out, const Molecule& atoms){
        for(auto i : atoms.atoms_){
            out << i;
        }
        return out;
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

                //if two atoms are inside the sum of their atomic radii
                //(plus a threshold that is provided by the user, default 1.3)
                //there is a bond between them
                if(result < bond_threshold_){
                    if(i % 2 == 1 && j % 2 == 1){
                        continue;
                    }
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


bool bond_is_in_unit_cell(bonds_t bonds_in_unit_cell,
        int first_atom, int second_atom){
    for(auto bond : bonds_in_unit_cell){
        if(bond.first == first_atom && bond.second == second_atom){
            return true;
        }
    }
    return false;
}


enum h_state_enum {h_unit_cell_enum,
    h_border_crossing_enum,
    h_both_enum,
    undefined_enum
};

std::ostream& operator<<(std::ostream& os, h_state_enum c)
{
    switch(c)
    {
        case h_both_enum: os << "Both Unit Cell & Border Crossing";    break;
        case h_unit_cell_enum: os << "Unit Cell"; break;
        case h_border_crossing_enum : os << "Border Crossing";  break;
        case undefined_enum  : os << "Undefined";   break;
        default    : os.setstate(std::ios_base::failbit);
    }
    return os;
}

int main(int argc, char **argv) {
    string periodicity_direction { "" };
    double periodicity_distance { 0.0 };
    complex<double> t { 1.0 };
    double threshold { 1.7 };
    double hubbard { 0.0 };

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
                hubbard = std::stod(argv[++i]);
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

    std::cout << "Parameter values: " << '\n';
    PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
    PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(hubbard);

    //Set the natural units for this calculation.
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

    size_t BRILLOUIN_ZONE_RESOLUTION = 1000;
	const int K_POINTS_PER_PATH = BRILLOUIN_ZONE_RESOLUTION/10;
	vector<unsigned int> numMeshPoints = {
		BRILLOUIN_ZONE_RESOLUTION,
        1,
        1
	};
	const size_t ENERGY_RESOLUTION = 1000;
	const double ENERGY_LOWER_BOUND = -10;
	const double ENERGY_UPPER_BOUND = 10;

    Molecule whole_molecule;
    // Add molecules from xyz file
    whole_molecule.construct_molecule_from_file(in);

    std::vector<int> atoms_in_unit_cell;
    double x_min = whole_molecule.get_x_min();
    double x_max = whole_molecule.get_x_min() + periodicity_distance;

    for(size_t i = 0; i < whole_molecule.size(); ++i){
        double mol_x_coord = whole_molecule.get_x_coords(i);
        if(mol_x_coord + k_eps < x_max && mol_x_coord >= x_min - k_eps){
            atoms_in_unit_cell.push_back(i);
        }
    }

    std::cout << "Selected atom in unit cell: " << '\n';
    for(size_t j = 0; j < atoms_in_unit_cell.size(); ++j){
        std::cout << atoms_in_unit_cell[j] << " ";
    }
    std::cout << '\n';

    std::cout << "The box spans in x-dir from " << whole_molecule.get_x_min()
              << " to " << (whole_molecule.get_x_min() + periodicity_distance) << '\n';

    Molecule molecule;
    molecule.construct_molecule_from_indexlist(whole_molecule,
            atoms_in_unit_cell, periodicity_distance);

    std::cout << "New (tiny) molecule: " << '\n';
    std::cout << molecule << std::endl;

    Bonds bonds;
    bonds.add_bonds(molecule, threshold, t);
    //bonds.print_bonds();

    bonds_t bonds_in_unit_cell = bonds.bonds_;
    std::cout << "Bonds in molecule: " << '\n';
    for(auto el : bonds_in_unit_cell){
        std::cout << "(" << el.first << "," << el.second << "), ";
    }
    std::cout << '\n';
    std::cout << std::boolalpha;

    double unit_cell_size_y = molecule.find_unit_cell_size();

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

    PRINTVAR(r_AB[0]); PRINTVAR(r_AB[1]); PRINTVAR(r_AB[2]);

    r[1] = Vector3d({0,	periodicity_distance,	0});
    r[1] = Vector3d({0,	periodicity_distance,	0});
    r[2] = Vector3d({0,	0, periodicity_distance});


    // TODO: Compare r_AB to distance vectors computed below.

    Vector3d k_orig[3];
	for(unsigned int n = 0; n < 3; n++){
		k_orig[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	}

    std::cout << "Not crashed yet -1" << '\n';

    std::cout << "y unit cell size: " << unit_cell_size_y << '\n';
    std::cout << "What is k[0].x: " << k_orig[0].x << '\n';
    std::cout << "What is k[0].y: " << k_orig[0].y << '\n';
    std::cout << "What is k[0].z: " << k_orig[0].z << '\n';

	//Setup the BrillouinZone.
    BrillouinZone brillouinZone(
		{
			{k_orig[0].x,k_orig[0].y,k_orig[0].z},
			{k_orig[1].x,k_orig[1].y,k_orig[1].z},
			{k_orig[2].x,k_orig[2].y,k_orig[2].z},
		},
		SpacePartition::MeshType::Nodal
	);

    std::cout << "Not crashed yet 0" << '\n';
    //Create mesh.
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(
		numMeshPoints
	);

    std::cout << "Not crashed yet 1" << '\n';
	//Setup model.
	//model = create_hamiltonian(molecule, bonds, t, hubbard);

    std::cout << "mesh size: " << mesh.size() << '\n' << '\n' << '\n';
    for(unsigned int m = 0; m < mesh.size(); m++){
        Index kIndex = brillouinZone.getMinorCellIndex(
            mesh[m],
            numMeshPoints
        );

        Vector3d kmesh({mesh[m][0], mesh[m][1], mesh[m][2]});

        complex<double> one(1, 0);
        // complex<double> h_unit_cell = -t;
        // complex<double> h_border_crossing = -t * exp(-i*Vector3d::dotProduct(kmesh, k_orig[0]/3.0));
        // complex<double> h_both = -t * (one + exp(-i*Vector3d::dotProduct(kmesh, k_orig[0])));

        h_state_enum h_state = undefined_enum;

        static int printer = 0;
        const int comp_val = 3;
        for(auto bond : bonds_in_unit_cell){
            bool first_is_odd = true;
            complex<double> h;

            double x_diff {
                molecule.get_x_coords(bond.second)
                - molecule.get_x_coords(bond.first)
            };
            double y_diff {
                molecule.get_y_coords(bond.second)
                - molecule.get_y_coords(bond.first)
            };

            if(printer == comp_val){
                std::cout << "Distance Vector between atoms " << bond.first
                          << ", " << bond.second << " in original bond: ("
                          << x_diff << ", " << y_diff << ", 0)" << '\n';
            }
            Vector3d dist_vec({x_diff, y_diff, 0});

            h = -t * exp(-i*Vector3d::dotProduct(kmesh, dist_vec));
            if(printer == comp_val)
                std::cout << "h for original bond: " << h << '\n';

            if((bond.first + bond.second) % 2 == 0){
                if(bond_is_in_unit_cell(bonds_in_unit_cell,
                        bond.first+1, bond.second)
                        || bond_is_in_unit_cell(bonds_in_unit_cell,
                        bond.first, bond.second+1)){

                    if(printer == comp_val){
                        std::cout << "Skipping hopping for bond (" << bond.first
                              << ", " << bond.second << ")." << '\n';
                    }
                    continue;
                }
                h_state = h_unit_cell_enum;
            }
            else{
                if(bond.second % 2 == 1){
                    first_is_odd = false;
                }
                if(bond_is_in_unit_cell(bonds_in_unit_cell,
                        bond.first-first_is_odd, bond.second-(1-first_is_odd))){
                    double second_x_diff { 0.0 };
                    double second_y_diff { 0.0 };
                    if(first_is_odd){
                        second_x_diff =
                            molecule.get_x_coords(bond.first-first_is_odd)
                            - molecule.get_x_coords(bond.second-(1-first_is_odd));
                        second_y_diff =
                            molecule.get_y_coords(bond.first-first_is_odd)
                            - molecule.get_y_coords(bond.second-(1-first_is_odd));
                    }
                    else{
                        second_x_diff =
                            molecule.get_x_coords(bond.second-(1-first_is_odd))
                            - molecule.get_x_coords(bond.first-first_is_odd);

                        second_y_diff =
                            molecule.get_y_coords(bond.second-(1-first_is_odd))
                            - molecule.get_y_coords(bond.first-first_is_odd);
                    }

                    Vector3d second_dist_vec({second_x_diff, second_y_diff, 0});

                    h = -t * (exp(-i*Vector3d::dotProduct(kmesh, dist_vec)) + exp(-i*Vector3d::dotProduct(kmesh, second_dist_vec)));
                    h_state = h_both_enum;
                    if(printer > comp_val && printer < comp_val + 3){
                        std::cout << "h both is: " << h << '\n';
                        std::cout << "dist vec: " << dist_vec << '\n';
                        std::cout << "dist vec2: " << second_dist_vec << '\n';
                    }
                }
                else{
                    if(first_is_odd){
                        if(printer == comp_val)
                            std::cout << "Spec Border hopping distance vector between atoms ("
                            << (bond.first - 1) << ", " << bond.second << "): "
                            << dist_vec << '\n';
                        h = -t * exp(-i*Vector3d::dotProduct(kmesh, dist_vec));
                        if(printer == comp_val){
                            std::cout << "h for spec border hopping:" << h << '\n';
                            std::cout << "Comparing dist vec with minus dist vec:"
                                      << -t * exp(-i*Vector3d::dotProduct(kmesh, dist_vec))
                                      << " against "
                                      << -t * exp(-i*Vector3d::dotProduct(kmesh, -dist_vec))
                                      << '\n';
                        }
                    }
                    h_state = h_border_crossing_enum;
                }
            }
            int first_atom { bond.first };
            int second_atom { bond.second };
            if(h_state == h_border_crossing_enum){
                double x_diff_test = molecule.get_x_coords(first_atom - first_is_odd) - molecule.get_x_coords(second_atom - (1-first_is_odd));
                double y_diff_test = molecule.get_y_coords(first_atom - first_is_odd) - molecule.get_y_coords(second_atom - (1-first_is_odd));
                Vector3d dist_vec_test({x_diff_test, y_diff_test, 0});
                h = -t * (one + exp(-i*Vector3d::dotProduct(kmesh, dist_vec_test)));
                model << HoppingAmplitude(
        			h,
        			{kIndex[0], kIndex[1], kIndex[2], (first_atom - first_is_odd)},
                    {kIndex[0], kIndex[1], kIndex[2], (second_atom - (1-first_is_odd))}
        		) + HC;
                if(printer == comp_val){
                    std::cout << "Added hoppings from " << (second_atom - (1-first_is_odd))
                            << " to " << (first_atom - first_is_odd)
                            << " as " << h_state << " with h value " << h << '\n';
                    if(h_state = h_border_crossing_enum)
                        std::cout << "CALLED BORDER CROSSING FROM 0" << '\n';
                }
                continue;
            }
            if(first_atom % 2 == 1){
                first_atom -= 1;

                model << HoppingAmplitude(
        			h,
        			{kIndex[0], kIndex[1], kIndex[2], second_atom},
        			{kIndex[0], kIndex[1], kIndex[2], first_atom}
        		) + HC;

                if(printer == comp_val){
                    std::cout << "Added hoppings from " << second_atom
                            << " to " << first_atom
                            << " as " << h_state << " with h value " << h << '\n';
                    if(h_state = h_border_crossing_enum)
                        std::cout << "CALLED BORDER CROSSING FROM 1" << '\n';
                }
            }
            else if(second_atom % 2 == 1){
                second_atom -= 1;
                model << HoppingAmplitude(
        			h,
        			{kIndex[0], kIndex[1], kIndex[2], first_atom},
        			{kIndex[0], kIndex[1], kIndex[2], second_atom}
        		) + HC;
                if(printer == comp_val){
                    std::cout << "Added hoppings from " << first_atom
                            << " to " << second_atom
                            << " as " << h_state << " with h value " << h << '\n';
                    if(h_state = h_border_crossing_enum)
                        std::cout << "CALLED BORDER CROSSING FROM 2" << '\n';
                }
            }
            else{
                model << HoppingAmplitude(
        			h,
        			{kIndex[0], kIndex[1], kIndex[2], first_atom},
        			{kIndex[0], kIndex[1], kIndex[2], second_atom}
        		) + HC;
                if(printer == comp_val){
                    std::cout << "Added hoppings from " << first_atom
                            << " to " << second_atom
                            << " as " << h_state << " with h value " << h << '\n';
                    if(h_state = h_border_crossing_enum)
                        std::cout << "CALLED BORDER CROSSING FROM 3" << '\n';
                }
            }

            h_state = undefined_enum;
        }
        printer++;
    }
    std::cout << '\n';
    std::cout << "Not crashed yet 2" << '\n';

    model.construct();

    //Setup the solver.
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
    solver.setVerbose(true);
	solver.run();

	//Setup the property extractor.
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	propertyExtractor.setEnergyWindow(
		ENERGY_LOWER_BOUND,
		ENERGY_UPPER_BOUND,
		ENERGY_RESOLUTION
	);

	//Calculate the density of states.
	Property::DOS dos = propertyExtractor.calculateDOS();


    std::cout << "Not crashed yet 3" << '\n';

    Vector3d lowest({-M_PI / periodicity_distance, 0, 0});
    Vector3d highest({M_PI / periodicity_distance, 0, 0});

    std::cout << "Not crashed yet 3.1" << '\n';

    vector<vector<Vector3d>> paths = {
        {lowest, highest}
    };

    std::cout << "Not crashed yet 3.2" << '\n';

    //Calculate the band structure along the path Gamma -> M -> K -> Gamma.
	//Array<double> bandStructure({K_POINTS_PER_PATH}, 0);
	Range interpolator(0, 1, K_POINTS_PER_PATH);

    std::cout << "Not crashed yet 3.3" << '\n';

    ofstream myfile;
	myfile.open("results.txt");
	std::cout << "Writing results to file" << '\n';
    for (size_t i = 0; i < atoms_in_unit_cell.size(); i++) {
        myfile << "value" << i;
        if(i != atoms_in_unit_cell.size()-1){
            myfile << ", ";
        }
        else{
            myfile << std::endl;
        }
    }
    //myfile << "value0, value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, value11, value12, value13" << std::endl;
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

            //std::cout << "Not crashed yet 5" << '\n';
			//bandStructure[{0, n}] =
            for(size_t i = 0; i < atoms_in_unit_cell.size(); ++i){
                myfile << propertyExtractor.getEigenValue(kIndex, atoms_in_unit_cell[i]);
                if(i != atoms_in_unit_cell.size()-1)
                    myfile << ", ";
                else{
                    myfile << std::endl;
                }
            }
			// //bandStructure[{1, n}] = propertyExtractor.getEigenValue(kIndex, 1);
		}
	}

    std::cout << "Not crashed yet 6" << '\n';

    int basisSize = model.getBasisSize();


	std::cout << "Writing done." << '\n';
	myfile.close();


    /**
	//Find max and min value for the band structure.
	double min = bandStructure[{0, 0}];
	double max = bandStructure[{1, 0}];
	for(unsigned int n = 0; n < K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0, n}])
			min = bandStructure[{0, n}];
		if(max < bandStructure[{1, n}])
			max = bandStructure[{1, n}];
	}
    **/
    return 0;
}

void print_help(bool full){
    if(full){
        printf("/**********************************************************************/\n");
        printf("// Tight-Binding and Mean-Field Hubbard approximation                 //\n");
        printf("// Bachelor's thesis, Spring Semester 2019, ETH Zurich                //\n");
        printf("// Authors: Robin Worreby                                             //\n");
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
