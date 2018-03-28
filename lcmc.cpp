#include <iostream>
#include <string>
#include "time.h"
#include "stdlib.h"
#include <fstream>
#include "random_mars.cpp"
#include "math.h"
using namespace std;

//I SHOULD CREATE A VECTOR CLASS

//test moving two particles close and far and see potential

class Atom {
	public:
		//positions
		double x, y, z;
		//forces
		double fx, fy, fz;
		//normalized orientation vector
		double ox, oy, oz;
		void set_pos(double new_x, double new_y, double new_z) {
			x=new_x;
			y=new_y;
			z=new_z;
		}
		void set_force(double new_x, double new_y, double new_z) {
			fx=new_x;
			fy=new_y;
			fz=new_z;
		}
		void set_orientvec(double new_x, double new_y, double new_z) {
			double mag = sqrt(new_x*new_x + new_y*new_y + new_z*new_z);
			ox = new_x/mag;
			oy = new_y/mag;
			oz = new_z/mag;
		}
		double get_orientvec(int comp) {
			if (comp==0) {
				return ox;
			}
			else if (comp==1) {
				return oy;
			}
			else {
				return oz;
			}
		}
};

//FUCNTION DECLARATIONS
///////////////////////
///////////////////////
inline double pbc(double pos);
inline double min_image(double dist);
inline double displacement(double x1, double x2, double y1, double y2, double z1, double z2);
inline double dot_product(double x1, double x2, double y1, double y2, double z1, double z2);
void initialization();
//double* lj_total_energy();
int mc_move();
void simulate(int cycles);
bool overlap();
bool in_box();
bool net_force();
double lj_accept(int particle, double new_x, double new_y, double new_z);
inline double mi_sign(double d);
void orient_norm();
double gb_total_energy();
inline double sig(Atom a, Atom b);
//////////////////////
//////////////////////

//GAY-BERNE VARIABLES
/////////////////////
////////////////////
//well-depth parameter
double eps_0 = 1.0;
//epsilon parameter ratio
double eps_es = 1.0/5.0;
//double eps_es = 1.0;
//length and breadth parameter ratio
//double length = 5.0;
//double breadth = 1.0;
double length = 3.0;
double breadth = 1.0;
double len_brh = length/breadth;
//parameter
//double sigma0 = sqrt(2)*breadth;
double sigma0 = 1.0;
//exponent parameters for epsilon functions
double u = 2;
double v = 1;
//parameter
double chi = (len_brh*len_brh - 1)/(len_brh*len_brh + 1);
//parameter
double chi_p = (1 - pow(eps_es,1.0/u))/(1 + pow(eps_es,1.0/u));



//LENNARD JONES POTENTIAL VARIABLES
///////////////////////////////////
//////////////////////////////////
//const double EPSILON = (1.654e-21)*(6.022E+23)/1000.0;
const double EPSILON = 1.0;
//units are in Angstroms and kJ/mol
//const double SIGMA = 3.40;
const double SIGMA = 1.0;
//const double KB = 0.0083144621;
const double KB = 1.0;
//have this around 2.5 - 3 sigma
const double R_CUT = 8.0*SIGMA;
//maximum displacement in simulation
const double MAX_MOVE = .11;
//maximum rotation for LC molecules
//const double MAX_ROT = .5;
//number of atoms
const int ATOMS = 125;
//density parameter
//const double D = .0127;
//const double D = 0.3;
//size of box 
//const double BL = cbrt(ATOMS*1.0/D);
//const double BL = 1000;
//const double BL = 10;
//snapshot parameter for measuring thermodynamic quantities
const int SNAP = 50;
//temperature
//const double TEMP = .95;

//MISCELLANEOUS PARAMETERS
/////////////////////////
/////////////////////////
//overlap parameter
const double D0 = 1.3;
//seed for RNG
int SEED=5;
//const int SEED = 5;
//random number generator
//RanMars *ran = new RanMars(SEED);
//number of cycles used to equilibrate system
const int EQ_CYCLES = 500;
//constants for RDF
//const double K = BL*.5; //range of distances included
const int M = 64; //number of bins
//for visualizing liquid crystals
const double SCALE = sqrt(3);
///////////////////////
//////////////////////


//Arrays for holding all of the data
//////////////////////////////////
/////////////////////////////////
//2D array containing the system
Atom box[ATOMS]={ };
//array used to calculate the RDF
double rdf[M] = { };
//tensor used to calculate the order parameter
double Q[3][3];
///////////////////////////////
//////////////////////////////

double D=.15;
double TEMP=.5;
double BL = cbrt(ATOMS*1.0/D);
double K=.5*BL;
RanMars *ran;

int main(int argc, char* argv[]) {
	//initialize system according to paper
	/*
	//for calculating gay-berne potential graph
	Atom mol1, mol2;
	mol1.set_pos(-5, 0, 0);
	mol2.set_pos(5, 0, 0);
	mol1.set_orientvec(0,1,0);
	mol2.set_orientvec(1,0,0);
	box[0]=mol1;
	box[1]=mol2;
	int t;
	ofstream twopot;
	twopot.open("mix2mix.csv");
	//twopot << "#POTENTIAL GRAPH OF TWO SIDE TO SIDE MOLECULES APPROACHING EACH OTHER" << endl;
	//twopot << "#Displacement Energy sig uidotuj" << endl;
	for (t=0; t<100; t++) {
		twopot << box[1].x*2.0 << "," << gb_total_energy() << endl; // " " << sig(box[0], box[1]) << " " << dot_product(box[0].ox, box[1].ox,box[0].oy,box[1].oy,box[0].oz,box[1].oz) << endl;
		box[0].set_pos(-5 + 5.0*double(t)/100.0, 0, 0);
		box[1].set_pos(5 - 5.0*double(t)/100.0, 0, 0);
	}
	twopot.close()
	*/
	/*
	//for running a bash script
	if(sscanf(argv[1], "%lf", &D) != 1){
		exit(1);
	}
	if(sscanf(argv[2], "%lf", &TEMP) != 1){
		exit(1);
	}
	if(sscanf(argv[3], "%d", &SEED) != 1){
		exit(1);
	}
	*/

	BL = cbrt(ATOMS*1.0/D);
	K = BL*0.5;
	ran = new RanMars(SEED);

	initialization();
	
	//cout << "Boxlength: " << BL << endl;
	//double* e = lj_total_energy();
	//cout << "No net force? " << net_force() << endl;
	
	simulate(5000);
	//ensures there are no overlaps and all particles are in the box at the end of run
	//cout << "All particles in box? " << in_box() << endl;
	//cout << "Any overlaps? " << overlap() << endl;
	//makes file with final positions
	ofstream xyz_final;
	xyz_final.open("final.dat");
	xyz_final << "#NOTE: (0,0) is at bottom left corner, with x values increasing right and y values increasing up." << endl;
	xyz_final << "#x y z" << endl;
	int i;
	for (i=0; i<ATOMS; i++) {
		xyz_final << box[i].x << " " << box[i].y << " " << box[i].z << endl;
	}
	xyz_final.close();
	//makes file for RDF plot
	ofstream rdf_out;
	rdf_out.open("rdf_mc.dat");
	rdf_out << "#Radial Distribution Function" << endl;
	rdf_out << "#x y" << endl;
	for (i=0; i<M; i++) {
		rdf_out << (i*K/M)/SIGMA << " " << rdf[i] << endl;
	}
	rdf_out.close();
	//orient_norm();
	
	return 0;
}

//dot product
inline double dot_product(double x1, double x2, double y1, double y2, double z1, double z2) {
	return x1*x2 + y1*y2 + z1*z2;
}

//GAY-BERNE POTENTIAL FUNCTIONS
inline double eps1(Atom a, Atom b) {
	double orient_prod = dot_product(a.ox, b.ox, a.oy, b.oy, a.oz, b.oz);
	return 1.0/sqrt(1.0 - chi*chi*orient_prod*orient_prod);
}

inline double eps2(Atom a, Atom b) {
	//SHOULD I SWITCH A AND B (I DONT THINK IT SHOULD MATTER)
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	double dz = a.z - b.z;
	double mi_dx, mi_dy, mi_dz;
	mi_dx=mi_sign(dx);
	mi_dy=mi_sign(dy);
	mi_dz=mi_sign(dz);
	double mag = sqrt(dot_product(mi_dx, mi_dx, mi_dy, mi_dy, mi_dz, mi_dz));
	double dx_norm, dy_norm, dz_norm;
	dx_norm = mi_dx/mag;
	dy_norm = mi_dy/mag;
	dz_norm = mi_dz/mag;
	double r_ua = dot_product(dx_norm, a.ox, dy_norm, a.oy, dz_norm, a.oz);
	double r_ub = dot_product(dx_norm, b.ox, dy_norm, b.oy, dz_norm, b.oz);
	double chip_ua_ub = chi_p*dot_product(a.ox, b.ox, a.oy, b.oy, a.oz, b.oz);
	return 1.0 - (chi_p/2.0)*((pow(r_ua + r_ub, 2)/(1+chip_ua_ub))+(pow(r_ua - r_ub, 2)/(1-chip_ua_ub)));
}

inline double eps(Atom a, Atom b) {
	return eps_0*pow(eps1(a,b), v)*pow(eps2(a,b),u);
}

inline double sig(Atom a, Atom b) {
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	double dz = a.z - b.z;
	double mi_dx, mi_dy, mi_dz;
	mi_dx=mi_sign(dx);
	mi_dy=mi_sign(dy);
	mi_dz=mi_sign(dz);
	double mag = sqrt(dot_product(mi_dx, mi_dx, mi_dy, mi_dy, mi_dz, mi_dz));
	double dx_norm, dy_norm, dz_norm;
	dx_norm = mi_dx/mag;
	dy_norm = mi_dy/mag;
	dz_norm = mi_dz/mag;
	double r_ua = dot_product(dx_norm, a.ox, dy_norm, a.oy, dz_norm, a.oz);
	double r_ub = dot_product(dx_norm, b.ox, dy_norm, b.oy, dz_norm, b.oz);
	double chi_ua_ub = chi*dot_product(a.ox, b.ox, a.oy, b.oy, a.oz, b.oz);
	return sigma0*(1.0/sqrt(1.0 - (chi/2.0)*((pow(r_ua + r_ub, 2)/(1+chi_ua_ub))+(pow(r_ua - r_ub, 2)/(1-chi_ua_ub)))));
}

//applies PBC to a given 1D position --
inline double pbc(double pos) {
	return pos-BL*floor(pos/BL);
}

//returns the 1D displacement according to the min image convention
//how to apply nonzero BL to this function??
inline double min_image(double dist) {
	return min(min(abs(dist + BL), abs(dist-BL)), abs(dist));
}

//returns displacement according to min image concention
inline double displacement(double x1, double x2, double y1, double y2, double z1, double z2) {
	return sqrt(pow(min_image(x1-x2),2) + pow(min_image(y1-y2),2) + pow(min_image(z1-z2), 2));
}

//creates file init.dat that contains the initialization of the system
//(0,0) is in the bottom left corner of the box, so should look exactly like config in paper
//configures files to be best used with gnuplot
//how should i initialize with the cube lattice
void initialization() {
	int i, j, k;
	int lim = cbrt(double(ATOMS));
	//cout << lim << endl;
	//cout << "lim"<<lim*lim*lim << endl;
	if (lim*lim*lim<ATOMS) {
		lim+=1;
	}
	double dist = (BL-1.0)/(lim*1.0-1);
	ofstream xyz_init, init_sim;
	xyz_init.open("init.dat");
	init_sim.open("init_sim.xyz");
	xyz_init << "#NOTE: (0,0) is at bottom left corner, with x values increasing right and y values increasing up." << endl;
	xyz_init << "#x y" << endl;
	init_sim << 3*ATOMS << "\n\n";
	i=0;
	int x,y,z;
	for (x=0; x<lim; x++) {
		for (y=0; y<lim; y++) {
			for (z=0; z<lim; z++) {
				double rand_ox = (((ran->uniform())*2)-1);
				double rand_oy = (((ran->uniform())*2)-1);
				double rand_oz = (((ran->uniform())*2)-1);
				box[i].set_pos(pbc((x+.5)*BL/lim), pbc((y+.5)*BL/lim), pbc((z+.5)*BL/lim));
				box[i].set_orientvec(1, 1, 1);
				i++;
				//cout << i << endl;
			}
		}
	}
	for (i=0; i<ATOMS; i++) {
	//	box[i].set_pos(dist*(i%lim), dist*((i/lim)%lim), dist*((i/(lim*lim))%lim));
		box[i].set_force(0.0,0.0,0.0);
		xyz_init << box[i].x << " " << box[i].y << " " << box[i].z << endl;
		init_sim << "He " << 10*box[i].x << " " << 10*box[i].y << " " << 10*box[i].z << endl;
		init_sim << "He " << 10*box[i].x + SCALE*box[i].ox << " " << 10*box[i].y + SCALE*box[i].oy << " " << 10*box[i].z + SCALE*box[i].oz << endl;
		init_sim << "He " << 10*box[i].x - SCALE*box[i].ox << " " << 10*box[i].y - SCALE*box[i].oy << " " << 10*box[i].z - SCALE*box[i].oz << endl;
	}
	init_sim.close();
	xyz_init.close();
}


/*
void calc_forces(int i, int j, double dx, double dy, double dz, double dt, double rep, double attr) {
	double frc = (48.0*EPSILON/(dt*dt))*(rep - .5*attr);
	double f_x = frc*dx;
	double f_y = frc*dy;
	double f_z = frc*dz;
	box[i].set_force(box[i].fx+f_x, box[i].fy+f_y, box[i].fz+f_z);
	box[j].set_force(box[j].fx-f_x, box[j].fy-f_y, box[j].fz-f_z);
}
*/

//does minimum image convention paying attention to sign
inline double mi_sign(double d) {
	if (d<-.5*BL) {
		return BL+d;
	}
	else if (d>.5*BL) {
		return d-BL;
	}
	return d;
}

//returns total energy of the system for lennard jones potential and updates forces
/*double* lj_total_energy() {
	double res = 0;
	double virial = 0;
	static double r[2];
	int i,j;
	for (i=0; i<ATOMS; i++) {
		for (j=i+1; j<ATOMS; j++) {
			double dx = box[i].x - box[j].x;
			double dy = box[i].y - box[j].y;
			double dz = box[i].z - box[j].z;
			double mi_dx, mi_dy, mi_dz;
			mi_dx=mi_sign(dx);
			mi_dy=mi_sign(dy);
			mi_dz=mi_sign(dz);
			double d_sqr = mi_dx*mi_dx + mi_dy*mi_dy + mi_dz*mi_dz;
			double sig_sqr = SIGMA*SIGMA;
			double repulsion = pow(sig_sqr/d_sqr, 6);
			double attraction = pow(sig_sqr/d_sqr, 3);
			if (d_sqr<(R_CUT*R_CUT)) {
				res += 4.0*EPSILON*(repulsion-attraction);
				double frc = (48.0*EPSILON/d_sqr)*(repulsion - .5*attraction);
				double f_x, f_y, f_z;
				f_x = frc*mi_dx;
				f_y = frc*mi_dy;
				f_z = frc*mi_dz;
				//calculate virial here
				box[i].set_force(box[i].fx+f_x, box[i].fy+f_y, box[i].fz+f_z);
				box[j].set_force(box[j].fx-f_x, box[j].fy-f_y, box[j].fz-f_z);
				virial+=f_x*dx+f_y*dy+f_z*dz;
			}
		}
	}
	//adding the tail correction due to R_CUT
	res += (8.0/3.0)*M_PI*D*ATOMS*EPSILON*pow(SIGMA,3)*((1.0/3.0)*pow(SIGMA/R_CUT, 9)-pow(SIGMA/R_CUT, 3));
	r[0]=res;
	r[1]=virial;
	return r;
}
*/
double gb_total_energy() {
	double res = 0;
	int i, j;
	for (i=0; i<ATOMS; i++) {
		for (j=i+1; j<ATOMS; j++) {
			double d = displacement(box[i].x, box[j].x, box[i].y, box[j].y, box[i].z, box[j].z);
			//cout << d << "\t"; 
			if (d < sig(box[i], box[j])-sigma0) {
				return 10000000000;
			}
			if (d<R_CUT) {
				//cout << sig(box[i],box[j]) << "\t" << sigma0 << endl;
				double attraction = pow((sigma0/(d-sig(box[i],box[j])+sigma0)), 6);
				res+=4.0*eps(box[i],box[j])*(attraction*attraction-attraction);
			}
			//cout << res << endl;
		}
	}
	//add tail correction
	return res;
}

double pressure_calc() {
	double res = 0;
	int i,j;
	for (i=0; i<ATOMS; i++) {
		for (j=i+1; j<ATOMS; j++) {
			double dx = box[i].x - box[j].x;
			double dy = box[i].y - box[j].y;
			double dz = box[i].z - box[j].z;
			double mi_dx, mi_dy, mi_dz;
			mi_dx=mi_sign(dx);
			mi_dy=mi_sign(dy);
			mi_dz=mi_sign(dz);
			double d_sqr = mi_dx*mi_dx + mi_dy*mi_dy + mi_dz*mi_dz;
			double sig_sqr = SIGMA*SIGMA;
			//double repulsion = pow(sig_sqr/d_sqr, 6);
			double attraction = pow((sigma0/(sqrt(d_sqr)-sig(box[i],box[j])+sigma0)), 6);
			if (d_sqr<(R_CUT*R_CUT)) {
				double frc = (48.0*eps(box[i],box[j])/d_sqr)*(attraction*attraction - .5*attraction);
				double f_x, f_y, f_z;
				f_x = frc*mi_dx;
				f_y = frc*mi_dy;
				f_z = frc*mi_dz;
				//calculate virial here
				res+=f_x*dx+f_y*dy+f_z*dz;
			}
		}
	}
	//adding the tail correction due to R_CUT
	res /=(3*BL*BL*BL);
	res+=D*KB*TEMP;
	//adding tail correction
	//res+=(16.0/3.0)*M_PI*D*D*EPSILON*SIGMA*SIGMA*SIGMA*((2.0/3.0)*pow(SIGMA/R_CUT, 9)-pow(SIGMA/R_CUT,3));
	return res;
}

double gb_force_update() {
	double virial = 0;
	int i,j;
	for (i=0; i<ATOMS; i++) {
		for (j=i+1; j<ATOMS; j++) {
			double dx = box[i].x - box[j].x;
			double dy = box[i].y - box[j].y;
			double dz = box[i].z - box[j].z;
			double mi_dx, mi_dy, mi_dz;
			mi_dx=mi_sign(dx);
			mi_dy=mi_sign(dy);
			mi_dz=mi_sign(dz);
			double d_sqr = mi_dx*mi_dx + mi_dy*mi_dy + mi_dz*mi_dz;
			double attraction = pow((sigma0/(sqrt(d_sqr)-sig(box[i],box[j])+sigma0)), 6);
			if (d_sqr<(R_CUT*R_CUT)) {
				double frc = (48.0*eps(box[i],box[j])/d_sqr)*(attraction*attraction - .5*attraction);
				double f_x, f_y, f_z;
				f_x = frc*mi_dx;
				f_y = frc*mi_dy;
				f_z = frc*mi_dz;
				//calculate virial here
				box[i].set_force(box[i].fx+f_x, box[i].fy+f_y, box[i].fz+f_z);
				box[j].set_force(box[j].fx-f_x, box[j].fy-f_y, box[j].fz-f_z);
				virial+=f_x*dx+f_y*dy+f_z*dz;
			}
		}
	}
	return virial;
}

//retur probability if a trial move is to be accepted for Gay_Berne potential
//REMEMBER THAT THE ORIENTATION ALSO NEEDS TO CHANGE
double gb_accept(int particle, Atom test) {
	int i;
	double energy_change = 0.0;
	for (i=0; i<ATOMS; i++) {
		if (i!=particle) {
			double dis_new = displacement(box[i].x, test.x, box[i].y, test.y, box[i].z, test.z);
			if (dis_new < sig(box[i], test)-sigma0) {
				return -1.0;
			}
			double dis_old = displacement(box[i].x, box[particle].x, box[i].y, box[particle].y, box[i].z, box[particle].z);
			double attraction_new = pow((sigma0/(dis_new-sig(test, box[i])+sigma0)), 6);
			double attraction_old = pow((sigma0/(dis_old-sig(box[particle], box[i])+sigma0)), 6);
			energy_change += 4.0*(eps(test,box[i])*(attraction_new*attraction_new-attraction_new)-eps(box[particle],box[i])*(attraction_old*attraction_old-attraction_old));
		}
	}
	return exp(-1.0*energy_change/(KB*TEMP));
}

//generates a random movement for random particle and moves if accepted (returns 1 if accepted, 0 if not)
int mc_move() {
	int rand_particle = (ran->uniform())*ATOMS;
	double move = ran->uniform();
	double accpt = ran->uniform();
	Atom test=box[rand_particle];
	if (move>.5) {
		double rand_x = pbc(box[rand_particle].x + MAX_MOVE * (((ran->uniform())*2)-1));
		double rand_y = pbc(box[rand_particle].y + MAX_MOVE * (((ran->uniform())*2)-1));
		double rand_z = pbc(box[rand_particle].z + MAX_MOVE * (((ran->uniform())*2)-1));
		test.set_pos(rand_x, rand_y, rand_z);
	}
	else {
		//rotate by picking random points on a sphere
		test.ox=0;
		test.oy=0;
		test.oz=0;
		while(dot_product(box[rand_particle].ox, test.ox, box[rand_particle].oy, test.oy, box[rand_particle].oz, test.oz)<.8) {
			double r1,r2;
			double r3 = 1;
			while(r3>=1) {
				r1 = (((ran->uniform())*2)-1);
				r2 = (((ran->uniform())*2)-1);
				r3 = r1*r1 + r2*r2;
			}
			double u = 2*sqrt(1.0-r3);
			test.set_orientvec(u*r1, u*r2, 1.0-2.0*r3);
		}
	}
	if (accpt<gb_accept(rand_particle, test)) {
		box[rand_particle].set_pos(test.x, test.y, test.z);
		box[rand_particle].set_orientvec(test.ox, test.oy, test.oz);
		return 1;
	}
	return 0;
}


//runs equilibrations cycles and specified number of data cycles
void simulate(int cycles) {
	//equilibration cycles
	double scale = sqrt(3);
	int tot = 0;
	double accepts = 0;
	double pressav = 0.0;
	int i, j;
	for (i=0; i<EQ_CYCLES; i++) {
		for (j=0; j<ATOMS; j++) {
			mc_move();
		}
	}
	//file for VMD visualization
	ofstream sim;
	sim.open("sim.xyz"); 
	ofstream pot_press;
	pot_press.open("pot_press.dat");
	pot_press << "#Data for potential energy and instantaneous pressure as function of mc cycle." << endl;
	pot_press << "#cycle potential pressure" << endl;
	//data cycles
	for (i=0; i<cycles; i++) {
		for (j=0; j<ATOMS; j++) {
			accepts+=mc_move();
			//uncomment the following few lines if visualization data is wanted
			//(it greatly slows down the program)
		}
		if (i%SNAP == 0) {
			tot+=1;
			/*
			sim << 3*ATOMS << "\n\n";
			int z;
			for(z=0; z<ATOMS; z++) {
				sim << "He " << 10*box[z].x << " " << 10*box[z].y << " " << 10*box[z].z << endl;
				sim << "He " << 10*box[z].x + SCALE*box[z].ox << " " << 10*box[z].y + SCALE*box[z].oy << " " << 10*box[z].z + SCALE*box[z].oz << endl;
				sim << "He " << 10*box[z].x - SCALE*box[z].ox << " " << 10*box[z].y - SCALE*box[z].oy << " " << 10*box[z].z - SCALE*box[z].oz << endl;

			}
			*/
			int z;
			//order parameter/Q tensor calculations
			int r,c;
			for (z=0; z<ATOMS; z++) {
				for (r=0; r<3; r++) {
					for (c=0; c<3; c++) {
						Q[r][c] += box[z].get_orientvec(r)*box[z].get_orientvec(c) * 1.5;
					}
					Q[r][r]-=.5;
				}
			}
			//radial distribution function
			int k,l;
			for (k=0; k<ATOMS-1; k++) {
				for (l=k+1; l<ATOMS; l++) {
					double disp = displacement(box[k].x, box[l].x, box[k].y, box[l].y, box[k].z, box[l].z);
					if (disp<=K) {
						//cout << "b";
						rdf[(int)floor(disp/(K/M))]+=1.0;
						//cout << "a";
					}
				}
			}
			//makes file with final positions
			double pot = gb_total_energy();
			double press = gb_force_update();
			//double press = pressure_calc();
			//int l;
			for (l=0; l<ATOMS; l++) {
				press+=box[l].fx * box[l].x + box[l].fy * box[l].y + box[l].fz * box[l].z;
			}
			press/=(3*BL*BL*BL);
			press+=D*KB*TEMP;
			pot_press << i << " " << pot << " " << press <<  endl;
			pressav+=press;
		}
		//radial distribution function is calculated after every cycle
		
	}
	cout << "press av: " << pressav/tot << endl;
	//normalize each component
	int r,c;
	for (r=0; r<3; r++) {
		for (c=0; c<3; c++) {
			Q[r][c]/= (tot*ATOMS);
		}
	}
	//normalization for the RDF bins
	
	int m;
	for (m=0; m<M; m++) {
		if (rdf[m] != 0) {
			rdf[m]/=tot;
			rdf[m]/=pow(ATOMS, 1)/2;
			rdf[m]/=(M_PI*(pow((m+1)*(K/M), 3) - pow(m*(K/M), 3)));
		}
	}
	
	ofstream q;
	q.open("Q.dat");
	for (r=0; r<3; r++) {
		for (c=0; c<3; c++) {
			q << Q[r][c];
			if (c!=2) {
				q << " ";
			}
		}
		q << endl;
	}
	sim.close();
	pot_press.close();
	ofstream acc;
	//FILE *fp = fopen("acc_per.dat", "a");
	//fprintf(fp, "%lf, %lf, %d, %lf, ", TEMP, D, SEED, 100*accepts/(double(ATOMS)*cycles));
	acc.open("acc_per2.dat");
	//print out the acceptance ratio to ensure it is a reasonable number
	acc << "accept % "<< 100*accepts/(double(ATOMS)*cycles) << "%" << endl;
}

//TESTING FUNCTIONS
//returns true if any particles overlap
bool overlap() {
	int i,j;
	for(i=0; i<ATOMS; i++) {
		for(j=i+1; j<ATOMS; j++) {
			if (displacement(box[i].x, box[j].x, box[i].y, box[j].y, box[i].z, box[j].z) < D0) {
				return true;
			}
		}
	}
	return false;
}

//returns true if all particles are in the box
bool in_box() {
	int i;
	for (i=0; i<ATOMS; i++) {
		if (box[i].x<0 || box[i].y<0 || box[i].z<0 || box[i].x>BL || box[i].y>BL || box[i].z>BL) {
			return false;
		}
	}
	return true;
}

//PRINTS ALL OF THE NET FORCES ON THE PARTICLES
bool net_force() {
	int i;
	for (i=0; i<ATOMS; i++) {
		cout << box[i].fx << " " << box[i].fy << " " << box[i].fz << endl;
		//if (box[i].fx>.1 || box[i].fx<-.1 || box[i].fy>.1 || box[i].fy<-.1 || box[i].fz>.1 || box[i].fz<-.1) {
		//if (box[i].fx!=0 || box[i].fy!=0 || box[i].fz!=0) {
		//	return false;
		//}
	}
	return true;
}

void orient_norm() {
	int i;
	for (i=0; i<ATOMS; i++) {
		cout << sqrt(dot_product(box[i].ox, box[i].ox, box[i].oy, box[i].oy, box[i].oz, box[i].oz)) << endl;
	}
}