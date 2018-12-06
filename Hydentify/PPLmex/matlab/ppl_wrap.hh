#ifndef PPL_WRAP_H_
#define PPL_WRAP_H_

#include <vector>
#include <gmpxx.h>
#include <ppl.hh>

using namespace std;
using namespace Parma_Polyhedra_Library;

struct pplw_hs {
	vector< vector<mpf_class> > C;
	vector< mpf_class> d;
	vector< vector<mpf_class> > ClosurePoints;
};

class ppl_wrap {
 public:
	ppl_wrap();
	ppl_wrap(double *C1_d, double *d1_d, int rows_C1, int cols_C1,
		 int rows_d1);
	ppl_wrap(double *P, int rows_P, int cols_P);

	struct pplw_hs compute_CD();
	struct pplw_hs compute_ClosurePoints();
	C_Polyhedron get_Polytope();
	void intersec_assign(C_Polyhedron P2);
	void remove_dimensions(double *dimensions);

 private:
	C_Polyhedron P;
	int rows_C;
	int cols_C;
	int rows_d;
	int rows_P;
	int cols_P;

	 vector< vector<mpf_class> >* vec2mat(vector<mpf_class>* C,
						 int *rows, int *cols);
	 vector<mpf_class> *double2mpzclass(double *D, int *rows, int *cols);
};

#endif				/* PPL_WRAP_H_ */
