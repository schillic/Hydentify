#include <iostream>
#include <ppl.hh>
#include <mex.h>
#include "ppl_wrap.hh"

using namespace std;
using namespace Parma_Polyhedra_Library;

struct cmd_dispatch_str {
  const char* cmd;
  int arg_in_min;
  int arg_in_max;
  int arg_out_min;
  int arg_out_max;
  void (*handler)(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
};

void do_points2polytope(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_mappingdim(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_isempty(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_extremepoints(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_intersection(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_printpolytope(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_contains(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_isdisjoint(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void do_equals(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);

struct cmd_dispatch_str cmd_tab[] =
{
  {
    "Points2Polytope", 2, 2, 2, 2, do_points2polytope
  },
  {
    "ExtremePoints", 3, 3, 1, 1, do_extremepoints
  },
  {
    "Intersection", 5, 5, 2, 2, do_intersection
  },
  {
    "IsEmpty", 3, 3, 1, 1, do_isempty
  },
  {
    "mapping_dim", 4, 4, 2, 2, do_mappingdim
  },
  {
    "PrintPolytope", 3, 3, 0, 0, do_printpolytope
  },
  {
    "Contains", 5, 5, 1, 1, do_contains
  },
  {
    "IsDisjoint", 5, 5, 1, 1, do_isdisjoint
  },
  {
    "Equals", 5, 5, 1, 1, do_equals
  },
};

struct pplw_hs
intersection(double *C1_d, double *d1_d, int rows_C1, int cols_C1,
	     int rows_d1, double *C2_d, double *d2_d, int rows_C2,
	     int cols_C2, int rows_d2)
{

	ppl_wrap P1_wrap(C1_d, d1_d, rows_C1, cols_C1, rows_d1);
	ppl_wrap P2_wrap(C2_d, d2_d, rows_C2, cols_C2, rows_d2);

	C_Polyhedron P2 = P2_wrap.get_Polytope();
	C_Polyhedron P1 = P1_wrap.get_Polytope();

	dimension_type dim_P1 = P1.space_dimension(), dim_P2 =
	    P2.space_dimension();
	struct pplw_hs hfspace;

	//    check if P1 & P2 have same dimension
	if (dim_P1 == dim_P2) {
		// P1 = intersection of P1 & P2
		P1_wrap.intersec_assign(P2);
    P1_wrap.get_Polytope().is_empty();
		hfspace = P1_wrap.compute_CD();
	} else {
		cout << "Polyhedron A and Polyhedron B must have the same dimension! \n";
	}

	return hfspace;
}


int check_sane_Cd(int rows_c, int rows_d)
{
  return rows_c == rows_d;
}

void do_points2polytope(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_Points = mxGetM(prhs[1]);
  int cols_Points = mxGetN(prhs[1]);

  /* Assign pointers to each input */
  double *Points = mxGetPr(prhs[1]);
  ppl_wrap P1_wrap(Points, rows_Points, cols_Points);
  P1_wrap.get_Polytope().is_empty();

  struct pplw_hs hfspace = P1_wrap.compute_CD();
  /* Create matrix for the return argument. */

  plhs[0] = mxCreateDoubleMatrix((int)hfspace.C.size(), (int)hfspace.C.
        at(0).size(), mxREAL);
  plhs[1] = mxCreateDoubleMatrix((int)hfspace.d.size(), 1, mxREAL);

  double *C_new = mxGetPr(plhs[0]), *d_new = mxGetPr(plhs[1]);

  int k = 0, i = 0, j = 0;
  for (i = 0; i < (int)hfspace.C.at(0).size(); i++) {
    for (j = 0; j < (int)hfspace.C.size(); j++) {
      C_new[k] = hfspace.C.at(j).at(i).get_d();
      if (k < (int)hfspace.d.size()) {
        d_new[k] = hfspace.d.at(j).get_d();
      }
      k++;
    }
  }
}

void do_extremepoints(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C = mxGetM(prhs[1]);
  int cols_C = mxGetN(prhs[1]);
  int rows_d = mxGetM(prhs[2]);

  if(!check_sane_Cd(rows_C, rows_d))
  {
    mexErrMsgTxt("C rows and d rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C_d = mxGetPr(prhs[1]);
  double *d_d = mxGetPr(prhs[2]);

  ppl_wrap P1_wrap(C_d, d_d, rows_C, cols_C,
      rows_d);
  P1_wrap.get_Polytope().is_empty(); /* minimize it */
  struct pplw_hs hfspace = P1_wrap.compute_ClosurePoints();

  /* Create matrix for the return argument. */
  if (hfspace.ClosurePoints.size() > 0)
  {
    plhs[0] = mxCreateDoubleMatrix((int) hfspace.ClosurePoints.size(),
        (int) hfspace.ClosurePoints.at(0).size(), mxREAL);
    double *ExtrPoints = mxGetPr(plhs[0]);
    int k = 0, i = 0, j = 0;

    for (i = 0; i < (int)hfspace.ClosurePoints.at(0).size(); i++) {
      for (j = 0; j < (int)hfspace.ClosurePoints.size(); j++) {
        ExtrPoints[k] = hfspace.ClosurePoints.at(j).at(i).get_d();
        k++;
      }
    }
  }
  else
  {
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
  }
}

void do_isempty(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C = mxGetM(prhs[1]);
  int cols_C = mxGetN(prhs[1]);
  int rows_d = mxGetM(prhs[2]);

  if(!check_sane_Cd(rows_C, rows_d))
  {
    mexErrMsgTxt("C rows and d rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C_d = mxGetPr(prhs[1]);
  double *d_d = mxGetPr(prhs[2]);
  mxLogical bool_val;


  ppl_wrap P1_wrap(C_d, d_d, rows_C,
      cols_C, rows_d);
  C_Polyhedron P = P1_wrap.get_Polytope();
  P.minimized_constraints();
  bool_val = P.is_empty();
  plhs[0] = mxCreateLogicalScalar(bool_val);
}

void do_intersection(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C1 = mxGetM(prhs[1]);
  int rows_C2 = mxGetM(prhs[2]);
  int rows_d1 = mxGetM(prhs[3]);
  int rows_d2 = mxGetM(prhs[4]);
  int cols_C1 = mxGetN(prhs[1]);
  int cols_C2 = mxGetN(prhs[2]);
  
  if(!check_sane_Cd(rows_C1, rows_d1))
  {
    mexErrMsgTxt("C1 rows and d1 rows are not equal: check polytope description!");
  }


  if(!check_sane_Cd(rows_C2, rows_d2))
  {
    mexErrMsgTxt("C2 rows and d2 rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C1_d = mxGetPr(prhs[1]);
  double *C2_d = mxGetPr(prhs[2]);
  double *d1_d = mxGetPr(prhs[3]);
  double *d2_d = mxGetPr(prhs[4]);

  struct pplw_hs hfspace =
    intersection(C1_d, d1_d, rows_C1, cols_C1,
        rows_d1, C2_d,
        d2_d, rows_C2, cols_C2,
        rows_d2);

  /* Create matrix for the return argument. */
  plhs[0] =
    mxCreateDoubleMatrix((int)hfspace.C.size(),
        (int)hfspace.C.
        at(0).size(), mxREAL);
  plhs[1] =
    mxCreateDoubleMatrix((int)hfspace.d.size(),
        1, mxREAL);

  double *C_new = mxGetPr(plhs[0]), *d_new =
    mxGetPr(plhs[1]);

  int k = 0, i = 0, j = 0;
  for (i = 0; i < (int)hfspace.C.at(0).size();
      i++) {
    for (j = 0; j < (int)hfspace.C.size();
        j++) {
      C_new[k] =
        hfspace.C.at(j).
        at(i).get_d();
      if (k < (int)hfspace.d.size()) {
        d_new[k] =
          hfspace.d.
          at(j).get_d();
      }
      k++;
    }
  }
}

void do_contains(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C1 = mxGetM(prhs[1]);
  int rows_C2 = mxGetM(prhs[2]);
  int rows_d1 = mxGetM(prhs[3]);
  int rows_d2 = mxGetM(prhs[4]);
  int cols_C1 = mxGetN(prhs[1]);
  int cols_C2 = mxGetN(prhs[2]);
  mxLogical ret;
  
  if(!check_sane_Cd(rows_C1, rows_d1))
  {
    mexErrMsgTxt("C1 rows and d1 rows are not equal: check polytope description!");
  }


  if(!check_sane_Cd(rows_C2, rows_d2))
  {
    mexErrMsgTxt("C2 rows and d2 rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C1_d = mxGetPr(prhs[1]);
  double *C2_d = mxGetPr(prhs[2]);
  double *d1_d = mxGetPr(prhs[3]);
  double *d2_d = mxGetPr(prhs[4]);

  ppl_wrap p1(C1_d, d1_d, rows_C1,
      cols_C1, rows_d1);
  ppl_wrap p2(C2_d, d2_d, rows_C2,
      cols_C2, rows_d2);

  ret = p1.get_Polytope().contains(p2.get_Polytope());

  //plhs[0] = mxCreateDoubleScalar(ret);
  plhs[0] = mxCreateLogicalScalar(ret);
}

void do_isdisjoint(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C1 = mxGetM(prhs[1]);
  int rows_C2 = mxGetM(prhs[2]);
  int rows_d1 = mxGetM(prhs[3]);
  int rows_d2 = mxGetM(prhs[4]);
  int cols_C1 = mxGetN(prhs[1]);
  int cols_C2 = mxGetN(prhs[2]);
  mxLogical ret;
  
  if(!check_sane_Cd(rows_C1, rows_d1))
  {
    mexErrMsgTxt("C1 rows and d1 rows are not equal: check polytope description!");
  }


  if(!check_sane_Cd(rows_C2, rows_d2))
  {
    mexErrMsgTxt("C2 rows and d2 rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C1_d = mxGetPr(prhs[1]);
  double *C2_d = mxGetPr(prhs[2]);
  double *d1_d = mxGetPr(prhs[3]);
  double *d2_d = mxGetPr(prhs[4]);

  ppl_wrap p1(C1_d, d1_d, rows_C1,
      cols_C1, rows_d1);
  ppl_wrap p2(C2_d, d2_d, rows_C2,
      cols_C2, rows_d2);

  ret = p1.get_Polytope().is_disjoint_from(p2.get_Polytope());

  //plhs[0] = mxCreateDoubleScalar(ret);
  plhs[0] = mxCreateLogicalScalar(ret);
}

void do_equals(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C1 = mxGetM(prhs[1]);
  int rows_C2 = mxGetM(prhs[2]);
  int rows_d1 = mxGetM(prhs[3]);
  int rows_d2 = mxGetM(prhs[4]);
  int cols_C1 = mxGetN(prhs[1]);
  int cols_C2 = mxGetN(prhs[2]);
  mxLogical ret;
  
  if(!check_sane_Cd(rows_C1, rows_d1))
  {
    mexErrMsgTxt("C1 rows and d1 rows are not equal: check polytope description!");
  }


  if(!check_sane_Cd(rows_C2, rows_d2))
  {
    mexErrMsgTxt("C2 rows and d2 rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C1_d = mxGetPr(prhs[1]);
  double *C2_d = mxGetPr(prhs[2]);
  double *d1_d = mxGetPr(prhs[3]);
  double *d2_d = mxGetPr(prhs[4]);

  ppl_wrap p1(C1_d, d1_d, rows_C1,
      cols_C1, rows_d1);
  ppl_wrap p2(C2_d, d2_d, rows_C2,
      cols_C2, rows_d2);

  ret = p1.get_Polytope()==p2.get_Polytope();
  /*
  ret = 0;
  if(p1.get_Polytope().space_dimension() == p2.get_Polytope().space_dimension())
    ret = p1.get_Polytope().contains(p2.get_Polytope()) &&
      p2.get_Polytope().contains(p1.get_Polytope());
*/
  //plhs[0] = mxCreateDoubleScalar(ret);
  plhs[0] = mxCreateLogicalScalar(ret);
}

void do_mappingdim(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C = mxGetM(prhs[1]);
  int cols_C = mxGetN(prhs[1]);
  int rows_d = mxGetM(prhs[2]);

  if(!check_sane_Cd(rows_C, rows_d))
  {
    mexErrMsgTxt("C rows and d rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C_d = mxGetPr(prhs[1]);
  double *d_d = mxGetPr(prhs[2]);
  double *dimensions = mxGetPr(prhs[3]);

  ppl_wrap P1_wrap(C_d, d_d, rows_C, cols_C,
      rows_d);
  P1_wrap.remove_dimensions(dimensions);

  struct pplw_hs hfspace = P1_wrap.compute_CD();

  /* Create matrix for the return argument. */
  plhs[0] =
    mxCreateDoubleMatrix((int)hfspace.C.size(),
        (int)hfspace.C.
        at(0).size(), mxREAL);
  plhs[1] =
    mxCreateDoubleMatrix((int)hfspace.d.size(),
        1, mxREAL);

  double *C_new = mxGetPr(plhs[0]), *d_new =
    mxGetPr(plhs[1]);

  int k = 0, i = 0, j = 0;
  for (i = 0; i < (int)hfspace.C.at(0).size();
      i++) {
    for (j = 0; j < (int)hfspace.C.size();
        j++) {
      C_new[k] =
        hfspace.C.at(j).
        at(i).get_d();
      if (k < (int)hfspace.d.size()) {
        d_new[k] =
          hfspace.d.
          at(j).get_d();
      }
      k++;
    }
  }
}

void do_printpolytope(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  /* Get numbers of rows and cols of input data */
  int rows_C = mxGetM(prhs[1]);
  int cols_C = mxGetN(prhs[1]);
  int rows_d = mxGetM(prhs[2]);
  
  if(!check_sane_Cd(rows_C, rows_d))
  {
    mexErrMsgTxt("C rows and d rows are not equal: check polytope description!");
  }

  /* Assign pointers to each input */
  double *C_d = mxGetPr(prhs[1]);
  double *d_d = mxGetPr(prhs[2]);


  ppl_wrap P1_wrap(C_d, d_d, rows_C,
      cols_C, rows_d);
  C_Polyhedron P = P1_wrap.get_Polytope();
  P.is_empty();
  P.ascii_dump(cout);
  cout << endl << " ********************************************* " << endl;
  P.constraints().print();
  cerr << endl;
}


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
#define MAX_INPUT_SIZE    50
  int i,j;
  char input[MAX_INPUT_SIZE];
  if (mxGetString(prhs[0], input, MAX_INPUT_SIZE-1) != 0) {
    mexErrMsgTxt("The first argument must be the command string.");
  }

  for(i = 0; i < sizeof(cmd_tab)/sizeof(struct cmd_dispatch_str); ++i)
  {
    if (strncmp(input, cmd_tab[i].cmd, MAX_INPUT_SIZE-1) == 0)
    {
      if (cmd_tab[i].arg_in_min < nrhs || cmd_tab[i].arg_in_max > nrhs)
      {
        mexErrMsgTxt("Check number of input parameters.");
      }
      if (cmd_tab[i].arg_out_min < nlhs || cmd_tab[i].arg_out_max > nlhs)
      {
        mexErrMsgTxt("Check number of output parameters.");
      }
      cmd_tab[i].handler(nlhs, plhs, nrhs, prhs);
      return;
    }
  }

  mexErrMsgTxt("Unknown command. Check first argument.");
}
