#include <iostream>
#include <vector>
#include "ppl_wrap.hh"


/*  Format: C * x + d >= 0  */

ppl_wrap::ppl_wrap()
{
}

#define FPOINT_FACTOR 100000000UL
ppl_wrap::ppl_wrap(double *Points, int rows_P, int cols_P):rows_P(rows_P),
cols_P(cols_P)
{

	int i, j, k;

	// Convert the standard double representation to a much
	// accurate representation using GMP library
	vector < mpf_class > *P1_row =
	    double2mpzclass(Points, &rows_P, &cols_P);
	vector < vector < mpf_class > >*P1 = vec2mat(P1_row, &rows_P, &cols_P);
	vector < mpz_class > single_point;

	list < Variable > VarList;
	list < Variable >::iterator it;
	Linear_Expression lin_exp[rows_P];
#ifdef DEBUG
		cout << "Cols num: " << cols_P << endl;
#endif

	for (i = 0; i < cols_P; i++)	// create a list of n_col_C "Variable" objects
		VarList.push_back(Variable(i));

	C_Polyhedron ph(cols_P, EMPTY);

#ifdef DEBUG
		cout << "Rows num: " << rows_P << endl;
#endif

	for (j = 0; j < rows_P; j++) {
#ifdef DEBUG
    cout << " j = " << j << "\n";
#endif
		single_point.clear();

		/*Create linear expression for expressing the point */
#ifdef DEBUG
    cout << " k = ";
#endif
		for (k = 0; k < cols_P; k++) {	// pick the needed row of matrix C
#ifdef DEBUG
      cout << "(" << k << ") " << (*P1)[k][j] << " ";
#endif
      /* single_point is integer only */
			//single_point.push_back((*P1)[k][j] * FPOINT_FACTOR);
			single_point.push_back((*P1)[k][j]);

		}
#ifdef DEBUG
		cout << endl;
#endif

		int m;
		for (it = VarList.begin(); it != VarList.end(); it++)	// creates linear expression "lin_exp_temp" by using the iterator of the Variable-List
		{
			m = it->id();
#ifdef DEBUG
      cout << " m = " << m << "\n";
#endif
			lin_exp[j] += single_point[m] * (*it);
#ifdef DEBUG
    cout << "single_point[ " << m << "] = " << single_point[m] << endl;
#endif
    }

		//ph.add_generator(point(lin_exp[j], FPOINT_FACTOR));
		ph.add_generator(point(lin_exp[j]));
	}
	P = ph;
  //P.drop_some_non_integer_points();
}

ppl_wrap::ppl_wrap(double *C_d, double *d_d, int rows_C, int cols_C,
		   int rows_d):rows_C(rows_C), cols_C(cols_C), rows_d(rows_d)
{

	int cols_d = 1;
	vector < mpf_class > *C1_row = double2mpzclass(C_d, &rows_C, &cols_C);
	vector < mpf_class > *d1 = double2mpzclass(d_d, &rows_d, &cols_d);
	vector < vector < mpf_class > >*C1 = vec2mat(C1_row, &rows_C, &cols_C);
	vector < mpf_class > C_row_temp;

	list < Variable > VarList;
	list < Variable >::iterator it;

	Constraint_System CS_P;
	Linear_Expression lin_exp[rows_C];

	int i, k, j;

	for (i = 0; i < cols_C; i++)	// create a list of n_col_C "Variable" objects
		VarList.push_back(Variable(i));

	for (j = 0; j < rows_C; j++)	// wrapping loop to create array of linear expressions "lin_exp"
	{
		C_row_temp.clear();	//C_row_temp has to be cleared in every iteration
		Linear_Expression lin_exp_temp;

		for (k = 0; k < cols_C; k++) {	// pick the needed row of matrix C
			C_row_temp.push_back((*C1)[k][j]);
		}

		int m;
		for (it = VarList.begin(); it != VarList.end(); it++)	// creates linear expression "lin_exp_temp" by using the iterator of the Variable-List
		{
			m = it->id();
			lin_exp_temp += C_row_temp[m] * (*it);
		}
		lin_exp[j] = lin_exp_temp;	// creates array of linear expressions "lin_exp"
	}

	for (int i = 0; i < rows_d; i++)	// creates constraint system "CS_P" with "lin_exp" and vector "d"
		CS_P.insert(lin_exp[i] + d1->at(i) >= 0);
	
  P = C_Polyhedron(CS_P);	// creates Polytop "P" with constraint system "CS_P"
  //P.drop_some_non_integer_points();
}

C_Polyhedron ppl_wrap::get_Polytope()
{
	return P;		// Returns the private property of the ppl_wrap class
}

void ppl_wrap::intersec_assign(C_Polyhedron P2)
{
	P.intersection_assign(P2);	// Computes an intersection of P and P2 and stores it in P.
  P.minimized_constraints();
}

struct pplw_hs ppl_wrap::compute_CD()
{

	unsigned int i = 0;
	struct pplw_hs temp_hfspace;
	vector < vector < mpf_class > >C;
	vector < mpf_class > d;

	Constraint_System ConSys_P;
	//ConSys_P = P.constraints();
	ConSys_P = P.minimized_constraints();
	Constraint_System::const_iterator cs_it;

	vector < mpf_class > C_row;
	mpf_class C_element, d_element;

	C.clear();		// C and d must be cleared before filled with new content
	d.clear();

	for (cs_it = ConSys_P.begin(); cs_it != ConSys_P.end(); cs_it++)	// create C matrix and d vector of constraint system
	{
		C_row.clear();
		for (i = 0; i < (*cs_it).space_dimension(); i++) {	// stores the coefficients of the linear expressoins in C_row
			C_element = (*cs_it).coefficient(Variable(i));
			C_row.push_back(C_element);
		}
		C.push_back(C_row);
		d_element = (*cs_it).inhomogeneous_term();	// the d vector is build of the inhomogeneous part of the constraint
		d.push_back(d_element);
		
		/* Christian: It can happen that the constrain is A <= b and A >= b.
		 * Then PPL wraps this into A = b. In this case we automatically add A <= b again.
		 * But we still need to add -A <= -b. This was fixed.
		 */
		switch ((*cs_it).type()) {
			case Constraint::EQUALITY: // "=" (also add -A <= -b)
				C_row.clear();
				for (i = 0; i < (*cs_it).space_dimension(); i++) {	// stores the coefficients of the linear expressoins in C_row
					C_element = (*cs_it).coefficient(Variable(i));
					C_row.push_back(-C_element);
				}
				C.push_back(C_row);
				d.push_back(-d_element);
				break;
			// case Constraint::NONSTRICT_INEQUALITY: // "<=" (everything is fine)
			case Constraint::STRICT_INEQUALITY: // "<" (unexpected in our setting)
				cout << "PPL Error: Unexpected strict inequality!" << endl;
				cerr << "PPL Error: Unexpected strict inequality!" << endl;
				break;
		}
	}

	temp_hfspace.C = C;
	temp_hfspace.d = d;

	return temp_hfspace;
}

struct pplw_hs ppl_wrap::compute_ClosurePoints()
{

	struct pplw_hs temp_hfspace;

	vector < vector < mpf_class > >ClosurePoints;
	vector < mpf_class > CP_row;
	mpf_class CP_element;

	Generator_System gs_P;
	gs_P = P.minimized_generators();	// minimized_generators are used
	Generator_System::const_iterator gs_it;

	for (gs_it = gs_P.begin(); gs_it != gs_P.end(); gs_it++) {	// create ClosurePoints matrix
		if ((*gs_it).is_point()) {
			CP_row.clear();
			for (unsigned int i = 0; i < (*gs_it).space_dimension(); i++) {	// stores the point of every generator in CP_row
				CP_element = (*gs_it).coefficient(Variable(i));
				CP_row.push_back(CP_element /
						 (*gs_it).divisor());
			}
			ClosurePoints.push_back(CP_row);
		} else
			cout << "Generator is not a point:" << endl;
	}

	temp_hfspace.ClosurePoints = ClosurePoints;

	return temp_hfspace;
}

vector < vector < mpf_class > >*ppl_wrap::vec2mat(vector < mpf_class > *C,
						  int *rows, int *cols)
{
	int i = 0;
	vector < mpf_class > C_element;
	vector < vector < mpf_class > >*C_mat =
	    new vector < vector < mpf_class > >;
	for (int spalte = 0; spalte < (*cols); spalte++) {
		for (int zeile = 0; zeile < (*rows); zeile++) {
			C_element.push_back(C->at(i));
			i++;
		}
		C_mat->push_back(C_element);
		C_element.clear();
	}

	return C_mat;
}

vector < mpf_class > *ppl_wrap::double2mpzclass(double *D, int *rows, int *cols)
{

	int i = 0;
	mpf_t dummy;
	mpf_init(dummy);
	mpf_class dummy_class;
	vector < mpf_class > *dummy_vector = new vector < mpf_class >;

	for (i = 0; i < ((*cols) * (*rows)); i++) {
		mpf_set_d(dummy, D[i]);
		dummy_class = mpf_class(dummy);
		dummy_vector->push_back(dummy_class);
	}

	mpf_clear(dummy);
	return dummy_vector;
}

void ppl_wrap::remove_dimensions(double *dimensions)
{
	Variables_Set vars;

	int i = 0;
	for (i = 0; i < cols_C; i++) {
		if ((i != dimensions[0] - 1) && (i != dimensions[1] - 1))
			vars.insert(Variable(i));
	}

	P.remove_space_dimensions(vars);
}
