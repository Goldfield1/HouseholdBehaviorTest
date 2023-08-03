#define MAIN
#include "myheader.h"


typedef struct {
    double power;
    double C_tot;

    par_struct *par;

} solver_precompute_struct;


/////////////
// 5. MAIN //
/////////////
EXPORT void solve(sol_struct *sol, par_struct *par){
    
    // loop backwards
    for (int t = par->T-1; t >= 0; t--){
        solution::solve_single(t,sol,par);
    }
    
}

// use Python for simulation




