#define MAIN
#include "myheader.h"


typedef struct {
    double power;
    double C_tot;

    par_struct *par;

} solver_precompute_struct;

/* 
double objfunc_precompute(unsigned n, const double *x, double *grad, void *solver_data_in){
    // unpack


}
*/

/////////////
// 5. MAIN //
/////////////


EXPORT void solve(sol_struct *sol, par_struct *par){
    
    
    // loop backwards
    //for (int t = par->T-1; t >= 0; t--){
    for (int t = par->T-1; t >= 0; t--){
        //par->r = 0.03;
        single::solve_single(t,sol,par);
        //single::solve_vfi(t,sol,par);
    }
    
}




EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}


double wwage_func(int t, double capital, par_struct* par) {

    return (1.0 - par->tau )* par->w_vec[t];// * (1.0 + par->alpha * capital);
}


EXPORT double test(sol_struct *sol, par_struct *par, double capital){
        // flow-utility

        float c = 0.5;
        float h = 1.0;
        int t = 6;
        float a = 1.0;

        double* V_next = &sol->V[index::index3(t+1,0,0,par->T,par->Na,par->Nk)];

        double Util = utils::util(c,h,par);
        
        // continuation value
        double *a_grid = par->a_grid;
        double *k_grid = par->k_grid;
        
        double pen = 0;
        if (c < 0) {
            pen += c * 1000;
            c = 0.0000001;
        }
        if (h < 0) {
            pen += h * 1000;
            h = 0.0;
        }

        double income = h * wwage_func(t, capital, par);
        double a_next = (1 + par->r) * (a + income - c);
        double k_next = capital + h;
        //k_next = 0;
    

        double _Vnext = tools::interp_2d(a_grid,k_grid, par->Na,par->Nk, V_next, a_next,k_next);

        //double _Vnext = tools::interp_1d(a_grid,k_grid, par->Na,par->Nk, V_next, a_next,k_next);
        

        // return discounted sum
        return _Vnext;
};
