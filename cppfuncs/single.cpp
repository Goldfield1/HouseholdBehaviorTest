#ifndef MAIN
#define SINGLE
#include "myheader.cpp"
#include <tuple>
#endif

namespace single {
    typedef struct {
        
        double a;
        double capital;
        int t;
        double *V_next;
        double *c_next;
        double last_c;
        double last_h;
        par_struct *par;
        int iIn;

    } solver_single_struct;
    
    double wage_func(int t, double capital, par_struct* par) {

        return (1.0 - par->tau )* par->w_vec[t]; //* (1.0 + par->alpha * capital);
    }

    double cons_last(double hours, double assets, double capital, par_struct* par) {

        double income = hours * wage_func(par->T-1, capital, par);

        return income + assets;
    }

    double value_of_choice(double c,double h, double a, double capital, double* V_next, int t,par_struct* par){

        
        // continuation value
        double *a_grid = par->a_grid;
        double *k_grid = par->k_grid;
        
        double pen = 0;
        if (c < 0) {
            pen += c * 1000;
            c = 0.00001;
        }
        if (h < 0) {
            pen += h * 1000;
            h = 0.0;
        }

        // flow-utility
        double Util = utils::util(c,h,par);

        double income = h * wage_func(t, capital, par);
        double a_next = (1 + par->r) * (a + income - c);
        double k_next = capital + h;
        //k_next = 0;
    

        double _Vnext = tools::interp_2d(a_grid,k_grid, par->Na,par->Nk, V_next, a_next,k_next);

        //double _Vnext = tools::interp_1d(a_grid,k_grid, par->Na,par->Nk, V_next, a_next,k_next);
        

        // return discounted sum
        return Util + par->rho*_Vnext + pen;
    }

    double objfunc_last(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;
        
        par_struct *par = solver_data->par;

        double h = x[0];
        double a = solver_data->a;
        int t = solver_data->t;
        double capital = solver_data->capital;
        double c = cons_last(h, a, capital, par);
        
        double pen = 0;
        if (c < 0) {
            pen += c * 1000;
            c = 0.00001;
        }

        //return 0.5;
        return - utils::util(c,h,par); // - pen;
        //return - value_of_choice(c,h,a,capital,solver_data->V_next,t,par);

    }

    double objfunc(unsigned n, const double *x, double *grad, void *solver_data_in){
        
        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;
        
        par_struct *par = solver_data->par;
        
        double c = x[0];
        double h = x[1];
        
        double a = solver_data->a;
        int t = solver_data->t;
        double capital = solver_data->capital;

        return - value_of_choice(c,h,a,capital,solver_data->V_next,t,par);

    }

    /*
    // double* objfunc_egm(unsigned n, const double *x, double *grad, double* c_next, void *solver_data_in){
    double egm_step(double h, double* c_next, void *solver_data_in, int iIn){
    
        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;
        
        par_struct *par = solver_data->par;
        
        //double h = x[0];

        int t = solver_data->t;
        double capital = solver_data->capital;


        double a;

        // loop over a grid
        double wage = wage_func(t, capital, par);

        double fac = pow((wage / par->beta) , (1.0/par->gamma));
        
        int na = par->Na;
        double* marg_u_plus = new double[na];
        double* c_endo = new double[na];

        double* ell_endo = new double[na];
        double* m_endo = new double[na];
        double* m_exo = new double[na];

        double* _c_sol = new double[na];
        double* _h_sol = new double[na];

        int idx;
        int idxH;
        for (int iA=0; iA<par->Na;iA++){ 

            idx = index::index3(t+1,iA,0,par->T,par->Na,par->Nk);
            idxH = index::index3(t+1,iA,0,par->T,par->Na,par->Nk);

            //auto h_now = h[idxH];

            a = par->a_grid[iA] ; // / (1 + par->r);

            marg_u_plus[iA] = pow(c_next[idx] , (par->eta));
            c_endo[iA] = pow((marg_u_plus[iA]) , (1 / par->eta));
            
            m_endo[iA] = (c_endo[iA] + a / (1 + par->r) - wage * ell_endo[iA]);

            m_exo[iA] = a; /// (1 + par->r); //* (1 + par->r);

            
        } 
        for (int iA=0; iA<par->Na;iA++){ 
            //_c_sol[iA] = tools::interp_1d(m_endo,par->Na, c_endo, m_exo[iA]);
            _h_sol[iA] = tools::interp_1d(m_endo,par->Na, ell_endo, m_exo[iA]);
            _c_sol[iA] = c_endo[iA];
        }

        double aa = _c_sol[iIn];
        double bb = _h_sol[iIn];
        return {aa, bb};

    } */


    double objfunc_egm(unsigned n, const double *x, double *grad, void *solver_data_in){

        double h = x[0];
        
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;

        par_struct *par = solver_data->par;

        double a = solver_data->a;
        int t = solver_data->t;

        double capital = solver_data->capital;
        double* c_next = solver_data->c_next;
        int iIn = solver_data->iIn;

        //double c = egm_step(h, c_next, solver_data_in, iIn);
        auto c = 0;
        
        return - value_of_choice(c,h,a,capital,solver_data->V_next,t,par);
    }


    void solve_single(int t,sol_struct *sol,par_struct *par){
        // terminal period
        //int idx = index::index2(2,2,par->T,par->Na);
        //sol->c[idx] = 1.0;

        //int idx = index::index2(1,1,par->T,par->Na);
        //sol->c[idx] = 1.0;
        
        #pragma omp parallel num_threads(par->threads) 
        {

        solver_single_struct* solver_data = new solver_single_struct;
        solver_single_struct* solver_data_egm = new solver_single_struct;

        
        if (t == (par->T-1)){

            //#pragma omp parallel num_threads(par->threads) 
            //{
                // 1. allocate objects for solver
                

                int dim = 1;
                double lb[1],ub[1],x[1];

                auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                //auto opt = nlopt_create(NLOPT_LD_MMA, dim);
                double minf=0.0;

                #pragma omp for
                for (int iA=0; iA<par->Na;iA++){
                    for (int ik=0; ik<par->Nk;ik++){
                        // resources
                        double a = par->a_grid[iA];
                        double k = par->k_grid[ik];

                        int idx = index::index3(t,iA,ik,par->T,par->Na,par->Nk);

                        solver_data->a = a;
                        solver_data->capital = k;
                        solver_data->t = t;
                        //solver_data->V_next = 0
                        solver_data->par = par;
                        nlopt_set_min_objective(opt, objfunc_last, solver_data);

                        // set nlopt tolerance
                        nlopt_set_ftol_abs(opt, 1.0e-12);
                        nlopt_set_xtol_abs1(opt, 1.0e-10);

                        lb[0] = 0.0;
                        ub[0] = 50;
                        nlopt_set_lower_bounds(opt, lb);
                        nlopt_set_upper_bounds(opt, ub);

                        // optimize
                        x[0] = - a / wage_func(t, k, par) + 1.0e-8;

                        if (x[0] < 0) {
                            x[0] = 1.0;
                        } else {
                            x[0] = x[0] * 1.0002;
                        }

                        nlopt_optimize(opt, x, &minf);

                        sol->h[idx] = x[0];
                        sol->V[idx] = -minf;     
                        //sol->V[idx] = 2.0;
                        sol->c[idx] = cons_last(sol->h[idx], a, k, par);

                        // set egm also
                        sol->h_egm[idx] = x[0];
                        sol->c_egm[idx] = sol->c[idx];
                        sol->V_egm[idx] = sol->V[idx];


                    }
                }   

                nlopt_destroy(opt);
        
            //}
        } else {

            //#pragma omp parallel num_threads(par->threads) 
            //{
                // 1. allocate objects for solver
                //solver_single_struct* solver_data = new solver_single_struct;

                int dim2 = 2;
                double lb[2],ub[2],x[2];

                //auto opt2 = nlopt_create(NLOPT_LD_LBFGS, dim2); // NLOPT_LD_MMA NLOPT_LD_LBFGS     NLOPT_LN_BOBYQA
                auto opt2 = nlopt_create(NLOPT_LN_BOBYQA, dim2);
                //auto opt2 = nlopt_create(NLOPT_LN_BOBYQA, dim2);

                double lb_egm[1],ub_egm[1],x_egm[1];
                auto opt_egm = nlopt_create(NLOPT_LN_BOBYQA, 1);

                double minf=0.0;
                double minf_egm=0.0;

                //#pragma omp for
                //solver_data->last_c = 0.001;
                //solver_data->last_h = 0.001;
                #pragma omp for
                for (int iA=0; iA<par->Na;iA++){
                    for (int ik=0; ik<par->Nk;ik++){
                        // resources
                        double a = par->a_grid[iA];
                        double k = par->k_grid[ik];

                        int idx = index::index3(t,iA,ik,par->T,par->Na,par->Nk);
                        int idx_last = index::index3(t,iA-1,ik,par->T,par->Na,par->Nk);

                        solver_data->a = a;
                        solver_data->capital = k;
                        solver_data->t = t;
                        //solver_data->V_next = &sol->V[index::index3(t+1,0,0,par->T,par->Na,par->Nk)];
                        solver_data->V_next = &sol->V[index::index3(t+1,0,0,par->T,par->Na,par->Nk)];
                        solver_data->par = par;

                        nlopt_set_min_objective(opt2, objfunc, solver_data);

                        // bounds
                        lb[0] = 0.001;
                        ub[0] = 50;
                        lb[1] = 0.001;
                        ub[1] = 50;
                        nlopt_set_lower_bounds(opt2, lb);
                        nlopt_set_upper_bounds(opt2, ub);

                        nlopt_set_ftol_abs(opt2, 1.0e-12);
                        nlopt_set_xtol_abs1(opt2, 1.0e-10);


                        //nlopt_set_initial_step1(opt2, 0.0000001);


                        //nlopt_set_vector_storage(opt2, 4);


                        
                        if (ik == 0) {
                            x[0] = 0.99;
                            //x[0] = 0.00008;
                            x[1] = 1.0;
                        } else {
                            x[0] = solver_data->last_c;
                            x[1] = solver_data->last_h;
                        }

                        nlopt_optimize(opt2, x, &minf);

                        sol->h[idx] = x[1];
                        sol->c[idx] =  x[0];
                        sol->V[idx] = -minf;

                        solver_data->last_h = x[1] / 1.001;
                        solver_data->last_c = x[0] * 1.001;

                        /*
                        auto c_next = sol->c[index::index3(t+1,iA,0,par->T,par->Na,par->Nk)];

                        // objfunc_egm(double h, double* c_next, void *solver_data_in, int iIn){
                        sol->c_egm[idx] = egm_step(sol->h[index::index3(t,iA,0,par->T,par->Na,par->Nk)] ,sol->c, solver_data, iA);
                        */ 
                        // EGM //
                        /* 
                        solver_data_egm->a = a;
                        solver_data_egm->capital = k;
                        solver_data_egm->t = t;
                        solver_data_egm->iIn = iA;
                        //solver_data_egm->V_next = &sol->V[index::index3(t+1,0,0,par->T,par->Na,par->Nk)];

                        solver_data_egm->V_next = &sol->V[index::index3(t+1,0,0,par->T,par->Na,par->Nk)];
                        solver_data_egm->c_next = &sol->c[index::index3(t+1,0,ik,par->T,par->Na,par->Nk)];;

                        solver_data_egm->par = par;
                        
                        nlopt_set_min_objective(opt_egm, objfunc_egm, solver_data_egm);

                        lb_egm[0] = 0.1;
                        ub_egm[0] = 10;
                        nlopt_set_lower_bounds(opt_egm, lb_egm);
                        nlopt_set_upper_bounds(opt_egm, ub_egm);
                        nlopt_set_ftol_abs(opt_egm, 1.0e-12);
                        nlopt_set_initial_step1(opt_egm, 0.00001);

                        
                        if (ik == 0) {
                            x_egm[0] = 0.01;
                        } else {
                            x_egm[0] = solver_data_egm->last_h ;
                        }

                        nlopt_optimize(opt_egm, x_egm, &minf_egm);
                        


                        //sol->h_egm[idx] = x_egm[0];
                        sol->h_egm[idx] = sol->h[idx];
                        //sol->V_egm[idx] = -minf_egm;
                        double [_c, _h] = egm_step(sol->h_egm[idx], sol->c, solver_data, iA);
                        sol->c_egm[idx] = _c;
                        sol->h_egm[idx] = _h;
                        // sol->c_egm[idx] = egm_step(x_egm[0] ,sol->c, solver_data, iA);

                        solver_data_egm->last_h = x_egm[1] / 1.01;
                        */

                    }


                }   

                nlopt_destroy(opt2);
                nlopt_destroy(opt_egm);
        } 
        } // pragma
    }
}