typedef struct par_struct
{
 int T;
 double rho;
 double beta;
 double eta;
 double gamma;
 double alpha;
 double w;
 double tau;
 double r;
 double a_max;
 double a_min;
 int Na;
 double k_max;
 int Nk;
 int simT;
 int simN;
 int threads;
 double* a_grid;
 double* k_grid;
 double* h_sol_grid;
 double* c_sol_grid;
 double* w_vec;
} par_struct;

int get_int_par_struct(par_struct* x, char* name){

 if( strcmp(name,"T") == 0 ){ return x->T; }
 else if( strcmp(name,"Na") == 0 ){ return x->Na; }
 else if( strcmp(name,"Nk") == 0 ){ return x->Nk; }
 else if( strcmp(name,"simT") == 0 ){ return x->simT; }
 else if( strcmp(name,"simN") == 0 ){ return x->simN; }
 else if( strcmp(name,"threads") == 0 ){ return x->threads; }
 else {return -9999;}

}


double get_double_par_struct(par_struct* x, char* name){

 if( strcmp(name,"rho") == 0 ){ return x->rho; }
 else if( strcmp(name,"beta") == 0 ){ return x->beta; }
 else if( strcmp(name,"eta") == 0 ){ return x->eta; }
 else if( strcmp(name,"gamma") == 0 ){ return x->gamma; }
 else if( strcmp(name,"alpha") == 0 ){ return x->alpha; }
 else if( strcmp(name,"w") == 0 ){ return x->w; }
 else if( strcmp(name,"tau") == 0 ){ return x->tau; }
 else if( strcmp(name,"r") == 0 ){ return x->r; }
 else if( strcmp(name,"a_max") == 0 ){ return x->a_max; }
 else if( strcmp(name,"a_min") == 0 ){ return x->a_min; }
 else if( strcmp(name,"k_max") == 0 ){ return x->k_max; }
 else {return NAN;}

}


double* get_double_p_par_struct(par_struct* x, char* name){

 if( strcmp(name,"a_grid") == 0 ){ return x->a_grid; }
 else if( strcmp(name,"k_grid") == 0 ){ return x->k_grid; }
 else if( strcmp(name,"h_sol_grid") == 0 ){ return x->h_sol_grid; }
 else if( strcmp(name,"c_sol_grid") == 0 ){ return x->c_sol_grid; }
 else if( strcmp(name,"w_vec") == 0 ){ return x->w_vec; }
 else {return NULL;}

}


